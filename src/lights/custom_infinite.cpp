
/*
pbrt source code is Copyright(c) 1998-2016
Matt Pharr, Greg Humphreys, and Wenzel Jakob.

This file is part of pbrt.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "lights/custom_infinite.h"

#include <array>
#include <queue>

#include "imageio.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

	// CustomInfiniteAreaLight Method Definitions
	// Implementation of Median Cut Algorithm
	CustomInfiniteAreaLight::CustomInfiniteAreaLight(const Transform &LightToWorld,
		const Spectrum &L, int nSamples,
		const std::string &texmap, int maxregionnum)
		: Light((int)LightFlags::Infinite, LightToWorld, MediumInterface(),
			nSamples), maxRegionNum(maxregionnum) {
		// Read texel data from _texmap_ and initialize _Lmap_
		Point2i resolution;
		std::unique_ptr<RGBSpectrum[]> texels(nullptr);
		if (texmap != "") {
			texels = ReadImage(texmap, &resolution);
			if (texels)
				for (int i = 0; i < resolution.x * resolution.y; ++i)
					texels[i] *= L.ToRGBSpectrum();
		}
		if (!texels) {
			resolution.x = resolution.y = 1;
			texels = std::unique_ptr<RGBSpectrum[]>(new RGBSpectrum[1]);
			texels[0] = L.ToRGBSpectrum();
		}
		Lmap.reset(new MIPMap<RGBSpectrum>(resolution, texels.get()));

		width = resolution.x;
		height = resolution.y;
		sumAreaTable.reset(new RGBSpectrum[(width+1) * (height+1)]);
		lightSources.reset(new LightSource[maxRegionNum]);

		CalculateSumAreaTable(texels.get(), width, height, sumAreaTable.get());
		ProcessMedianCut(texels.get(), sumAreaTable.get(), width, height, lightSources.get());
	

		std::unique_ptr<Float[]> img(new Float[width * height]);
		float fwidth = 0.5f / std::min(width, height);
		ParallelFor(
			[&](int64_t v) {
				Float vp = (v + .5f) / (Float)height;
				Float sinTheta = std::sin(Pi * (v + .5f) / height);
				for (int u = 0; u < width; ++u) {
					Float up = (u + .5f) / (Float)width;
					img[u + v * width] = Lmap->Lookup(Point2f(up, vp), fwidth).y();
					img[u + v * width] *= sinTheta;
				}
			},
			height, 32);
		distribution.reset(new Distribution2D(img.get(), width, height));
	}
	
	RGBSpectrum CustomInfiniteAreaLight::Query(
		int r, int t) {
		int index = height * t + r;
		return sumAreaTable[index];
	}

	//static 
	RGBSpectrum CustomInfiniteAreaLight::Query(
		int l, int r, int b, int t) {
		// Note: l and b are not included
		RGBSpectrum b0 = Query(r, t);
		RGBSpectrum b1 = Query(l, t);
		RGBSpectrum b2 = Query(r, b);
		RGBSpectrum b3 = Query(l, b);

		return b0 - b1 - b2 + b3;
	}

	int CustomInfiniteAreaLight::FindMedianCut(
		int l, int r, int b, int t) {
		bool cutVertical = (t - b) > (r - l);

		Float totalEnergy = Query(l, r, b, t).y();
		Float targetEnergy = 0.5f * totalEnergy;
	
		int start = cutVertical ? b : l;
		int end = cutVertical ? t : r;
		while (start < end) {
			int mid = start + ((end - start) >> 1);
			Float v = cutVertical ? Query(l, r, b, mid).y() : Query(l, mid, b, t).y();

			if (v < targetEnergy) {
				start = mid + 1;
			}
			else {
				end = mid;
			}
		}
		return start;
	}

	// static 
	void CustomInfiniteAreaLight::CalculateSumAreaTable(
		RGBSpectrum* texmap, 
		int width, int height, 
		RGBSpectrum* sumAreaTable) {
		sumAreaTable[0] = texmap[0].y();
		std::cout << "Start Init Table" << std::endl;
		ParallelFor(
			[sumAreaTable](int64_t i) mutable {
				sumAreaTable[i] = RGBSpectrum();
			},
			width, 32);
		ParallelFor(
			[sumAreaTable, &width](int64_t i) mutable {
				sumAreaTable[i * width] = RGBSpectrum();
			},
			height, 32);
		std::cout << "Finish Init Table" << std::endl;

		std::cout << "Start Fill in Table" << std::endl;
		// TODO: Make this parallel
		for (int i = 1; i <= width; i++) {
			for (int j = 1; j <= height; j++) {
				int index = j * width + i;
				int originalIndex = (j - 1)*width + (i - 1);
				sumAreaTable[index] = 
					sumAreaTable[(j - 1)*width + i]
					+ sumAreaTable[j * width + (i - 1)]
					+ texmap[originalIndex]
					- sumAreaTable[originalIndex];
			}
		}
		std::cout << "Finish Fill in Table" << std::endl;
	}

	void CustomInfiniteAreaLight::ProcessMedianCut(
		RGBSpectrum* texmap,
		RGBSpectrum* sumAreaTable,
		int width, int height,
		LightSource* lightSources) {
		std::queue<std::array<int, 4>> blocks;
		blocks.push(std::array<int, 4>{ 0, width, 0, height });
		
		while (blocks.size() < maxRegionNum) {
			int l = blocks.front()[0];
			int r = blocks.front()[1];
			int b = blocks.front()[2];
			int t = blocks.front()[3];
			blocks.pop();

			if (r - l < 2 && t - b < 2) {
				blocks.push(std::array<int, 4>{ l, r, b, t });
				continue;
			}

			int cut = FindMedianCut(l, r, b, t);
			if (r - l < t - b) {
				blocks.push(std::array<int, 4>{ l, r, cut, t });
				blocks.push(std::array<int, 4>{ l, r, b, cut });
			}
			else {
				blocks.push(std::array<int, 4>{ l, cut, b, t });
				blocks.push(std::array<int, 4>{ cut, r, b, t });
			}
		}

		std::cout << "Blocks Size: " << blocks.size() << std::endl;

		for (int i = 0; i < blocks.size(); i++) {
			std::array<int, 4> current = blocks.front();
			blocks.pop();

			int num = (current[1] - current[0])*(current[3] - current[2]);
			RGBSpectrum color = Query(current[0], current[1], current[2], current[3]);

			if (color.y() < 0.0f) {
				std::cout << current[0] << " " << current[1] << " " << current[2] << " " << current[3] << std::endl;
				//std::cout << b0.y() << " " << b1.y() << " " << b2.y() << " " << b3.y() << std::endl;
			}
			std::cout << num << std::endl;
			//color /= (float)num;
			CHECK(color.y() >= 0.0f);

			Float u = ((Float)current[1] - 0.5f) / (Float)width;
			Float v = ((Float)current[3] - 0.5f) / (Float)height;

			//transform coordinate to theta, phi
			lightSources[i] = LightSource{ color, 2 * Pi * u,  Pi * v };
		}
	}

	Spectrum CustomInfiniteAreaLight::Power() const {
		return Pi * worldRadius * worldRadius *
			Spectrum(sumAreaTable[width + height * width],
				SpectrumType::Illuminant);
	}

	Spectrum CustomInfiniteAreaLight::Le(const RayDifferential &ray) const {
		Vector3f w = Normalize(WorldToLight(ray.d));
		Point2f st(SphericalPhi(w) * Inv2Pi, SphericalTheta(w) * InvPi);
		return Spectrum(Lmap->Lookup(st), SpectrumType::Illuminant);
	}

	Spectrum CustomInfiniteAreaLight::Sample_Li(const Interaction &ref, const Point2f &u,
		Vector3f *wi, Float *pdf,
		VisibilityTester *vis) const {
		ProfilePhase _(Prof::LightSample);
		
		// Pick which light source to sample
		CHECK(u[0] < 1 && u[0] >= 0);
		int sampleLightIndex = u[0] * maxRegionNum;

		return Spectrum(lightSources[sampleLightIndex].spectrum, 
			SpectrumType::Illuminant);
	}

	Float CustomInfiniteAreaLight::Pdf_Li(const Interaction &, const Vector3f &w) const {
		ProfilePhase _(Prof::LightPdf);
		return 1.0f / maxRegionNum;
	}

	Spectrum CustomInfiniteAreaLight::Sample_Le(const Point2f &u1, const Point2f &u2,
		Float time, Ray *ray, Normal3f *nLight,
		Float *pdfPos, Float *pdfDir) const {
		ProfilePhase _(Prof::LightSample);
		// Compute direction for infinite light sample ray
		Point2f u = u1;

		// Find $(u,v)$ sample coordinates in infinite light texture
		Float mapPdf;
		Point2f uv = distribution->SampleContinuous(u, &mapPdf);
		if (mapPdf == 0) return Spectrum(0.f);
		Float theta = uv[1] * Pi, phi = uv[0] * 2.f * Pi;
		Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
		Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
		Vector3f d =
			-LightToWorld(Vector3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta));
		*nLight = (Normal3f)d;

		// Compute origin for infinite light sample ray
		Vector3f v1, v2;
		CoordinateSystem(-d, &v1, &v2);
		Point2f cd = ConcentricSampleDisk(u2);
		Point3f pDisk = worldCenter + worldRadius * (cd.x * v1 + cd.y * v2);
		*ray = Ray(pDisk + worldRadius * -d, d, Infinity, time);

		// Compute _InfiniteAreaLight_ ray PDFs
		*pdfDir = sinTheta == 0 ? 0 : mapPdf / (2 * Pi * Pi * sinTheta);
		*pdfPos = 1 / (Pi * worldRadius * worldRadius);
		return Spectrum(Lmap->Lookup(uv), SpectrumType::Illuminant);
	}

	void CustomInfiniteAreaLight::Pdf_Le(const Ray &ray, const Normal3f &, Float *pdfPos,
		Float *pdfDir) const {
		ProfilePhase _(Prof::LightPdf);
		Vector3f d = -WorldToLight(ray.d);
		Float theta = SphericalTheta(d), phi = SphericalPhi(d);
		Point2f uv(phi * Inv2Pi, theta * InvPi);
		Float mapPdf = distribution->Pdf(uv);
		*pdfDir = mapPdf / (2 * Pi * Pi * std::sin(theta));
		*pdfPos = 1 / (Pi * worldRadius * worldRadius);
	}

	std::shared_ptr<CustomInfiniteAreaLight> CreateCustomInfiniteLight(
		const Transform &light2world, const ParamSet &paramSet) {
		Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
		Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
		std::string texmap = paramSet.FindOneFilename("mapname", "");
		int nSamples = paramSet.FindOneInt("samples",
			paramSet.FindOneInt("nsamples", 1));
		int maxRegionNum = paramSet.FindOneInt("maxregion", 64);
		if (PbrtOptions.quickRender) nSamples = std::max(1, nSamples / 4);
		return std::make_shared<CustomInfiniteAreaLight>(light2world, L * sc, nSamples,
			texmap, maxRegionNum);
	}

}  // namespace pbrt
