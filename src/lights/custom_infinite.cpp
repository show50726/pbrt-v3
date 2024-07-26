
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

	// Implementation of Median Cut Algorithm
	CustomInfiniteAreaLight::CustomInfiniteAreaLight(const Transform &LightToWorld,
		const Spectrum &L, int nSamples,
		const std::string &texmap, int maxdepth)
		: Light((int)LightFlags::Infinite, LightToWorld, MediumInterface(),
			nSamples), maxdepth(maxdepth) {
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
		sumAreaTable.reset(new RGBSpectrum[(width + 1) * (height + 1)]);

		CalculateSumAreaTable(texels.get(), width, height, sumAreaTable.get());
		ProcessMedianCut(texels.get(), 0, width - 1, 0, height - 1, maxdepth);
	
		// std::cout << "Light Numbers: " << lightSources.size() << std::endl;

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
		int x, int y) {
		int index = (width + 1) * (y + 1) + (x + 1);
		return sumAreaTable[index];
	}

	RGBSpectrum CustomInfiniteAreaLight::Query(
		int l, int r, int b, int t) {
		return Query(r, t)
			- Query(l - 1, t)
			- Query(r, b - 1)
			+ Query(l - 1, b - 1);
	}

	int CustomInfiniteAreaLight::FindMedianCut(
		int l, int r, int b, int t,
		bool cutVertical) {
		Float totalEnergy = Query(l, r, b, t).y();
		Float targetEnergy = 0.5f * totalEnergy;
		// std::cout << "Target energy: " << targetEnergy << std::endl;

		int min_b = (cutVertical ? b : l) + 1;
		int max_b = cutVertical ? t : r;
		while (min_b < max_b) {
			int mid = min_b + ((max_b - min_b) >> 1);
			Float v = cutVertical ? Query(l, r, b, mid).y() : Query(l, mid, b, t).y();

			if (v < targetEnergy) {
				min_b = mid + 1;
			}
			else {
				max_b = mid;
			}
		}
		return min_b;
	}

	// static 
	void CustomInfiniteAreaLight::CalculateSumAreaTable(
		RGBSpectrum* texmap, 
		int width, int height, 
		RGBSpectrum* sumAreaTable) {
		ParallelFor(
			[sumAreaTable](int64_t i) mutable {
				sumAreaTable[i] = RGBSpectrum();
			},
			width, 32);
		ParallelFor(
			[sumAreaTable, &width](int64_t i) mutable {
				sumAreaTable[i * (width + 1)] = RGBSpectrum();
			},
			height, 32);

		float solidAngle = ((2.f * Pi) / (width - 1)) * ((Pi) / (height - 1));
		// TODO: Make this parallel
		for (int i = 1; i <= width; i++) {
			for (int j = 1; j <= height; j++) {
				int index = j * (width + 1) + i;
				int originalIndex = (j - 1)*width + (i - 1);

				Float sinTheta = std::sin(Pi * (Float)j / (Float)height);

				sumAreaTable[index] = texmap[originalIndex] * solidAngle * sinTheta
					+ sumAreaTable[j * (width + 1) + (i - 1)]
					+ sumAreaTable[(j - 1) * (width + 1) + i]
					- sumAreaTable[(j - 1) * (width + 1) + (i - 1)];
			}
		}
	}

	bool CustomInfiniteAreaLight::ShouldCutVertical(int l, int r, int b, int t) {
		// Float sinTheta = std::sin(Pi * (t + b) * 0.5f / (Float)height);
		// bool cutVertical = (t - b) > ((r - l) * sinTheta);
		return t - b > r - l;
	}

	void CustomInfiniteAreaLight::ProcessMedianCut(
		RGBSpectrum* texmap,
		int l, int r, int b, int t,
		int depth) {
		if (depth == 0) {
			RGBSpectrum color = Query(l, r, b, t);

			CHECK(color.y() >= 0.0f);
			// std::cout << l << " " << r << " " << b << " " << t << std::endl;

			Float u = (Float)(l + r) * 0.5f / (Float)width;
			Float v = (Float)(b + t) * 0.5f / (Float)height;

			//transform coordinate to theta, phi
			lightSources.push_back(LightSource{ color, 2 * Pi * u,  Pi * v });
			return;
		}

		if (r - l <= 1 && t - b <= 1)
			return;
		
		bool cutVertical = ShouldCutVertical(l, r, b, t);
		int cut = FindMedianCut(l, r, b, t, cutVertical);

		if (cut == l || cut == r || cut == b || cut == t) {
			return;
		}

		if (cutVertical) {
			ProcessMedianCut(texmap, l, r, b, cut, depth - 1);
			ProcessMedianCut(texmap, l, r, cut + 1, t, depth - 1);
		}
		else {
			ProcessMedianCut(texmap, l, cut, b, t, depth - 1);
			ProcessMedianCut(texmap, cut + 1, r, b, t, depth - 1);
		}
	}

	Spectrum CustomInfiniteAreaLight::Power() const {
		return Pi * worldRadius * worldRadius *
			Spectrum(Lmap->Lookup(Point2f(.5f, .5f), .5f),
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
		int sampleLightIndex = u[0] * lightSources.size();
		*pdf = 1.0f / lightSources.size();

		// std::cout << "y: " << lightSources[sampleLightIndex].spectrum.y() << std::endl;
		CHECK(lightSources[sampleLightIndex].spectrum.y() >= 0.0f);
		
		Float cosTheta = std::cos(lightSources[sampleLightIndex].theta), sinTheta = std::sin(lightSources[sampleLightIndex].theta);
		Float sinPhi = std::sin(lightSources[sampleLightIndex].phi), cosPhi = std::cos(lightSources[sampleLightIndex].phi);
		*wi =
			LightToWorld(Vector3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta));
		*vis = VisibilityTester(ref, Interaction(ref.p + *wi * (2 * worldRadius),
			ref.time, mediumInterface));
		return Spectrum(lightSources[sampleLightIndex].spectrum,
			SpectrumType::Illuminant);
	}

	Float CustomInfiniteAreaLight::Pdf_Li(const Interaction &, const Vector3f &w) const {
		ProfilePhase _(Prof::LightPdf);
		return 1.0f / lightSources.size();
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
		int maxDepth = paramSet.FindOneInt("maxdepth", 6);
		if (PbrtOptions.quickRender) nSamples = std::max(1, nSamples / 4);
		return std::make_shared<CustomInfiniteAreaLight>(light2world, L * sc, nSamples,
			texmap, maxDepth);
	}

}  // namespace pbrt
