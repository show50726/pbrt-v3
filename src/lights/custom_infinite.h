
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_LIGHTS_CUSTOM_INFINITE_H
#define PBRT_LIGHTS_CUSTOM_INFINITE_H

// lights/infinite.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

namespace pbrt {

	// InfiniteAreaLight Declarations
	class CustomInfiniteAreaLight : public Light {
	public:
		// InfiniteAreaLight Public Methods
		CustomInfiniteAreaLight(const Transform &LightToWorld, const Spectrum &power,
			int nSamples, const std::string &texmap, int maxregionnum);
		
		void Preprocess(const Scene &scene) {
			scene.WorldBound().BoundingSphere(&worldCenter, &worldRadius);
		}
		
		Spectrum Power() const;
		Spectrum Le(const RayDifferential &ray) const;
		Spectrum Sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wi,
			Float *pdf, VisibilityTester *vis) const;
		Float Pdf_Li(const Interaction &, const Vector3f &) const;
		Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
			Ray *ray, Normal3f *nLight, Float *pdfPos,
			Float *pdfDir) const;
		void Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos,
			Float *pdfDir) const;

	private:
		struct LightSource {
			RGBSpectrum spectrum;
			Float theta;
			Float phi;
		};

		static void CalculateSumAreaTable(RGBSpectrum* texmap, int width, int height, RGBSpectrum* sumAreaTable);
		RGBSpectrum Query(int r, int t);
		RGBSpectrum Query(int l, int r, int b, int t);
		int FindMedianCut(int l, int r, int b, int t);
		void ProcessMedianCut(
			RGBSpectrum* texmap, 
			int l, int r, int b, int t,
			int depth);

		std::unique_ptr<MIPMap<RGBSpectrum>> Lmap;
		std::unique_ptr<RGBSpectrum[]> sumAreaTable;
		std::vector<LightSource> lightSources;
		std::unique_ptr<Distribution2D> distribution;
		int maxdepth;
		Point3f worldCenter;
		Float worldRadius;
		int width;
		int height;
	};

	std::shared_ptr<CustomInfiniteAreaLight> CreateCustomInfiniteLight(
		const Transform &light2world, const ParamSet &paramSet);

}  // namespace pbrt

#endif  // PBRT_LIGHTS_CUSTOM_INFINITE_H
