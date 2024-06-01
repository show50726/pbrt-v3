
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CAMERAS_CUSTOM_REALISTIC_H
#define PBRT_CAMERAS_CUSTOM_REALISTIC_H


// cameras/custom_realistic.h*
#include "pbrt.h"
#include "camera.h"
#include "paramset.h"
#include "film.h"

namespace pbrt {

// RealisticCamera Declarations
class CustomRealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	CustomRealisticCamera(const AnimatedTransform &cam2world,
		float hither, float yon, float sopen,
		float sclose, float filmdistance, float aperture_diameter, std::string specfile,
		float filmdiag, Film *film, const Medium *medium);
	float GenerateRay(const CameraSample &sample, Ray *) const;

private:
	// RealisticCamera Public Methods

};


CustomRealisticCamera *CreateCutomRealisticCamera(const ParamSet &params,
	const AnimatedTransform &cam2world, Film *film, const Medium *medium);

}

#endif	// PBRT_CAMERAS_CUSTOM_REALISTIC_H