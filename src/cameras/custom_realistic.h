
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
		float hither, float yon, float sopen, float sclose, 
		// filmdistance: the distance between the film and the lens closet to the film
		float filmdistance, 
		float aperture_diameter, std::string specfile,
		// filmdiag: distance from the lower left of the film image to the upper right
		float filmdiag, 
		Film *film, const Medium *medium);
	float GenerateRay(const CameraSample &sample, Ray *) const;

private:
	// RealisticCamera Public Methods
	struct LensElementInterface {
		Float curvatureRadius;
		// axPos: Distance from the current surface to next surface
		Float axPos;
		// nd: index of refraction
		Float nd;
		Float apertureRadius;
	};

	std::vector<LensElementInterface> _ReadLensFromFile(std::string filename);
	
	// Do the fully simulation to trace the ray throught the lens system
	// TODO: Implement thick lens approximation
	bool _CastRayFromFilm(
		/* camera space */ const Ray& input, 
		/* camera space */ Ray* output) const;

	Float _CalculateTotalThickness() const;
	
	std::vector<LensElementInterface> lens_system_;
	Float shutter_opon_;
	Float shutter_close_;
	Float film_distance_;
	Float aperture_radius_;
	Float film_diag_;
};


CustomRealisticCamera *CreateCutomRealisticCamera(const ParamSet &params,
	const AnimatedTransform &cam2world, Film *film, const Medium *medium);

}

#endif	// PBRT_CAMERAS_CUSTOM_REALISTIC_H