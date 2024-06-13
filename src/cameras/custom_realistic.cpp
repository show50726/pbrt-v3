
#include "cameras/custom_realistic.h"

#include "core/parser.h"
#include "core/sampling.h"
#include "floatfile.h"
#include "reflection.h"

namespace pbrt {
namespace {

static bool _IntersectLensTest(
	/* camera space */ const Ray& ray,
	Float curvatureCenterZ,
	Float curvatureRadius,
	Float apertureRadius,
	Float* t,
	Normal3f* normal) {
	// Transform the ray into lens space
	Point3f newRayOrigin = ray.o - Vector3f(0, 0, curvatureCenterZ);
	Float A = ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z;
	Float B = 2 * (ray.d.x * newRayOrigin.x + ray.d.y * newRayOrigin.y + ray.d.z * newRayOrigin.z);
	Float C = newRayOrigin.x * newRayOrigin.x
		+ newRayOrigin.y * newRayOrigin.y
		+ newRayOrigin.z * newRayOrigin.z
		- curvatureRadius * curvatureRadius;
	Float t0, t1;
	if (!Quadratic(A, B, C, &t0, &t1))
		return false;

	bool useCloserT = (ray.d.z > 0) ^ (curvatureRadius < 0);
	*t = useCloserT ? std::min(t0, t1) : std::max(t0, t1);
	if (*t < 0) 
		return false;

	*normal = Normal3f(Vector3f(newRayOrigin + *t * ray.d));
	*normal = Faceforward(Normalize(*normal), -ray.d);
	return true;
}

static Vector3f _CalculateRefractedRay(
	Vector3f incidentDirection,
	Normal3f normal,
	Float n1, Float n2) {
	// Heckbert's method
	CHECK_NE(n2, 0);
	Float n = n1 / n2;
	Float c1 = Dot(incidentDirection, normal);
	Float c2 = std::sqrt(1.0f - n * n*(1.0f - c1 * c1));
	return n * -incidentDirection + (Vector3f)normal * (n * c1 - c2);
}

} // namespace

std::vector<CustomRealisticCamera::LensElementInterface> CustomRealisticCamera::_ReadLensFromFile(std::string filename) {
	std::vector<Float> lensData;
	if (!ReadFloatFile(filename.c_str(), &lensData)) {
		Error("Error reading lens specification file \"%s\".",
			filename.c_str());
		return {};
	}

	if (lensData.size() % 4 != 0) {
		Error(
			"Excess values in lens specification file \"%s\"; "
			"must be multiple-of-four values, read %d.",
			filename.c_str(), (int)lensData.size());
		return {};
	}

	std::vector<LensElementInterface> result;
	for (int i = 0; i < lensData.size(); i += 4) {
		// The lens units in dat files are measured in mm (millimeters), but in pbrt the unit is m (meters). 
		LensElementInterface lensElement;
		lensElement.curvatureRadius = 0.001f * lensData[i];
		lensElement.axPos = 0.001f * lensData[i + 1];
		lensElement.nd = lensData[i + 2];
		lensElement.apertureRadius = 0.001f * 0.5f * lensData[i + 3];
		result.push_back(lensElement);
	}

	return result;
}

CustomRealisticCamera::CustomRealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, std::string specfile, 
				 float filmdiag, Film *f, const Medium *medium)
	: Camera(cam2world, sopen, sclose, f, medium), shutter_opon_(sopen), shutter_close_(sclose),
	film_distance_(filmdistance * 0.001f), aperture_radius_(aperture_diameter * 0.001f * 0.5f), film_diag_(filmdiag * 0.001f),
	lens_system_(std::move(_ReadLensFromFile(specfile)))
{
}

Float CustomRealisticCamera::_CalculateTotalThickness() const {
	Float total = 0;
	for (int i = 0; i < lens_system_.size(); i++) {
		total += lens_system_[i].axPos;
	}
	return total;
}


bool CustomRealisticCamera::_CastRayFromFilm(
	/* camera space */ const Ray& input, 
	/* camera space */ Ray* output) const {
	static const Transform CameraToLens = Scale(1, 1, -1);
	Float currentZ = -film_distance_;
	Ray currentRay = CameraToLens(input);
	for (int i = lens_system_.size() - 1; i >= 0; i--) {
		const LensElementInterface& lens = lens_system_[i];

		currentZ -= lens.axPos;

		bool isStop = (lens.curvatureRadius == 0);
		Float t;
		Normal3f normal;
		if (isStop) {
			// Aperture Stop.
			// Propergate the ray without refraction
			if (currentRay.d.z >= 0.0f)
				return false;
			t = (currentZ - currentRay.o.z) / currentRay.d.z;
		}
		else {
			if (!_IntersectLensTest(
				currentRay,
				currentZ + lens.curvatureRadius,
				lens.curvatureRadius,
				lens.apertureRadius,
				&t,
				&normal)) {
				return false;
			}
		}
		
		if (isnan(t))
			return false;

		Point3f pHit = currentRay(t);
		Float r2 = pHit.x * pHit.x + pHit.y * pHit.y;
		if (r2 > lens_system_[i].apertureRadius * lens_system_[i].apertureRadius) 
			return false;
		currentRay.o = pHit;

		if (!isStop)
		{
			Float n1 = lens_system_[i].nd;
			Float n2 = (i > 0 && lens_system_[i - 1].nd != 0)
				? lens_system_[i - 1].nd
				: 1;
			Vector3f refracted = _CalculateRefractedRay(Normalize(-currentRay.d), normal, n1, n2);

			currentRay.d = refracted;
		}

	}

	if (output != nullptr) {
		static const Transform LensToCamera = Scale(1, 1, -1);
		*output = LensToCamera(currentRay);
	}
	return true;
}

Bounds2f GetPhysicalExtent(const Film& film, Float filmdiag) {
	Float aspect = (Float)film.fullResolution.y / (Float)film.fullResolution.x;
	Float x = std::sqrt(filmdiag * filmdiag / (1 + aspect * aspect));
	Float y = aspect * x;
	return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}

float CustomRealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	ProfilePhase prof(Prof::GenerateCameraRay);

	// 1. Sample from exit pupil
	// 2. Cast the ray, do the full simulation
	// 3. If the ray is blocked, fill the weight with 0, 
	//    or fill the cos term to the weight

	Float totalThickness = _CalculateTotalThickness();
	// Convert to sample to [0, 1]
	Point2f s(sample.pFilm.x / film->fullResolution.x,
		sample.pFilm.y / film->fullResolution.y);
	// Map the [0, 1] from screen to film
	Point2f pFilm = GetPhysicalExtent(*film, film_diag_).Lerp(s);
	//std::cout << "s: " << sample.pFilm << std::endl;
	Point3f filmPoint = Point3f(-pFilm.x, pFilm.y, 0);
	//std::cout << filmPoint * 1000 << std::endl;

	// Sample from the exit pupil
	// Use first lens as exit pupil
	Float exitPupilRadius = lens_system_.back().apertureRadius;
	Point2f sampleOnDisk = exitPupilRadius * ConcentricSampleDisk(sample.pLens);
	Point3f exitPupilPoint = Point3f(sampleOnDisk.x, sampleOnDisk.y, film_distance_);
	
	// Camera space ray
	Ray currentRay(filmPoint, exitPupilPoint - filmPoint, Infinity, 1.0f);
	// std::cout << currentRay.d.x << " " << currentRay.d.y << " " << currentRay.d.z << std::endl;
	if (!_CastRayFromFilm(currentRay, ray)) {
		return 0;
	}

	ray->o -= Vector3f(0, 0, totalThickness + film_distance_);
	*ray = CameraToWorld(*ray);
	ray->d = Normalize(ray->d);
	ray->medium = medium;

	// Calculate the weight of the ray
	Float cosTheta = Normalize(currentRay.d).z;
	Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
	Float sampleArea = Pi * exitPupilRadius * exitPupilRadius;

	Float a2 = filmPoint.x * filmPoint.x + filmPoint.y * filmPoint.y;
	Float z2 = film_distance_ * film_distance_;
	Float r2 = exitPupilRadius * exitPupilRadius;
	Float x1 = a2 + z2 - r2;
	Float x2 = std::sqrt((a2 + z2 + r2)*(a2 + z2 + r2) - 4.0f*r2*a2);
	Float f = 0.5f * (1.0f - x1 / x2);

	//return f;
	return cos4Theta * sampleArea / (film_distance_ * film_distance_);
}


CustomRealisticCamera *CreateCutomRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film, const Medium *medium) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	std::string lensFile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (lensFile == "") {
	    Error("No lens spec file supplied!\n");
	}
	return new CustomRealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   lensFile, filmdiag, film, medium);
}

}