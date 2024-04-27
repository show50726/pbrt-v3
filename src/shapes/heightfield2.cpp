
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

#include <assert.h>

// shapes/heightfield.cpp*
#include "paramset.h"
#include "primitive.h"
#include "sampling.h"
#include "shapes/heightfield2.h"
#include "stats.h"
#include "error.h"
#include "triangle.h"


namespace pbrt {

bool Heightfield::Voxel::Intersect(
	const Ray &ray, 
	SurfaceInteraction *isect) {
	Float minT = INT_FAST32_MAX;
	for (size_t i = 0; i < triangles.size(); ++i) {
		SurfaceInteraction surf;
		if (!Intersect(
			ray, 
			triangles[i], 
			&surf))
			continue;

		// TODO: Use primitive with smallest T
		*isect = surf;
		return true;
	}
	return false;
}

bool Heightfield::Voxel::Intersect(
	const Ray &ray, 
	std::vector<Point3f> triangle, 
	SurfaceInteraction *isect) {
	const Point3f &p1 = triangle[0];
	const Point3f &p2 = triangle[1];
	const Point3f &p3 = triangle[2];
	Vector3f e1 = p2 - p1;
	Vector3f e2 = p3 - p1;
	Vector3f s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);

	if (divisor == 0.)
		return false;
	float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	Vector3f s = ray.o - p1;
	float b1 = Dot(s, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;

	// Compute second barycentric coordinate
	Vector3f s2 = Cross(s, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
		return false;

	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < 0 || t > ray.tMax)
		return false;

	Vector3f dpdu, dpdv;
	CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
	*isect = (*objToWorld)(SurfaceInteraction(ray(t), Vector3f(0, 0, 0), Point2f(0, 0),
		-ray.d, dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0),
		ray.time, nullptr));
	return true;
}

bool Heightfield::Voxel::IntersectP(
	const Ray &ray, 
	std::vector<Point3f> triangle) {
	return true;
}

bool Heightfield::Voxel::IntersectP(const Ray &ray) {
	for (size_t i = 0; i < triangles.size(); ++i) {
		SurfaceInteraction surf;
		if (!IntersectP(
			ray,
			triangles[i]))
			continue;

		return true;
	}
	return false;
}

Heightfield::GridAccel::GridAccel(Heightfield* hf) 
	: heightfield(hf) {
	nVoxels[0] = hf->nx - 1;
	nVoxels[1] = hf->ny - 1;
	nVoxels[2] = 1;
	bounds = hf->ObjectBound();
	width = invWidth = Vector3f(1, 1, 1);

	int nv = nVoxels[0] * nVoxels[1] * nVoxels[2];
	voxels = AllocAligned<Voxel *>(nv);
	// TODO: need memset?
	memset(voxels, 0, nv * sizeof(Voxel *));

	int totalT = 0;
	for (int i = 0; i < nVoxels[0]; i++) {
		for (int j = 0; j < nVoxels[1]; j++) {
			int index = offset(i, j, 0);
			Point3f bottomLeft = Point3f(
				1.0f / nVoxels[0] * i,
				1.0f / nVoxels[1] * j,
				heightfield->z[index]);
			int bottomRightIndex = offset(i + 1, j, 0);
			Point3f bottomRight = Point3f(
				1.0f / nVoxels[0] * (i + 1),
				1.0f / nVoxels[1] * j,
				heightfield->z[bottomRightIndex]);
			int topLeftIndex = offset(i, j + 1, 0);
			Point3f topLeft = Point3f(
				1.0f / nVoxels[0] * i,
				1.0f / nVoxels[1] * (j + 1),
				heightfield->z[topLeftIndex]);
			int topRightIndex = offset(i + 1, j + 1, 0);
			Point3f topRight = Point3f(
				1.0f / nVoxels[0] * (i + 1),
				1.0f / nVoxels[1] * (j + 1),
				heightfield->z[topRightIndex]);
			std::vector<Point3f> triangle1 =
			{
				bottomLeft,
				topLeft,
				bottomRight
			};
			std::vector<Point3f> triangle2 =
			{
				bottomLeft,
				topRight,
				bottomRight
			};
			totalT += 2;
			// Make two triangles and add them into voxels
			/*auto leftTriangle = std::make_shared<Triangle>(
				hf->ObjectToWorld, hf->WorldToObject,
				false, mesh, i)*/
			if (!voxels[index]) {
				// Allocate new voxel and store primitive in it
				voxels[index] = voxelArena.Alloc<Voxel>();
				*voxels[index] = Voxel(heightfield->ObjectToWorld, std::move(triangle1));
				(*voxels[index]).AddTriangle(std::move(triangle2));
			}
			else {
				// Add primitive to already-allocated voxel
				voxels[index]->AddTriangle(std::move(triangle1));
				voxels[index]->AddTriangle(std::move(triangle2));
			}
		}
	}
	std::cout << totalT << std::endl;
}

Heightfield::GridAccel::~GridAccel() {
	for (int i = 0; i < nVoxels[0] * nVoxels[1] * nVoxels[2]; ++i)
		if (voxels[i]) voxels[i]->~Voxel();
	FreeAligned(voxels);
}

bool Heightfield::GridAccel::Intersect(
	const Ray &ray, 
	SurfaceInteraction *isect) const {
	int n = nVoxels[0] * nVoxels[1] * nVoxels[2];
	SurfaceInteraction currSect;
	for (int i = 0; i < n; i++) {
		if (!voxels[i]->Intersect(
			ray, &currSect))
			continue;
		(*isect) = currSect;
		return true;
	}
	return false;
}

bool Heightfield::GridAccel::IntersectP(const Ray &ray) const {
	// TODO: Impl
	return false;
}

Heightfield::Heightfield(const Transform *ObjectToWorld,
                         const Transform *WorldToObject,
                         bool reverseOrientation, int nu, int nv,
                         const float *zs)
    : Shape(ObjectToWorld, WorldToObject, reverseOrientation), nx(nu), ny(nv) {
    z = new float[nx * ny];
    memcpy(z, zs, nx * ny * sizeof(float));
	grid = std::make_unique<GridAccel>(this);
}

Heightfield::~Heightfield() {
    delete[] z;
}

Bounds3f Heightfield::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx * ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return Bounds3f(Point3f(0, 0, minz), Point3f(1, 1, maxz));
}

bool Heightfield::Intersect(const Ray &ray, Float *tHit,
                            SurfaceInteraction *isect,
                            bool testAlphaTexture) const {
	Vector3f oErr, dErr;
	Ray rayOS = (*WorldToObject)(ray, &oErr, &dErr);
	
	Float hitt0, hitt1;
    bool hitBound = ObjectBound().IntersectP(rayOS, &hitt0, &hitt1);
    if (!hitBound) return false;

	bool hitGrid = grid->Intersect(rayOS, isect);
	if (!hitGrid) return false;
	return true;

	/*tHit = &hitt0;
	Point3f pHit = rayOS(*tHit);
	Vector3f pError = Vector3f(0, 0, 0);
	*isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(0, 0),
		-rayOS.d, Vector3f(1, 0, 0), Vector3f(0, 0, 1), Normal3f(0, 0, 0), Normal3f(0, 0, 0),
		rayOS.time, this));*/
	/*return true;*/
}

bool Heightfield::IntersectP(const Ray &ray, bool testAlphaTexture) const {
	/*float rayStartZ = ray.o.z;
	if (fabsf(ray.d.z) < 1e-7) return false;
	float t = ray.d.z == 0
		? -1.0f
		: (0.0f - rayStartZ) / ray.d.z;
	if (t < 0.0f || t > ray.tMax)
		return false;*/
	// TODO: Impl
	return false;
}

Float Heightfield::Area() const { return 1.0f; }

Interaction Heightfield::Sample(const Point2f &u, Float *pdf) const {
    Error("Not implemented 1");
    return Interaction();
}

std::shared_ptr<Heightfield> CreateHeightfield2(
    const Transform *o2w, const Transform *w2o,
                                bool reverseOrientation,
                                const ParamSet &params){
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    assert(nitems == nu * nv);
    assert(nu != -1 && nv != -1 && Pz != NULL);
    return std::make_shared<Heightfield>(
		o2w, 
		w2o, 
		reverseOrientation, 
		nu, nv, Pz);
}

}  // namespace pbrt
