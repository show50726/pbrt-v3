
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
	Float* tHit,
	SurfaceInteraction *isect) {
	for (size_t i = 0; i < triangles.size(); ++i) {
		if (!Intersect(
			ray, 
			tHit,
			triangles[i], 
			normals[i],
			isect))
			continue;

		return true;
	}
	return false;
}

bool Heightfield::Voxel::Intersect(
	const Ray &ray, 
	Float* tHit,
	std::vector<Point3f> triangle, 
	std::vector<Normal3f> normal,
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

	float b0 = 1.0f - b1 - b2;

	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < 0 || t > ray.tMax || t > *tHit)
		return false;

	Normal3f interpolatedNormal = Normalize(normal[0] * b0 + normal[1] * b1 + normal[2] * b2);

	// Calculate tangent and bitangent
	Vector3f dpdu, dpdv;
	Point3f p = ray(t);
	CoordinateSystem(Vector3f(interpolatedNormal), &dpdu, &dpdv);

	// Calculate dndu, dndv
	Normal3f dndu, dndv;
	float du1 = p1.x - p3.x;
	float du2 = p2.x - p3.x;
	float dv1 = p1.y - p3.y;
	float dv2 = p2.y - p3.y;
	Normal3f dn1 = normal[0] - normal[2];
	Normal3f dn2 = normal[1] - normal[2];
	float determinant = du1 * dv2 - dv1 * du2;
	if (determinant == 0.f)
		dndu = dndv = Normal3f(0, 0, 0);
	else {
		float invdet = 1.f / determinant;
		dndu = (dv2 * dn1 - dv1 * dn2) * invdet;
		dndv = (-du2 * dn1 + du1 * dn2) * invdet;
	}

	// Calculate error bounds
	Float xAbsSum =
		(std::abs(b0 * p1.x) + std::abs(b1 * p2.x) + std::abs(b2 * p3.x));
	Float yAbsSum =
		(std::abs(b0 * p1.y) + std::abs(b1 * p2.y) + std::abs(b2 * p3.y));
	Float zAbsSum =
		(std::abs(b0 * p1.z) + std::abs(b1 * p2.z) + std::abs(b2 * p3.z));
	Vector3f pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);

	*tHit = t;
	*isect = (*objToWorld)(SurfaceInteraction(p, pError, Point2f(p.x, p.y),
		-ray.d, dpdu, dpdv, dndu, dndv,
		ray.time, nullptr));

	
	interpolatedNormal = Normalize((*objToWorld)(interpolatedNormal));
	isect->n = isect->shading.n = interpolatedNormal;
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

	Vector3f delta = bounds.pMax - bounds.pMin;
	for (int axis = 0; axis < 3; ++axis) {
		width[axis] = delta[axis] / nVoxels[axis];
		invWidth[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
	}

	int vertexNum = hf->nx * hf->ny;
	vertexPositions = AllocAligned<Point3f>(vertexNum);
	normals = AllocAligned<Normal3f>(vertexNum);
	memset(vertexPositions, 0, vertexNum * sizeof(Point3f));
	memset(normals, 0, vertexNum * sizeof(Normal3f));
	float invNVoxels[2] = { 1.0f / nVoxels[0],1.0f / nVoxels[1] };
	// Fill in vertex positions
	for (int i = 0; i < hf->nx; i++) {
		for (int j = 0; j < hf->ny; j++) {
			int vertexIndex = vertexOffset(i, j);
			vertexPositions[vertexIndex] = Point3f(
				invNVoxels[0] * i,
				invNVoxels[1] * j,
				heightfield->z[vertexIndex]);
		}
	}

	// Fill in normals
	for (int i = 0; i < hf->nx - 1; i++) {
		for (int j = 0; j < hf->ny - 1; j++) {
			int bottomLeftIndex = vertexOffset(i, j);
			Point3f bottomLeft = vertexPositions[bottomLeftIndex];
			int bottomRightIndex = vertexOffset(i + 1, j);
			Point3f bottomRight = vertexPositions[bottomRightIndex];
			int topLeftIndex = vertexOffset(i, j + 1);
			Point3f topLeft = vertexPositions[topLeftIndex];
			int topRightIndex = vertexOffset(i + 1, j + 1);
			Point3f topRight = vertexPositions[topRightIndex];

			Normal3f leftNormal = Normal3f(Normalize(Cross(topRight - bottomLeft, topLeft - bottomLeft)));
			Normal3f rightNormal = Normal3f(Normalize(Cross(bottomRight - bottomLeft, topRight - bottomLeft)));
			normals[bottomLeftIndex] += leftNormal;
			normals[bottomLeftIndex] += rightNormal;
			normals[bottomRightIndex] += rightNormal;
			normals[topLeftIndex] += leftNormal;
			normals[topRightIndex] += leftNormal;
			normals[topRightIndex] += rightNormal;
		}
	}


	for (int i = 0; i < vertexNum; i++) {
		normals[i] = Normalize(normals[i]);
	}

	int nv = nVoxels[0] * nVoxels[1] * nVoxels[2];
	voxels = AllocAligned<Voxel *>(nv);
	memset(voxels, 0, nv * sizeof(Voxel *));
	// Put triangles into each cell
	for (int i = 0; i < nVoxels[0]; i++) {
		for (int j = 0; j < nVoxels[1]; j++) {
			int bottomLeftIndex = vertexOffset(i, j);
			int bottomRightIndex = vertexOffset(i + 1, j);
			int topLeftIndex = vertexOffset(i, j + 1);
			int topRightIndex = vertexOffset(i + 1, j + 1);

			std::vector<Point3f> triangle1 =
			{
				vertexPositions[bottomLeftIndex],
				vertexPositions[topLeftIndex],
				vertexPositions[topRightIndex]
			};
			std::vector<Normal3f> normal1 =
			{
				normals[bottomLeftIndex],
				normals[topLeftIndex],
				normals[topRightIndex],
			};
			std::vector<Point3f> triangle2 =
			{
				vertexPositions[bottomLeftIndex],
				vertexPositions[topRightIndex],
				vertexPositions[bottomRightIndex],
			};
			std::vector<Normal3f> normal2 =
			{
				normals[bottomLeftIndex],
				normals[topRightIndex],
				normals[bottomRightIndex],
			};

			int index = offset(i, j, 0);
			if (!voxels[index]) {
				// Allocate new voxel and store primitive in it
				voxels[index] = voxelArena.Alloc<Voxel>();
				*voxels[index] = Voxel(
					heightfield->ObjectToWorld, 
					std::move(triangle1),
					std::move(normal1));
				(*voxels[index]).AddTriangle(std::move(triangle2), std::move(normal2));
			}
			else {
				// Add primitive to already-allocated voxel
				voxels[index]->AddTriangle(std::move(triangle1), std::move(normal1));
				voxels[index]->AddTriangle(std::move(triangle2), std::move(normal2));
			}
		}
	}
}

Heightfield::GridAccel::~GridAccel() {
	for (int i = 0; i < nVoxels[0] * nVoxels[1] * nVoxels[2]; ++i)
		if (voxels[i]) voxels[i]->~Voxel();
	FreeAligned(voxels);
}

bool Heightfield::GridAccel::Intersect(
	const Ray &ray, 
	SurfaceInteraction *isect) const {
	float rayT = 0.0f;
	Ray curRay = ray;
	
	Float hitt0, hitt1;
	if (!bounds.IntersectP(curRay, &hitt0, &hitt1)) {
		return false;
	}
	rayT = hitt0;
	Point3f gridIntersect = curRay(rayT);

	Float NextCrossingT[3], DeltaT[3];
	int Step[3], Out[3], Pos[3];
	for (int axis = 0; axis < 3; ++axis) {
		if (curRay.d[axis] == -0.f) curRay.d[axis] = 0.f;
		// Compute current voxel for axis
		Pos[axis] = posToVoxel(gridIntersect, axis);
		if (curRay.d[axis] >= 0) {
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT +
				(voxelToPos(Pos[axis] + 1, axis) - gridIntersect[axis]) / curRay.d[axis];
			DeltaT[axis] = width[axis] / curRay.d[axis];
			Step[axis] = 1;
			Out[axis] = nVoxels[axis];
		}
		else {
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT +
				(voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / curRay.d[axis];
			DeltaT[axis] = -width[axis] / curRay.d[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}

	bool hitSomething = false;
	float minT = MaxFloat;
	for (;;) {
		// Check for intersection in current voxel and advance to next
		Voxel *voxel = voxels[offset(Pos[0], Pos[1], Pos[2])];
		// PBRT_GRID_RAY_TRAVERSED_VOXEL(Pos, voxel ? voxel->size() : 0);
		if (voxel != NULL) {
			hitSomething |= voxel->Intersect(curRay, &minT, isect);
		}

		// Advance to next voxel

		// Find _stepAxis_ for stepping to next voxel
		int bits = ((NextCrossingT[0] < NextCrossingT[1]) << 2) +
			((NextCrossingT[0] < NextCrossingT[2]) << 1) +
			((NextCrossingT[1] < NextCrossingT[2]));
		const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
		int stepAxis = cmpToAxis[bits];
		if (curRay.tMax < NextCrossingT[stepAxis])
			break;
		Pos[stepAxis] += Step[stepAxis];
		if (Pos[stepAxis] == Out[stepAxis])
			break;
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}
	return hitSomething;
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

	bool hitGrid = grid->Intersect(rayOS, isect);
	if (!hitGrid) return false;
	*tHit = (isect->p.x - ray.o.x) / ray.d.x;
	
	return true;
}

bool Heightfield::IntersectP(const Ray &ray, bool testAlphaTexture) const {
	// TODO: Impl
	return false;
}

Float Heightfield::Area() const { return 1.0f; }

Interaction Heightfield::Sample(const Point2f &u, Float *pdf) const {
    Error("Sample of Heightfield is not implemented");
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
