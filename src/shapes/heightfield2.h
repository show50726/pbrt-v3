
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

#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield2.h*
#include "shape.h"

namespace pbrt {

// Heightfield Declarations
class Heightfield : public Shape {
  public:
    // Heightfield Public Methods
    Heightfield(const Transform *ObjectToWorld, const Transform *WorldToObject,
                bool reverseOrientation,
                int nu,
                int nv, const float *zs);
    ~Heightfield();
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    Float Area() const;
    Interaction Sample(const Point2f &u, Float *pdf) const;

	struct Voxel {
		// Voxel Public Methods
		uint32_t size() const { return triangles.size(); }
		Voxel() { }
		Voxel(const Transform * objToWorld, 
			std::vector<Point3f> op,
			std::vector<Normal3f> nm) : 
			objToWorld (objToWorld) {
			triangles.push_back(op);
			normals.push_back(nm);
		}
		void AddTriangle(std::vector<Point3f> prim, std::vector<Normal3f> nm) {
			triangles.push_back(prim);
			normals.push_back(nm);
		}
		bool Intersect(const Ray &ray, Float* tHit, SurfaceInteraction *isect);
		bool IntersectP(const Ray &ray);
		bool Intersect(const Ray &ray, Float* tHit, std::vector<Point3f> triangle, std::vector<Normal3f> normal, SurfaceInteraction *isect);
		bool IntersectP(const Ray &ray, std::vector<Point3f> triangle);
	private:
		std::vector<std::vector<Point3f>> triangles;
		std::vector<std::vector<Normal3f>> normals;
		const Transform * objToWorld;
	};

	// GridAccel Declarations
	class GridAccel {
	public:
		// GridAccel Public Methods
		GridAccel(Heightfield* hf);
		~GridAccel();
		bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
		bool IntersectP(const Ray &ray) const;

	private:
		// GridAccel Private Methods
		int posToVoxel(const Point3f &P, int axis) const {
			int v = static_cast<int>((P[axis] - bounds.pMin[axis]) *
				invWidth[axis]);
			return Clamp(v, 0, nVoxels[axis] - 1);
		}
		float voxelToPos(int p, int axis) const {
			return bounds.pMin[axis] + p * width[axis];
		}
		inline int offset(int x, int y, int z) const {
			return z * nVoxels[0] * nVoxels[1] + y * nVoxels[0] + x;
		}
		inline int vertexOffset(int x, int y) const {
			return x + y * heightfield->nx;
		}

		int nVoxels[3];
		Heightfield* heightfield;
		Bounds3f bounds;
		Vector3f width, invWidth;
		Voxel **voxels;
		Point3f *vertexPositions;
		Normal3f *normals;
		MemoryArena voxelArena;
	};


  private:
    // Heightfield Private Data
    float *z;
    int nx, ny;
	std::unique_ptr<GridAccel> grid;
};

std::shared_ptr<Heightfield> CreateHeightfield2(const Transform *o2w,
                                                  const Transform *w2o,
                                    bool reverseOrientation,
                                    const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_HEIGHTFIELD2_H
