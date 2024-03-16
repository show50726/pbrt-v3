
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
#include "sampling.h"
#include "shapes/heightfield2.h"
#include "stats.h"
#include "error.h"

namespace pbrt {

Heightfield::Heightfield(const Transform *ObjectToWorld,
                         const Transform *WorldToObject,
                         bool reverseOrientation, int nu, int nv,
                         const float *zs)
    : Shape(ObjectToWorld, WorldToObject, reverseOrientation), nx(nu), ny(nv) {
    z = new float[nx * ny];
    memcpy(z, zs, nx * ny * sizeof(float));
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
    Error("Not implemented");
    return false;
}

bool Heightfield::IntersectP(const Ray &ray, bool testAlphaTexture) const {
    Error("Not implemented");
    return false;
}

Float Heightfield::Area() const { return 1.0f; }

Interaction Heightfield::Sample(const Point2f &u, Float *pdf) const {
    Error("Not implemented");
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
