// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// The liquid's free surface as a triangle mesh, for drawing and for export:
// a smooth density field splatted from the particles onto a sparse grid, and
// its iso-surface by marching tetrahedra (six per cube, no lookup tables), with
// shared vertices and normals from the field gradient.
#pragma once

#include <cstdint>
#include <vector>

#include "lqd/math.hpp"

namespace lqd {

struct SurfaceMesh {
    std::vector<float> pos, nrm;      // xyz per vertex
    std::vector<std::uint32_t> tri;
    std::size_t vertices() const { return pos.size() / 3; }
    std::size_t triangles() const { return tri.size() / 3; }
};

struct MesherSettings {
    double cell = 0;     // grid cell, m (default: 0.6 x spacing)
    double radius = 0;   // splat radius, m (default: 1.8 x spacing)
    double iso = 0.45;   // fraction of the field of a resting lattice
};

SurfaceMesh meshLiquid(const std::vector<Vec3>& particles, double spacing, const MesherSettings& s = {});

}  // namespace lqd
