// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// A glass vessel as a revolved wall: the mid-surface profile of a VesselDef
// turned about +y. From it come the exact signed distance to the glass (for
// the fluid and for contacts), the cavity and its fill level, the mass and
// inertia (Pappus' theorems on the wall section), a render mesh with real
// thickness, and the shards when it breaks: Voronoi cells of the developed
// wall, crowded round the impact.
#pragma once

#include <cstdint>
#include <vector>

#include "lqd/library.hpp"
#include "lqd/math.hpp"
#include "lqd/rng.hpp"

namespace lqd {

struct WallPoint {
    double dist = 0;  // distance in the meridian plane to the mid-surface
    double s = 0;     // arc length from the rim along the profile
    double theta = 0; // azimuth
    double t = 0;     // wall thickness there
    double side = 1;  // +1 on the outer side (away from the cavity), -1 inside
    Vec3 normal;      // unit, from the mid-surface towards the query point
};

// Contact points on the glass surface, in the vessel frame.
struct WallSample {
    Vec3 p;
    double s = 0, theta = 0;
    double radius = 0.01;  // Hertz-equivalent radius of the surface there
    double reach = 0.005;  // bending reach of the wall, 0.76 sqrt(R t)
    double t = 0.001;
    double shellRadius = 0.05;  // the wall's larger principal radius (shell stiffness)
    Vec3 wallNormal{1, 0, 0};   // unit normal of the wall there (vessel frame)
    bool edge = false;     // a free edge (rim, foot rim, lip)
};

// A triangle mesh in the vessel frame with the wall cell each triangle belongs to.
struct VesselMesh {
    std::vector<Vec3> pos, nrm;
    std::vector<std::uint32_t> tri;
};

class VesselGeometry {
public:
    VesselGeometry() = default;
    explicit VesselGeometry(const VesselDef& def, int columns = 64);

    const VesselDef& def() const { return *def_; }
    double length() const { return s_.back(); }
    double height() const { return height_; }
    double maxRadius() const { return maxR_; }
    double rimHeight() const { return rimY_; }

    // Nearest point of the wall's mid-surface to p (vessel frame).
    WallPoint nearest(const Vec3& p) const;
    // Signed distance to the glass (negative inside the glass) and its normal.
    double sdf(const Vec3& p, Vec3* normal = nullptr, WallPoint* wp = nullptr) const;

    // The cavity in the meridian plane: inside test with a clearance from the wall.
    bool inCavity(double r, double y, double margin) const;
    // Cavity volume below height y (m^3), and the height holding volume v.
    double cavityVolume(double y) const;
    double cavityCapacity() const { return cavityVolume(rimY_); }
    double fillHeight(double volume) const;

    // Glass mass, centre of mass (on the axis) and inertia about it (vessel axes).
    double mass() const { return mass_; }
    Vec3 centreOfMass() const { return {0, comY_, 0}; }
    Mat3 inertia() const { return inertia_; }

    // Profile resampling used by the mesh and the shards: rows along s,
    // columns round the axis.
    int rows() const { return static_cast<int>(rowS_.size()); }
    int columns() const { return columns_; }
    double rowS(int i) const { return rowS_[i]; }
    // Mid-surface (r, y), unit outward normal (nr, ny) and thickness at arc length s.
    void at(double s, double& r, double& y, double& nr, double& ny, double& t) const;

    // Wall cells (row, column) -> render mesh of the cells for which keep(cell) is
    // true, with crack faces where a kept cell meets a dropped one.
    template <class Keep>
    VesselMesh mesh(Keep&& keep) const;
    VesselMesh mesh() const {
        return mesh([](int) { return true; });
    }
    int cellOf(int row, int col) const { return row * columns_ + ((col % columns_) + columns_) % columns_; }
    int cells() const { return (rows() - 1) * columns_; }
    // Centre of a cell in (s, theta) and its wall volume.
    void cellCentre(int cell, double& s, double& theta) const;
    double cellVolume(int cell) const;
    Vec3 cellPoint(int cell) const;

    const std::vector<WallSample>& samples() const { return samples_; }

private:
    void resample();
    void massProperties();
    void buildSamples();
    void cavityPolygon();
    VesselMesh buildMesh(const std::vector<char>& kept) const;

    const VesselDef* def_ = nullptr;
    std::vector<ProfilePoint> p_;
    std::vector<double> s_;
    std::vector<double> rowS_;
    std::vector<double> rowR_, rowY_, rowNr_, rowNy_, rowT_;
    std::vector<double> polyR_, polyY_;  // cavity polygon (closed along the axis)
    std::vector<WallSample> samples_;
    int columns_ = 64;
    double height_ = 0, maxR_ = 0, rimY_ = 0;
    double mass_ = 0, comY_ = 0;
    Mat3 inertia_;
};

template <class Keep>
VesselMesh VesselGeometry::mesh(Keep&& keep) const {
    std::vector<char> kept(cells());
    for (int c = 0; c < cells(); ++c) kept[c] = keep(c) ? 1 : 0;
    return buildMesh(kept);
}

// A broken vessel: every wall cell belongs to one shard.
struct ShardPattern {
    std::vector<int> cellShard;  // per wall cell
    int count = 0;
    // Seeds in the developed wall (s, theta) for nearest-seed lookups.
    std::vector<double> seedS, seedTheta;
    int shardAt(const VesselGeometry& g, double s, double theta) const;
};

// Voronoi pattern: `count` seeds, most of them crowded round the impact point
// (s0, theta0) within `spread` metres of developed wall, the rest anywhere.
ShardPattern makeShards(const VesselGeometry& g, double s0, double theta0, int count, double spread, Rng& rng);

}  // namespace lqd
