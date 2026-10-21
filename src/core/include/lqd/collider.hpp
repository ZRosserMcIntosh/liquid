// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Solid shapes the liquid, the droplets and the glass meet: floors, tables,
// walls, glass panes (with the holes rounds punch through them), the crash
// dummy's segments, vessels and their shards, and the shards of a broken pane.
// Every shape answers one question: signed distance and outward normal at a
// world point.
#pragma once

#include <string>
#include <vector>

#include "lqd/library.hpp"
#include "lqd/math.hpp"
#include "lqd/vessel.hpp"

namespace lqd {

enum class ShapeKind { Plane, Box, Capsule, Shell, Slab };

enum class Role { Floor, Table, Wall, Pane, Dummy, Vessel, Shard, Other };
const char* roleName(Role r);

struct Hole {
    double x = 0, y = 0, radius = 0;  // pane frame, metres from the pane centre
};

// A pane in pieces (from Glass Dynamics): which piece of glass each cell of
// the pane holds, for the glass ply nearest each face, and whether that piece
// is still in the frame. A monolithic pane is open wherever its piece has
// gone; a laminate is held by its interlayer and stays closed but for the
// holes punched in it.
struct PaneMap {
    int nx = 0, ny = 0;
    std::vector<int> front, back;      // piece per cell, rows from the bottom; -1 no glass
    std::vector<unsigned char> held;   // per piece
    bool laminated = false;
    double cell = 0, x0 = 0, y0 = 0;   // cell size; the lower left corner in the pane frame
    // The piece at a point in the pane frame (from its centre) on the given face, or -1.
    int pieceAt(double x, double y, bool frontFace) const;
    bool open(double x, double y) const;  // no glass there on either face
};

struct Collider {
    int id = -1;
    std::string name;
    Role role = Role::Other;
    ShapeKind kind = ShapeKind::Plane;
    // Placement (world). Plane: the surface through pos with normal rot*(0,1,0),
    // solid below. Box: centred at pos, half extents `half` along rot's axes.
    // Capsule: segment a-b (world), radius. Shell: a vessel or shard, posed by its body.
    // Slab: a flat piece of a pane, the polygon `outline` (x, y pairs in its own
    // frame about pos) extruded by +-half.z along its local z; the host moves it.
    Vec3 pos;
    Quat rot;
    Vec3 half{0.5, 0.5, 0.5};
    Vec3 a, b;
    double radius = 0;
    const VesselGeometry* shell = nullptr;
    const ShardPattern* shards = nullptr;
    int shard = -1;
    Vec3 comLocal;  // vessel-frame point the body's position refers to
    int body = -1;  // rigid body that moves it, or -1 for static
    const SurfaceDef* surface = nullptr;
    bool enabled = true;
    bool stains = true;  // keeps the drops that land on it as stains
    std::vector<Hole> holes;  // pane: through-holes along its local z
    std::string glazing;      // pane: Glass Dynamics glazing id
    const PaneMap* map = nullptr;  // pane in pieces
    std::vector<double> outline;   // slab
    std::vector<std::vector<double>> outlineHoles;  // slab: loops cut out of it
    int pane = -1, piece = -1;     // slab: the pane and piece it is
    Vec3 vel;                      // slab: m/s
    double boundRadius = 1e30;
    Vec3 color{0.7, 0.7, 0.7};

    // Signed distance (negative inside) and the outward unit normal.
    double sdf(const Vec3& p, Vec3* normal = nullptr) const;
    // World point to the vessel frame (shells).
    Vec3 toVessel(const Vec3& p) const { return comLocal + rot.unrotate(p - pos); }
    Vec3 fromVessel(const Vec3& q) const { return pos + rot.rotate(q - comLocal); }
    // Pane frame (local x, y in the plane, z the normal), world point.
    Vec3 toLocal(const Vec3& p) const { return rot.unrotate(p - pos); }
    bool inHole(const Vec3& local) const;
};

}  // namespace lqd
