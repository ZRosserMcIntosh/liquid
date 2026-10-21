// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Getting a simulation into Unreal, Unity, Blender or Houdini.
//
//   - Binary glTF 2.0 (.glb): the liquid's surface as a flipbook (one mesh per
//     recorded frame, each shown for its frame by a STEP scale animation), the
//     droplets with it, every glass body and shard as a mesh with its recorded
//     translation and rotation, the pieces of a broken pane the same way with
//     the stains on them as their children, the stains as flat decal meshes that
//     appear when they land, and optionally the static scene. glTF is metres, right-handed,
//     +y up, the core's own frame; the engines' glTF importers convert it.
//   - PLY point caches (one per frame): bulk particles and droplets with
//     velocity and radius, for Houdini, Blender or a Niagara/VFX Graph importer.
//   - Stains as JSON: centre, normal, tail direction, width, length, impact
//     angle, time and surface, for placing decals procedurally.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lqd/mesher.hpp"
#include "lqd/world.hpp"

namespace lqd {

struct RecordedFrame {
    double time = 0;
    std::vector<float> particles;  // xyz
    std::vector<float> drops;      // xyz + diameter
    std::vector<float> poses;      // per body: x y z, qx qy qz qw
    std::vector<float> dummy;      // per ragdoll capsule: a xyz, b xyz
    std::vector<std::vector<float>> pieces;  // per broken pane, per piece: x y z, qx qy qz qw, shown
    std::size_t stains = 0;
};

class Recorder {
public:
    explicit Recorder(double fps = 60, std::size_t maxFrames = 600) : fps_(fps), max_(maxFrames) {}
    // Records a frame when the world's clock has reached the next frame time.
    bool capture(const World& w);
    void captureNow(const World& w);
    double fps() const { return fps_; }
    std::vector<RecordedFrame> frames;
    double spacing = 0.005;

private:
    double fps_;
    std::size_t max_;
    double next_ = -1;
};

struct ExportOptions {
    bool scene = true;      // floor, table, walls, pane and dummy as simple meshes
    bool drops = true;
    bool stains = true;
    int stride = 1;         // every n-th recorded frame
};

std::vector<std::uint8_t> exportGlb(const World& w, const Recorder& r, const ExportOptions& o = {});
std::vector<std::uint8_t> exportPly(const RecordedFrame& f, double spacing);
std::string exportObj(const SurfaceMesh& m);
std::string exportStainsJson(const World& w);

// A body's mesh in its own frame (centre of mass at the origin).
VesselMesh bodyMesh(const World& w, int body);

}  // namespace lqd
