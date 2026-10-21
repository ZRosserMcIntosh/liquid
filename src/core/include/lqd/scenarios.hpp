// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Ready-made scenes, each configurable from a JSON object (the browser lab
// and lqdemo use the same builder):
//   wine_glass   a filled glass knocked off a table onto the floor
//   bottle_drop  an open bottle of wine dropped from a height
//   gunshot      a round through the crash dummy, with an optional glass pane
//                in front of it or behind it, the spatter on the wall and floor
//   drip         passive drops from a row of heights (drip-stain card)
//   angles       single drops thrown at a wall from 10 to 90 degrees (impact-angle card)
//   dam_break    a water column released in a tank (Martin & Moyce 1952)
#pragma once

#include <cstdint>
#include <string>

#include "lqd/world.hpp"

namespace lqd {

struct SceneConfig {
    std::string scenario = "wine_glass";
    int tier = 2;
    std::uint64_t seed = 1;
    std::string liquid;                 // default: the vessel's or the dummy's
    std::string vessel;                 // default per scenario
    std::string floor = "ceramic_tile";
    std::string table = "oak_table";
    std::string wall = "painted_drywall";
    double fill = -1;                   // m^3, default from the vessel
    double knock = 0.9;                 // m/s, wine_glass: speed of the knock at the bowl
    double height = 1.0;                // m, bottle_drop: base height; drip: highest drop
    double tilt = 25;                   // deg, bottle_drop
    // gunshot
    std::string projectile = "pistol_9mm_fmj";
    std::string dummy = "crash_dummy";
    std::string target = "chest";       // chest, abdomen, shoulder, thigh
    double distance = 3.0;              // m, muzzle to the dummy
    double shooterHeight = 1.4;         // m
    double speed = -1;                  // m/s, default from the projectile
    std::string pane = "none";          // none, front, behind
    std::string glazing = "annealed_6";
    double paneWidth = 1.0, paneHeight = 1.2, paneThickness = 0.006;  // m
    bool deferPanes = false;            // the host resolves pane strikes (Glass Dynamics)
};

// Parses a scene config from JSON (unknown keys are errors).
bool parseSceneConfig(const std::string& json, SceneConfig& out, std::string* error);

// Builds the scene into an empty world. Returns false with a message on error.
bool buildScene(World& w, const SceneConfig& c, std::string* error);

// What lies on the path of the shot: the muzzle, the aim point and the target zone's name.
struct ShotLine {
    Vec3 muzzle, aim;
};
ShotLine shotLine(const SceneConfig& c);

}  // namespace lqd
