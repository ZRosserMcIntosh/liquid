// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Droplets in flight and the stains they leave.
//
// A droplet is a sphere of the liquid flying through still air: gravity and
// the sphere drag curve (Schiller & Naumann 1933 below Re 1000, Newton's 0.44
// above), integrated semi-implicitly so the smallest mist is stable at any
// step. Above an air Weber number of 12 a drop breaks up after the
// Pilch & Erdman (1987) breakup time into drops at a Weber number of 6.
//
// Where a drop meets a surface it leaves a stain:
//   - width from the maximum spreading diameter of Laan et al. (2014), fitted
//     on blood and water: D/d = Re^(1/5) P^(1/2) / (1.24 + P^(1/2)),
//     P = We Re^(-2/5), on the normal component of the impact velocity;
//   - length from the impact angle, L = W / sin(alpha), the relation bloodstain
//     pattern analysts invert (alpha = asin(W/L));
//   - the tail pointing the way the drop travelled;
//   - splashing above Mundo et al.'s (1995) threshold K = We^(1/2) Re^(1/4)
//     > 57.7 (lowered on rough surfaces), with spines and satellite droplets.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "lqd/collider.hpp"
#include "lqd/library.hpp"
#include "lqd/math.hpp"
#include "lqd/rng.hpp"

namespace lqd {

class Fluid;

inline constexpr double kAirDensity = 1.204;      // kg/m^3 at 20 C, 101.325 kPa (published, ISO 2533)
inline constexpr double kAirViscosity = 1.813e-5; // Pa s at 20 C (published)

struct Droplet {
    Vec3 x, v;
    double d = 1e-3;  // diameter, m
    double t = 0;     // its own clock (it may have been born inside a step)
    double breakup = 0;
    int source = -1;
    bool fromBulk = false;  // a particle of the bulk liquid gone airborne
    bool satellite = false; // thrown from another drop's splash
};

struct Stain {
    Vec3 pos, normal, dir;  // centre, surface normal, unit tail direction in the surface
    double width = 0, length = 0;
    double alpha = 0;       // true impact angle to the surface, rad
    double diameter = 0;    // of the drop
    double speed = 0;       // of the drop at impact
    double volume = 0;
    double time = 0;
    int collider = -1;
    int source = -1;
    int spines = 0;
    bool satellite = false, splashed = false;
};

struct Source {
    std::string kind;  // "back spatter", "forward spatter", "drip", "splash", "pour"
    Vec3 origin;
    double time = 0;
    double volume = 0;
    int drops = 0;
};

// Sphere drag coefficient at Reynolds number Re.
double dragCoefficient(double Re);
// Terminal speed of a drop of diameter d (m) and density rho (kg/m^3) in air.
double terminalSpeed(double d, double rho);

struct StainShape {
    double width = 0, length = 0, alpha = 0;
    double K = 0, Kc = 0;
    bool splash = false;
    int spines = 0;
    double satelliteFraction = 0;
    int satellites = 0;
};
// The stain a drop of diameter d leaves arriving with normal speed vn and
// tangential speed vt on the surface.
StainShape stainShape(const LiquidDef& liquid, double d, double vn, double vt, const SurfaceDef* surface);

class Spray {
public:
    void configure(const LiquidDef& liquid) { liquid_ = &liquid; }
    const LiquidDef* liquid() const { return liquid_; }

    int addSource(const std::string& kind, const Vec3& origin, double time);
    void emit(const Droplet& d);

    // Advances every drop to time tEnd. Drops that land on liquid or that
    // came from the bulk and land anywhere are handed back in `rejoin`.
    void step(double tEnd, const std::vector<Collider>& colliders, const Fluid* fluid, std::vector<Droplet>* rejoin,
              Rng& rng);

    std::vector<Droplet> drops;
    std::vector<Stain> stains;
    std::vector<Source> sources;
    // Volume bookkeeping (m^3).
    double stainedVolume = 0, absorbedVolume = 0, lostVolume = 0, bulkReturned = 0;
    // Drops that landed in the bulk liquid, waiting to make up a particle.
    double toBulk = 0;
    Vec3 toBulkAt, toBulkV;
    // A surface holds a thin film where drops land; what lands beyond it pools
    // into bulk liquid (a puddle on a floor, a run on a wall). The world sets
    // the particle volume; new particles to make are listed in `pooled`.
    double particleVolume = 0;
    double filmThickness = 3e-4;  // m (estimated: a wet film before it flows)
    std::vector<std::pair<Vec3, Vec3>> pooled;  // position, surface normal
    int breakups = 0, splashes = 0;
    int maxDrops = 20000;
    double airVolume() const;

private:
    void land(const Droplet& d, const Collider& c, const Vec3& at, const Vec3& n, double t, std::vector<Droplet>& born,
              Rng& rng);
    void pool(const Collider& c, const Vec3& at, const Vec3& n, double volume);
    const LiquidDef* liquid_ = nullptr;
    struct Cell {
        double volume = 0;
    };
    std::unordered_map<std::uint64_t, Cell> cells_;
};

}  // namespace lqd
