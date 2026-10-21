// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// The bulk liquid: smoothed-particle hydrodynamics kept incompressible by
// position-based density constraints (Macklin & Mueller 2013), with
//   - walls by the analytic kernel volume Phi(d) (see kernel.hpp), so floors,
//     tables, glass walls and shards need no wall particles;
//   - the liquid's real viscosity, shear-thinning where the liquid is
//     (Carreau-Yasuda on the local strain rate), as the SPH Laplacian of
//     Morris et al. (1997);
//   - surface tension as the continuum surface force sigma kappa n
//     (Brackbill et al. 1992; Morris 2000) with the liquid's real sigma;
//   - two-way coupling: every push the liquid takes from a moving wall goes
//     back to that wall's body as an equal and opposite impulse (for bodies
//     at least ~200 particles heavy; lighter ones, below the liquid's
//     resolution, blend to one-way: they move the liquid, it barely moves them).
// Particles that end up alone in the air leave the bulk as droplets (spray.hpp).
#pragma once

#include <vector>

#include "lqd/collider.hpp"
#include "lqd/grid.hpp"
#include "lqd/kernel.hpp"
#include "lqd/library.hpp"
#include "lqd/rigid.hpp"

namespace lqd {

struct FluidSettings {
    double spacing = 0.005;      // m, particle spacing at rest
    int iterations = 4;          // density-constraint iterations per step
    double xsph = 0.01;          // numerical velocity smoothing (dimensionless)
    double surfaceTension = 1;   // multiplier on the liquid's sigma (1: the real value)
    double wallFriction = 0;     // tangential damping on wall contact per step (0: free slip)
    double pushSpeed = 5;        // m/s, the most one constraint iteration may move a particle per second
    double sprayDensity = 0.1;   // a particle whose density is within this fraction of rest of its own
                                 // kernel's (nearly alone) leaves the bulk as a droplet
};

struct WallTouch {
    int collider = -1;
    double d = 0;
    Vec3 n;
};

class Fluid {
public:
    void configure(const LiquidDef& liquid, const FluidSettings& s);
    const LiquidDef* liquid() const { return liquid_; }
    const FluidSettings& settings() const { return set_; }

    int add(const Vec3& x, const Vec3& v);
    void remove(const std::vector<int>& sortedIndices);
    std::size_t size() const { return x.size(); }

    // One step of dt. Moving walls get their reaction impulses in
    // bodies[...].fluidImpulse/fluidAngular. Indices of particles that should
    // become droplets are returned in `spray` (ascending).
    void step(double dt, const std::vector<Collider>& colliders, std::vector<RigidBody>& bodies, std::vector<int>* spray);

    // The largest stable step for the current state.
    double maxStep() const;
    double maxSpeed() const;
    double particleMass() const { return mass_; }
    double restDensity() const { return rho0_; }
    double support() const { return H_; }
    double volume() const { return mass_ * static_cast<double>(x.size()) / rho0_; }
    // Mean and largest relative compression in the last step.
    double meanCompression() const { return meanErr_; }
    double maxCompression() const { return maxErr_; }
    double kineticEnergy() const;
    // How much of the liquid's push a body takes back (0 one-way .. 1 two-way).
    double coupling(const RigidBody& b) const;
    // Liquid density at a point (fluid particles only), from the last step's grid.
    double densityAt(const Vec3& p) const;
    const NeighbourGrid& grid() const { return grid_; }
    const Kernel& kernel() const { return kernel_; }

    std::vector<Vec3> x, v;
    std::vector<double> rho;

private:
    void findWalls(const std::vector<Collider>& colliders);

    const LiquidDef* liquid_ = nullptr;
    FluidSettings set_;
    Kernel kernel_;
    double H_ = 0.01, mass_ = 0, rho0_ = 1000;
    std::vector<Vec3> xp_, dx_, normal_;
    std::vector<double> lambda_, mu_, fluidRho_;
    std::vector<WallTouch> touch_;  // kMaxTouch per particle
    std::vector<int> nTouch_;
    std::vector<int> off_, nb_;
    NeighbourGrid grid_;
    double meanErr_ = 0, maxErr_ = 0;
    double selfDensity_ = 0.3;
    double wallShift_ = 0;
};

}  // namespace lqd
