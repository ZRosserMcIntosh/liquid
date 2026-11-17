// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Rigid glass bodies (a whole vessel, or one shard of it) and the impact test
// that decides whether a strike breaks the glass.
//
// Contacts are sequential impulses on points of the glass surface against the
// static scene, with speculative contacts so that a fast body does not pass
// through a thin floor, Coulomb friction, and restitution from the surface.
//
// The impact test follows the strike the way a contact mechanic would:
// Hertz's peak force for the effective mass at the struck point on the local
// radius of the glass (a rim or edge is half its thickness), against a hard
// floor, or a compliant layer (carpet pile) taken as a Winkler foundation over
// a hard substrate; the bending
// stress that force raises in the wall (concentrated load on a plate, with the
// logarithmic peak cut off at the shell's bending reach 0.76 sqrt(R t)); and a
// two-parameter Weibull probability of failure at that stress.
#pragma once

#include <vector>

#include "lqd/collider.hpp"
#include "lqd/library.hpp"
#include "lqd/math.hpp"
#include "lqd/vessel.hpp"

namespace lqd {

struct RigidBody {
    int id = -1;
    std::string name;
    double mass = 0.1;
    Mat3 inertia = Mat3::diag(1e-4, 1e-4, 1e-4);  // body axes, about the centre of mass
    Vec3 x, v, w;                                  // centre of mass, velocity, angular velocity (world)
    Quat q;
    const VesselGeometry* geom = nullptr;
    int shard = -1;
    Vec3 comLocal;  // centre of mass in the vessel frame
    std::vector<WallSample> samples;
    int collider = -1;   // its shape in the world's collider list
    bool breakable = false, broken = false;
    bool sleeping = false;
    double still = 0;    // seconds it has been nearly still
    double restitution = 1;
    double born = 0, brokenAt = -1;  // when it appeared, and broke (-1: never)
    Vec3 fluidImpulse, fluidAngular;  // accumulated from the liquid since the last step
    double fluidMass = 0;             // liquid in contact with it in the last step

    Vec3 fromVessel(const Vec3& p) const { return x + q.rotate(p - comLocal); }
    Vec3 velocityAt(const Vec3& p) const { return v + cross(w, p - x); }
    Mat3 invInertiaWorld() const {
        const Mat3 R = Mat3::fromQuat(q);
        return R * inertia.inverse() * R.transposed();
    }
    // Effective mass for an impulse along n at world point p.
    double effectiveMass(const Vec3& p, const Vec3& n) const {
        const Vec3 rn = cross(p - x, n);
        const double k = 1.0 / mass + dot(rn, invInertiaWorld() * rn);
        return 1.0 / k;
    }
    void applyImpulse(const Vec3& J, const Vec3& at) {
        v += J / mass;
        w += invInertiaWorld() * cross(at - x, J);
        sleeping = false;
        still = 0;
    }
    double kineticEnergy() const { return 0.5 * mass * length2(v) + 0.5 * dot(w, Mat3::fromQuat(q) * (inertia * q.unrotate(w))); }
};

// One surface point arriving on a collider fast enough to count as a strike.
struct Strike {
    int body = -1, collider = -1, sample = -1;
    Vec3 point, normal;
    double speed = 0;  // normal approach speed, m/s
};

struct ImpactAssessment {
    double speed = 0, effectiveMass = 0;
    double force = 0, contactRadius = 0, indentation = 0;  // N, m, m
    double stress = 0;                                      // Pa, peak bending stress in the wall
    double probability = 0;                                 // of the wall breaking
    bool edge = false;
};

// `alongNormal` is |cos| of the angle between the strike direction and the
// wall's normal (1: struck square on the face, 0: struck end-on along the wall).
ImpactAssessment assessImpact(const VesselDef& v, const WallSample& s, const SurfaceDef& surface, double meff,
                              double speed, double alongNormal = 1);

// Advances one body by dt against the static colliders (those with body < 0).
// Strikes faster than `strikeSpeed` are appended to `strikes`.
void stepBody(RigidBody& b, const std::vector<Collider>& colliders, double dt, double strikeSpeed,
              std::vector<Strike>* strikes);

}  // namespace lqd
