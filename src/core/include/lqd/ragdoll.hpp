// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// The crash dummy as a ragdoll: every capsule is two point masses (its share
// of the dummy's mass, from its volume at the simulant's density) held at its
// length, the capsules joined where they meet, and every joint kept within a
// human range of motion: the spine and neck flex a little, hips and shoulders
// swing within elliptical cones (flexion, extension, abduction, adduction),
// knees and elbows are hinges that bend one way and do not hyperextend, and
// the ankles rock within a range. The limbs collide with each other and with
// the trunk. It is integrated by position-based dynamics (Verlet with
// constraint projection, Jakobsen 2001; Mueller et al. 2007) against the static
// scene with friction. The ranges are the usual adult ones (AAOS, Norkin &
// White), a little narrowed for a dummy.
//
// It stands on a release stand until it is hit (a stunt rig), then falls
// under gravity with the momentum the rounds left in it. Wounds and stains on
// it ride along: each is attached to its capsule by its position along the
// axis and its angle round it, measured from the trunk's facing.
//
// Heavy blunt things (impactors: a hammer, a wrecking ball, a car) push it
// through inelastic contacts that conserve momentum between the impactor and
// the point mass they strike. It keeps an energy book: what it was given,
// what it absorbed (in the blow, and afterwards in landing, sliding and its
// joints), and what still moves it.
#pragma once

#include <string>
#include <vector>

#include "lqd/collider.hpp"
#include "lqd/library.hpp"
#include "lqd/math.hpp"

namespace lqd {

struct Attachment {
    int segment = -1;
    double t = 0;       // along the capsule axis, 0 at a, 1 at b
    double angle = 0;   // round the axis from the trunk's facing
    double radial = 0;  // distance from the axis, m
    double normalAngle = 0, normalAxial = 0;  // the surface normal, the same way
};

// An impactor as the ragdoll sees it: its collider, mass and velocity, which
// the contacts change.
struct Pusher {
    int collider = -1;
    double mass = 0;
    Vec3 vel;
    bool touched = false;  // it touched the dummy this step
};

// What one segment absorbed in contacts during the last step, and where most of it went in.
struct SegmentHit {
    double energy = 0;      // J
    double best = 0;        // J, the contact that took the most
    Vec3 point, normal;     // that contact on the dummy's surface, the normal out of what it struck
    double speed = 0;       // m/s, the fastest closing speed
    int collider = -1;
};

// The dummy's energy since it left its stand: start + delivered = mechanical + absorbed.
struct EnergyBook {
    double start = 0;       // J, kinetic and potential (above y = 0) when it was released
    double fromHits = 0;    // J, kinetic energy the rounds' pushes gave it
    double fromImpactors = 0;  // J, work the impactors did on it
    double inBlows = 0;     // J, absorbed while an impactor was in contact
    double afterwards = 0;  // J, absorbed otherwise: landing, sliding, its joints
    double peakKinetic = 0; // J
    double delivered() const { return fromHits + fromImpactors; }
    double absorbed() const { return inBlows + afterwards; }
};

class Ragdoll {
public:
    void build(const DummyDef& def, const std::vector<int>& colliders, const std::vector<Collider>& all);
    bool active() const { return !seg_.empty(); }
    bool released() const { return released_; }
    void release(double time);
    double releaseTime() const { return releaseTime_; }
    double mass() const { return mass_; }

    // Advances by dt against the colliders that are not the dummy's own, and
    // the impactors in `pushers` (their velocities change by what they hit).
    void step(double dt, const std::vector<Collider>& colliders, std::vector<Pusher>* pushers = nullptr);
    // Moves the dummy's capsule colliders to the ragdoll.
    void pose(std::vector<Collider>& colliders) const;
    // An impulse J (N s) at world point `at` on segment s.
    void impulse(int segment, const Vec3& at, const Vec3& J);
    int segmentOf(int collider) const;

    Attachment attach(int segment, const Vec3& point, const Vec3& normal) const;
    void place(const Attachment& a, Vec3& point, Vec3& normal) const;
    Vec3 velocityOf(const Attachment& a) const;
    bool moving() const;

    // Every point moving at v (a stunt fall's run-up); and the whole body
    // turned by q about `pivot` (before it is released).
    void setVelocity(const Vec3& v);
    void turn(const Quat& q, const Vec3& pivot);
    // Its parts never move faster than this (default 15 m/s; see step()).
    void setSpeedCap(double v) { speedCap_ = v; }
    double kinetic() const;
    double potential() const;
    const EnergyBook& energy() const { return book_; }
    const std::vector<SegmentHit>& hits() const { return hits_; }
    double segmentMass(int s) const;
    const std::string& segmentName(int s) const { return seg_[s].name; }

    // Per segment: end points and radius (for drawing and export).
    int segments() const { return static_cast<int>(seg_.size()); }
    Vec3 a(int s) const { return x_[seg_[s].i]; }
    Vec3 b(int s) const { return x_[seg_[s].j]; }
    double radius(int s) const { return seg_[s].r; }

private:
    struct Seg {
        int i = 0, j = 0, collider = -1;
        double r = 0;
        std::string name;
    };
    struct Link {
        int i = 0, j = 0;
        double lo = 0, hi = 0;  // kept within [lo, hi] apart (lo == hi: rigid)
        double stiffness = 1;   // of the way back into range each pass (below 1: a soft limit)
    };
    // A joint's range: the child's direction (near -> far point) in its
    // parent's frame. Cone: flexion phi (forward +) and abduction psi
    // (outward +) about the rest direction -up. Hinge: bend about the
    // parent's axis towards bendDir (forward + or backward -), out of plane c.
    struct Limit {
        bool hinge = false;
        int near = 0, far = 0;      // the child's points
        int parentA = 0, parentB = 0;  // hinge: the parent segment's points (A to B ends at the joint)
        bool pelvisFrame = false;   // cone: the pelvis's frame, else the thorax's
        double outward = 1;         // +1 right, -1 left
        double bendSign = 1;        // hinge: +1 bends forward (elbow), -1 backward (knee)
        double lo1 = 0, hi1 = 0, lo2 = 0, hi2 = 0;
    };
    struct Pair {
        int a = 0, b = 0;           // segments that may collide
        double gap = 0;             // how close they may come (their radii, or less if closer at rest)
    };
    void thorax(Vec3& up, Vec3& lat, Vec3& fwd) const;
    void pelvisFrame(Vec3& up, Vec3& lat, Vec3& fwd) const;
    void solveLimits();
    void selfCollide();
    std::vector<Limit> limits_;
    std::vector<Pair> pairs_;
    void frame(int s, Vec3& e, Vec3& ref, Vec3& side) const;
    Vec3 facing() const;
    void solveLinks();
    void collide(const std::vector<Collider>& colliders, std::vector<Pusher>* pushers, double& work, bool exchange);

    std::vector<Vec3> x_, prev_, vel_;
    std::vector<double> w_;  // inverse masses
    std::vector<Seg> seg_;
    std::vector<Link> links_;
    int chestL_ = -1, chestR_ = -1, pelvis_ = -1, neck_ = -1;
    std::vector<SegmentHit> hits_;
    EnergyBook book_;
    bool released_ = false;
    static constexpr double kStep = 1.0 / 960;  // s
    double releaseTime_ = -1, mass_ = 0, h_ = kStep, still_ = 0, speedCap_ = 15, acc_ = 0;
};

}  // namespace lqd
