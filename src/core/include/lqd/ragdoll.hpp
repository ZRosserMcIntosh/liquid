// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// The crash dummy as a ragdoll: every capsule is two point masses (its share
// of the dummy's mass, from its volume at the simulant's density) held at its
// length, the capsules joined where they meet (ball joints at the knees and
// elbows, pinned shoulders and hips), the chest, pelvis, neck and head held as
// one rigid trunk, and knees and elbows kept from folding flat. It is
// integrated by position-based dynamics (Verlet with constraint projection,
// Jakobsen 2001; Mueller et al. 2007) against the static scene with friction.
//
// It stands on a release stand until it is hit (a stunt rig), then falls
// under gravity with the momentum the rounds left in it. Wounds and stains on
// it ride along: each is attached to its capsule by its position along the
// axis and its angle round it, measured from the trunk's facing.
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

class Ragdoll {
public:
    void build(const DummyDef& def, const std::vector<int>& colliders, const std::vector<Collider>& all);
    bool active() const { return !seg_.empty(); }
    bool released() const { return released_; }
    void release(double time);
    double releaseTime() const { return releaseTime_; }
    double mass() const { return mass_; }

    // Advances by dt against the colliders that are not the dummy's own.
    void step(double dt, const std::vector<Collider>& colliders);
    // Moves the dummy's capsule colliders to the ragdoll.
    void pose(std::vector<Collider>& colliders) const;
    // An impulse J (N s) at world point `at` on segment s.
    void impulse(int segment, const Vec3& at, const Vec3& J);
    int segmentOf(int collider) const;

    Attachment attach(int segment, const Vec3& point, const Vec3& normal) const;
    void place(const Attachment& a, Vec3& point, Vec3& normal) const;
    Vec3 velocityOf(const Attachment& a) const;
    bool moving() const;

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
        double rest = 0;
        bool minOnly = false;  // only keeps them apart (fold limits)
    };
    void frame(int s, Vec3& e, Vec3& ref, Vec3& side) const;
    Vec3 facing() const;
    void solveLinks();
    void collide(const std::vector<Collider>& colliders);

    std::vector<Vec3> x_, prev_, vel_;
    std::vector<double> w_;  // inverse masses
    std::vector<Seg> seg_;
    std::vector<Link> links_;
    int chestL_ = -1, chestR_ = -1, pelvis_ = -1, neck_ = -1;
    bool released_ = false;
    double releaseTime_ = -1, mass_ = 0, h_ = 1.0 / 480, still_ = 0;
};

}  // namespace lqd
