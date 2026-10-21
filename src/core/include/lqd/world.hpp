// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// A scene: static surfaces, glass bodies, the bulk liquid, droplets and
// stains, rounds in flight and the wounds they leave. advance() runs whole
// internal steps whose size depends only on the state, so the result does not
// depend on how the host chops its frames.
#pragma once

#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lqd/collider.hpp"
#include "lqd/fluid.hpp"
#include "lqd/library.hpp"
#include "lqd/rigid.hpp"
#include "lqd/rng.hpp"
#include "lqd/spray.hpp"
#include "lqd/vessel.hpp"

namespace lqd {

// Resolution per device tier. Lower tiers lower resolution, not physics.
struct Budget {
    const char* name = "mid";
    double spacing = 0.005;   // bulk liquid particle spacing, m
    int iterations = 4;
    int maxDrops = 12000;
    int columns = 64;         // vessel mesh columns
};
Budget budgetForTier(int tier);  // 0 potato .. 4 ultra

struct Event {
    double time = 0;
    std::string text;
};

// A round reaching a glass pane. The host may resolve it with Glass
// Dynamics (see docs/integration.md); otherwise the built-in estimate does.
struct PaneStrike {
    int projectile = -1, collider = -1;
    Vec3 point;               // world
    double x = 0, y = 0;      // pane frame from its lower left corner, m (Glass Dynamics' frame)
    double speed = 0;         // m/s
    double yaw = 0, pitch = 0;  // rad, of the path from the pane normal
    bool fromBack = false;    // arriving on the pane's -z face
    double time = 0;
    std::string glazing;
};
struct PaneResult {
    bool penetrated = true;
    double speedOut = 0;
    double holeRadius = 0.005;  // m
    bool shattered = false;     // the pane is gone (tempered glass crazes and falls)
};

struct Projectile {
    int id = -1;
    const ProjectileDef* def = nullptr;
    Vec3 pos, dir;
    double speed = 0, time = 0;
    enum class State { Waiting, Flying, AtPane, Stopped } state = State::Waiting;
    double launch = 0;          // when it leaves the muzzle
    bool inside = false;        // crossing the dummy
    double exitDistance = 0;    // remaining chord through the dummy, m
    double a = 0, b = 0;        // Poncelet resistance F = a + b v^2 there (N, kg/m)
    double entrySpeed = 0;
    std::vector<Vec3> path;     // corners: muzzle, entries, exits, panes, end
    std::string fate;
};

struct Wound {
    Vec3 pos, normal;
    bool exit = false;
    double time = 0, bled = 0, owed = 0;
    int source = -1;
};

struct VolumeBook {
    double initial = 0, bulk = 0, air = 0, stained = 0, absorbed = 0, lost = 0, bled = 0;
    double total() const { return bulk + air + stained + absorbed + lost; }
};

class World {
public:
    explicit World(const Library& lib, Budget budget = budgetForTier(2), std::uint64_t seed = 1);

    const Library& library() const { return lib_; }
    const Budget& budget() const { return budget_; }
    Rng& rng() { return rng_; }

    // ---- Building ----
    int addCollider(Collider c);
    int addPlane(Role role, const std::string& name, const Vec3& point, const Vec3& normal, const SurfaceDef* s);
    int addBox(Role role, const std::string& name, const Vec3& centre, const Vec3& half, const Quat& rot,
               const SurfaceDef* s);
    int addCapsule(Role role, const std::string& name, const Vec3& a, const Vec3& b, double r, const SurfaceDef* s);
    // A vessel standing with its base centre at `base`, rotated by `rot`; returns its body.
    int addVessel(const VesselDef& def, const Vec3& base, const Quat& rot);
    // Fills a vessel (or tops it up) with `volume` m^3 of liquid, moving with it.
    double fillVessel(int body, const LiquidDef& liquid, double volume);
    void setLiquid(const LiquidDef& liquid);
    // A block of liquid (dam break, pours): fills an axis-aligned box.
    int addLiquidBox(const LiquidDef& liquid, const Vec3& lo, const Vec3& hi, const Vec3& velocity);
    void addDummy(const DummyDef& def, const Vec3& feet, double yaw);
    int shoot(const ProjectileDef& def, const Vec3& muzzle, const Vec3& dir, double speed, double delay = 0);
    // Releases single drops (drip tests, impact-angle cards).
    int drop(const Vec3& at, const Vec3& velocity, double diameter, const std::string& kind);

    // ---- Running ----
    int advance(double t, int maxSteps = 1000000);
    void step(double dt);
    double nextStep() const;
    double time() const { return time_; }
    bool quiet() const;
    int steps() const { return steps_; }

    // ---- Panes ----
    std::function<PaneResult(const PaneStrike&)> paneResolver;  // optional, synchronous
    bool deferPanes = false;  // wait for resolveStrike() (a host that runs Glass Dynamics itself)
    const std::vector<PaneStrike>& pendingStrikes() const { return pending_; }
    void resolveStrike(int index, const PaneResult& r);
    static PaneResult estimatePane(const PaneStrike& s, const Collider& pane, const ProjectileDef& p);
    void punchHole(int collider, double x, double y, double radius);  // pane frame from the centre

    // ---- State ----
    Fluid fluid;
    Spray spray;
    std::vector<Collider> colliders;
    std::vector<RigidBody> bodies;
    std::vector<Projectile> projectiles;
    std::vector<Wound> wounds;
    std::vector<Event> events;
    const DummyDef* dummy = nullptr;
    VolumeBook volumes() const;
    // Floor area the bulk liquid has wetted (m^2) and covers now.
    double wettedArea() const { return wetCells_ * wetCell_ * wetCell_; }
    double coveredArea() const;
    int floorCollider() const { return floor_; }
    const ShardPattern* shardsOf(int body) const;
    int breaks() const { return breaks_; }
    double strikeSpeed = 0.25;  // m/s, slower contacts are not tested for fracture
    void log(const std::string& text);

private:
    void stepProjectiles(double tEnd);
    void spatter(const Vec3& at, const Vec3& normal, const Vec3& axis, const SpatterModel& m, double energy,
                 double time, const std::string& kind, double speed);
    void breakVessel(int body, const Strike& s, const ImpactAssessment& a);
    void bleed(double tEnd);
    void updateWet();
    void syncColliders();

    const Library& lib_;
    Budget budget_;
    Rng rng_;
    double time_ = 0;
    int steps_ = 0;
    std::deque<VesselGeometry> geoms_;
    std::deque<ShardPattern> patterns_;
    std::vector<int> bodyPattern_;
    std::vector<PaneStrike> pending_;
    double initialVolume_ = 0, bledVolume_ = 0;
    int floor_ = -1;
    int breaks_ = 0;
    std::vector<unsigned char> wet_;
    int wetN_ = 400, wetCells_ = 0;
    double wetCell_ = 0.01;
    bool liquidSet_ = false;
};

}  // namespace lqd
