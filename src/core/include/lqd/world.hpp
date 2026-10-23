// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// A scene: static surfaces, glass bodies, the bulk liquid, droplets and
// stains, rounds in flight, impactors and the wounds they leave. advance() runs whole
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
#include "lqd/ragdoll.hpp"
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
    int spurts = 0;  // arterial pulses thrown so far (spatterBoost above 5)
    Attachment at;  // where it is on the dummy (it moves with the ragdoll)
    bool rupture = false;  // the skin burst under a blow (not a round's wound)
};

// Something heavy and blunt in motion (see ImpactorDef): its collider moves
// with it, the dummy's contacts slow it, its brakes stop it.
struct Impactor {
    const ImpactorDef* def = nullptr;
    int collider = -1;
    Vec3 vel;
    double speed0 = 0;          // m/s it set off at
    double contactAt = -1;      // when it first touched the dummy
    bool stopped = false;
};

// The dummy's energy (J): its ragdoll's book, and what the rounds left in its
// tissue. Everything it was given is now either moving it or absorbed.
struct DummyEnergy {
    EnergyBook book;
    double inTissue = 0;        // deposited by rounds crossing it
    double kinetic = 0, potential = 0;
    int bursts = 0;
    double burstVolume = 0;     // m^3 thrown by bursts
    double given() const { return book.delivered() + inTissue; }
    double absorbed() const { return book.absorbed() + inTissue; }
};

// A pane in pieces, as Glass Dynamics hands it over (see docs/integration.md).
struct PanePiece {
    std::vector<double> outline;  // x, y pairs about its centroid, m
    std::vector<std::vector<double>> holes;  // loops cut out of it (a hole punched through)
    double thickness = 0, radius = 0;
    bool chip = false, posed = false, moved = false;
    Vec3 pos, vel;                // world
    Quat rot;
    int state = 0;                // 0 held, 1 pending, 2 free, 3 resting, -1 retired
    int collider = -1;            // its slab once it is loose
};
struct PanePieces {
    int pane = -1;
    PaneMap map;
    std::vector<PanePiece> pieces;
    int loose = 0;                // slabs made
};
struct PieceStain {
    std::size_t stain = 0;
    int pane = -1, piece = -1;    // -1: waiting for its piece
    Vec3 p, n, d;                 // in the piece's frame
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
    // The dummy; as a ragdoll it stands on a release stand until a round hits it.
    void addDummy(const DummyDef& def, const Vec3& feet, double yaw, bool ragdollOn = true);
    int shoot(const ProjectileDef& def, const Vec3& muzzle, const Vec3& dir, double speed, double delay = 0);
    // An impactor centred at `centre`, setting off along `dir` at `speed`; returns its index.
    int addImpactor(const ImpactorDef& def, const Vec3& centre, const Vec3& dir, double speed);
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

    // ---- Panes in pieces ----
    // A host running Glass Dynamics hands over a broken pane: the map of which
    // piece each cell holds (Glass Dynamics' fragment map of the glass ply
    // nearest each face), each piece's outline, and every piece's pose as it
    // moves. All are in Glass Dynamics' pane frame: x from the left edge, y up
    // from the bottom edge, z into the glass from the front face (the pane's
    // +z face here). Drops land on the pieces, and their stains go where the
    // pieces go; loose pieces between the host's updates fly on ballistically.
    void setPaneMap(int pane, int nx, int ny, double cell, const std::vector<int>& front, const std::vector<int>& back,
                    bool laminated);
    void addPanePiece(int pane, int piece, const std::vector<double>& outline, double thickness, bool chip,
                      const std::vector<std::vector<double>>& holes = {});
    // rot is (w, x, y, z); state 0 held, 1 pending, 2 free, 3 resting, -1 retired.
    void setPanePiecePose(int pane, int piece, const Vec3& pos, const Quat& rot, const Vec3& vel, int state);
    const PanePieces* panePieces(int pane) const;
    const std::deque<PanePieces>& brokenPanes() const { return broken_; }
    const std::vector<PieceStain>& pieceStains() const { return pieceStains_; }
    // Stains that have moved since the last call (on the dummy, on pieces of glass).
    std::vector<std::size_t> takeMovedStains();
    int maxLoosePieces = 256;

    // ---- State ----
    Fluid fluid;
    Spray spray;
    std::vector<Collider> colliders;
    std::vector<RigidBody> bodies;
    std::vector<Projectile> projectiles;
    std::vector<Wound> wounds;
    std::vector<Impactor> impactors;
    std::vector<Event> events;
    const DummyDef* dummy = nullptr;
    Ragdoll ragdoll;
    bool releaseOnHit = true;
    // More than life (1 is as measured, up to 11): the spatter's volume scales
    // by it, its drops' size by its square root and their speed by its 0.15
    // power, the bleeding by its square root. Above 5 the wounds also spurt
    // like a cut artery in the films: three pulses a heartbeat apart, each a
    // jet of up to 10 ml at up to 7 m/s at 11. For films and games; the stain
    // analysis still works on what lands.
    double spatterBoost = 1;
    // A hit's push on the ragdoll. 0 is as in life: the momentum the round
    // loses, which barely moves a body (a 9 mm round stopping in 78 kg moves
    // it about 3 cm/s). Above 0, at least this fraction of the energy the
    // round deposits goes into the body's motion, at the point it strikes,
    // so a hit off the centre spins it (the cinema's knockback).
    double knockback = 0;
    VolumeBook volumes() const;
    DummyEnergy dummyEnergy() const;
    // Floor area the bulk liquid has wetted (m^2) and covers now.
    double wettedArea() const { return wetCells_ * wetCell_ * wetCell_; }
    double coveredArea() const;
    int floorCollider() const { return floor_; }
    const ShardPattern* shardsOf(int body) const;
    int breaks() const { return breaks_; }
    double strikeSpeed = 0.25;  // m/s, slower contacts are not tested for fracture
    void log(const std::string& text);
    // The collider of ragdoll segment s.
    int dummyCollider(int s) const { return s >= 0 && s < static_cast<int>(dummyColliders_.size()) ? dummyColliders_[s] : -1; }
    // Stains riding on the dummy: (stain index, where on it).
    const std::vector<std::pair<std::size_t, Attachment>>& carried() const { return carried_; }

private:
    void stepProjectiles(double tEnd);
    void spatter(const Vec3& at, const Vec3& normal, const Vec3& axis, const SpatterModel& m, double energy, bool boosted,
                 double time, const std::string& kind, double speed, const Vec3& carrier = {},
                 double atLeast = 48);
    void breakVessel(int body, const Strike& s, const ImpactAssessment& a);
    void bleed(double tEnd);
    void stepImpactors(double dt);
    void carryImpactorStains();
    void ruptures(double dt);
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
    std::vector<std::pair<std::size_t, Attachment>> carried_;  // stains riding on the dummy
    std::size_t stainsSeen_ = 0;
    std::vector<int> dummyColliders_;
    // Per ragdoll segment: what it has absorbed in the blow under way, whether
    // it has passed its rupture threshold (since when), where, and its bursts.
    struct Bruise {
        double energy = 0, best = 0, speed = 0, since = -1;
        Attachment at;
        int bursts = 0;
    };
    std::vector<Bruise> bruises_;
    double tissue_ = 0, burstVolume_ = 0, speedCap_ = 15;
    std::vector<std::pair<std::size_t, Vec3>> onImpactor_;  // stains riding on an impactor: offset from its centre
    std::size_t impactorSeen_ = 0;
    void carryStains();
    void carryPieceStains();
    PanePieces& piecesOf(int pane);
    void markMoved(std::size_t stain);
    std::deque<PanePieces> broken_;
    std::vector<PieceStain> pieceStains_;
    std::vector<std::size_t> moved_;
    std::vector<unsigned char> movedFlag_;
    bool piecesDirty_ = false;
    std::size_t pieceSeen_ = 0;
    double initialVolume_ = 0, bledVolume_ = 0;
    int floor_ = -1;
    int breaks_ = 0;
    std::vector<unsigned char> wet_;
    int wetN_ = 400, wetCells_ = 0;
    double wetCell_ = 0.01;
    bool liquidSet_ = false;
};

}  // namespace lqd
