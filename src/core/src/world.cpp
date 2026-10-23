// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/world.hpp"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>

namespace lqd {

Budget budgetForTier(int tier) {
    switch (std::clamp(tier, 0, 4)) {
        case 0: return {"potato", 0.008, 3, 3000, 40};
        case 1: return {"low", 0.0065, 3, 6000, 48};
        case 2: return {"mid", 0.005, 4, 12000, 64};
        case 3: return {"high", 0.004, 4, 20000, 80};
        default: return {"ultra", 0.003, 5, 40000, 96};
    }
}

namespace {

// Bulk liquid leaving faster than this relative to the rest tears into drops
// (estimated: splash sheets and ligaments break up within centimetres), of this
// median diameter (estimated: the ligament drops of mm-scale splashes, 1-4 mm).
constexpr double kFragmentSpeed = 2.0;
constexpr double kFragmentMedian = 2.2e-3;
// A wound oozes at the dummy's bleed rate decaying over this (estimated: a
// stage-blood squib empties in a second or two), and stops after kBleedFor.
constexpr double kBleedDecay = 1.5;
constexpr double kBleedFor = 6.0;
// A blow is what a segment absorbs within about this long (estimated: a car
// bumper's contact with a leg lasts 10-20 ms, a fall's landing about as long);
// a segment past its rupture threshold bursts half of it later, when the blow
// is mostly over, and bursts at most kMaxBursts times.
constexpr double kBlowTime = 0.03;
constexpr int kMaxBursts = 2;

std::string fmt(const char* f, ...) {
    char buf[512];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof buf, f, ap);
    va_end(ap);
    return buf;
}

}  // namespace

World::World(const Library& lib, Budget budget, std::uint64_t seed) : lib_(lib), budget_(budget), rng_(seed) {
    wet_.assign(static_cast<std::size_t>(wetN_) * wetN_, 0);
    spray.maxDrops = budget.maxDrops;
}

void World::log(const std::string& text) { events.push_back({time_, text}); }

// ---- Building ----------------------------------------------------------------

int World::addCollider(Collider c) {
    c.id = static_cast<int>(colliders.size());
    switch (c.kind) {
        case ShapeKind::Plane: c.boundRadius = 1e30; break;
        case ShapeKind::Box: c.boundRadius = length(c.half); break;
        case ShapeKind::Capsule:
            c.pos = (c.a + c.b) * 0.5;
            c.boundRadius = 0.5 * length(c.b - c.a) + c.radius;
            break;
        case ShapeKind::Shell:
            c.boundRadius = c.shell ? std::hypot(c.shell->height(), c.shell->maxRadius()) : 1;
            break;
        case ShapeKind::Slab: {
            double r = 0;
            for (std::size_t i = 0; i + 1 < c.outline.size(); i += 2) r = std::max(r, std::hypot(c.outline[i], c.outline[i + 1]));
            c.boundRadius = std::hypot(r, c.half.z);
            break;
        }
    }
    if (c.surface) c.color = c.surface->color;
    if (c.role == Role::Floor && floor_ < 0) floor_ = c.id;
    colliders.push_back(c);
    return c.id;
}

int World::addPlane(Role role, const std::string& name, const Vec3& point, const Vec3& normal, const SurfaceDef* s) {
    Collider c;
    c.name = name;
    c.role = role;
    c.kind = ShapeKind::Plane;
    c.pos = point;
    const Vec3 n = normalize(normal);
    const Vec3 up{0, 1, 0};
    const Vec3 ax = cross(up, n);
    const double ang = std::acos(clamp(dot(up, n), -1.0, 1.0));
    c.rot = length(ax) > 1e-9 ? Quat::axisAngle(ax, ang) : (dot(up, n) > 0 ? Quat{} : Quat::axisAngle({1, 0, 0}, kPi));
    c.surface = s;
    return addCollider(c);
}

int World::addBox(Role role, const std::string& name, const Vec3& centre, const Vec3& half, const Quat& rot,
                  const SurfaceDef* s) {
    Collider c;
    c.name = name;
    c.role = role;
    c.kind = ShapeKind::Box;
    c.pos = centre;
    c.half = half;
    c.rot = rot;
    c.surface = s;
    return addCollider(c);
}

int World::addCapsule(Role role, const std::string& name, const Vec3& a, const Vec3& b, double r,
                      const SurfaceDef* s) {
    Collider c;
    c.name = name;
    c.role = role;
    c.kind = ShapeKind::Capsule;
    c.a = a;
    c.b = b;
    c.radius = r;
    c.surface = s;
    return addCollider(c);
}

int World::addVessel(const VesselDef& def, const Vec3& base, const Quat& rot) {
    geoms_.emplace_back(def, budget_.columns);
    const VesselGeometry& g = geoms_.back();
    RigidBody b;
    b.id = static_cast<int>(bodies.size());
    b.name = def.name;
    b.geom = &g;
    b.mass = g.mass();
    b.inertia = g.inertia();
    b.comLocal = g.centreOfMass();
    b.q = rot;
    b.x = base + rot.rotate(b.comLocal);
    b.samples = g.samples();
    b.breakable = true;
    b.born = time_;
    Collider c;
    c.name = def.name;
    c.role = Role::Vessel;
    c.kind = ShapeKind::Shell;
    c.shell = &g;
    c.comLocal = b.comLocal;
    c.body = b.id;
    c.pos = b.x;
    c.rot = b.q;
    c.stains = false;
    c.surface = lib_.surface("glassware");
    b.collider = addCollider(c);
    bodies.push_back(b);
    bodyPattern_.push_back(-1);
    char buf[256];
    std::snprintf(buf, sizeof buf, "%s: %.0f g of glass, cavity %.0f ml to the brim", def.name.c_str(),
                  g.mass() * 1e3, g.cavityCapacity() * 1e6);
    log(buf);
    return b.id;
}

void World::setLiquid(const LiquidDef& liquid) {
    if (liquidSet_) return;
    FluidSettings fs;
    fs.spacing = budget_.spacing;
    fs.iterations = budget_.iterations;
    fluid.configure(liquid, fs);
    spray.configure(liquid);
    spray.particleVolume = fluid.particleMass() / fluid.restDensity();
    liquidSet_ = true;
}

double World::fillVessel(int body, const LiquidDef& liquid, double volume) {
    setLiquid(liquid);
    const RigidBody& B = bodies[body];
    const VesselGeometry& g = *B.geom;
    const double h0 = fluid.settings().spacing;
    const double R = g.maxRadius();
    const double vp = fluid.particleMass() / fluid.restDensity();
    // The lattice points inside the cavity below `top`.
    auto lattice = [&](double top, std::vector<Vec3>* out) {
        int n = 0;
        for (double y = 0.5 * h0; y < top; y += h0)
            for (double xx = -R; xx <= R; xx += h0)
                for (double zz = -R; zz <= R; zz += h0) {
                    if (!g.inCavity(std::hypot(xx, zz), y, 0.5 * h0)) continue;
                    if (out) out->push_back({xx, y, zz});
                    ++n;
                }
        return n;
    };
    // Whole lattice layers up to the volume asked for, then an evenly thinned
    // top layer for the rest, so the particles hold the volume asked for.
    std::vector<Vec3> all;
    lattice(g.rimHeight(), &all);
    const std::size_t want = std::min(all.size(), static_cast<std::size_t>(std::llround(volume / vp)));
    std::vector<Vec3> pts;
    std::size_t i = 0;
    while (i < all.size() && pts.size() < want) {
        std::size_t j = i;
        while (j < all.size() && all[j].y == all[i].y) ++j;
        const std::size_t layer = j - i, need = want - pts.size();
        if (layer <= need) pts.insert(pts.end(), all.begin() + static_cast<long>(i), all.begin() + static_cast<long>(j));
        else
            for (std::size_t k = 0; k < need; ++k) pts.push_back(all[i + (k * layer + layer / 2) / need]);
        i = j;
    }
    const int n = static_cast<int>(pts.size());
    const double top = pts.empty() ? 0 : pts.back().y;
    for (const auto& q : pts) {
        const Vec3 pw = B.fromVessel(q);
        fluid.add(pw, B.velocityAt(pw));
    }
    const double got = n * fluid.particleMass() / fluid.restDensity();
    initialVolume_ += got;
    char buf[256];
    std::snprintf(buf, sizeof buf, "filled with %s: %d particles at %.1f mm, %.0f ml (asked %.0f ml, level %.0f mm)",
                  liquid.name.c_str(), n, h0 * 1e3, got * 1e6, volume * 1e6, top * 1e3);
    log(buf);
    return got;
}

int World::addLiquidBox(const LiquidDef& liquid, const Vec3& lo, const Vec3& hi, const Vec3& velocity) {
    setLiquid(liquid);
    const double h0 = fluid.settings().spacing;
    int n = 0;
    for (double y = lo.y + 0.5 * h0; y < hi.y; y += h0)
        for (double x = lo.x + 0.5 * h0; x < hi.x; x += h0)
            for (double z = lo.z + 0.5 * h0; z < hi.z; z += h0) {
                fluid.add({x, y, z}, velocity);
                ++n;
            }
    initialVolume_ += n * fluid.particleMass() / fluid.restDensity();
    return n;
}

void World::addDummy(const DummyDef& def, const Vec3& feet, double yaw, bool ragdollOn) {
    dummy = &def;
    const Quat q = Quat::axisAngle({0, 1, 0}, yaw);
    const SurfaceDef* skin = lib_.surface("dummy_vinyl");
    std::vector<int> ids;
    for (const auto& s : def.segments) {
        // Stains on the dummy are kept too: spatter on clothing and shoes is evidence.
        ids.push_back(addCapsule(Role::Dummy, s.name, feet + q.rotate(s.a), feet + q.rotate(s.b), s.radius, skin));
    }
    dummyColliders_ = ids;
    if (ragdollOn) {
        for (int id : ids) colliders[id].movable = true;
        ragdoll.build(def, ids, colliders);
        char buf[160];
        std::snprintf(buf, sizeof buf, "%s: %.0f kg ragdoll of %zu capsules on a release stand", def.name.c_str(),
                      ragdoll.mass(), ids.size());
        log(buf);
    }
    if (const LiquidDef* l = lib_.liquid(def.liquid)) {
        if (!liquidSet_) setLiquid(*l);
    }
}

int World::shoot(const ProjectileDef& def, const Vec3& muzzle, const Vec3& dir, double speed, double delay) {
    // A shot load is its pellets, each its own round: scattered about the aim
    // (a Gaussian of rms angle `spread`) and strung out along the path.
    const Vec3 aim = normalize(dir);
    const Vec3 u = perpendicular(aim), v = cross(aim, u);
    int first = -1;
    for (int k = 0; k < std::max(1, def.pellets); ++k) {
        Projectile p;
        p.id = static_cast<int>(projectiles.size());
        p.def = &def;
        p.pos = muzzle;
        Vec3 d = aim;
        double lag = 0;
        if (def.pellets > 1) {
            const double s = def.spread / std::sqrt(2.0);
            d = normalize(aim + u * (s * rng_.normal()) + v * (s * rng_.normal()));
            lag = def.shotString * rng_.uniform() / std::max(speed, 1.0);
        }
        p.dir = d;
        p.speed = speed;
        p.launch = time_ + delay + lag;
        p.time = p.launch;
        p.path.push_back(muzzle);
        projectiles.push_back(p);
        if (first < 0) first = p.id;
    }
    return first;
}

int World::addImpactor(const ImpactorDef& def, const Vec3& centre, const Vec3& dir, double speed) {
    const Vec3 d = normalize(dir, {0, 0, -1});
    // Its frame: z along its travel (level), y up, x across.
    const Quat q = Quat::axisAngle({0, 1, 0}, std::atan2(d.x, d.z));
    const SurfaceDef* s = def.surface.empty() ? nullptr : lib_.surface(def.surface);
    int id = -1;
    if (def.shape == ImpactorShape::Box) {
        id = addBox(Role::Impactor, def.name, centre, def.half, q, s);
    } else {
        // A sphere is a capsule too short to matter.
        const Vec3 h = q.rotate({def.shape == ImpactorShape::Capsule ? 0.5 * def.length : 1e-4, 0, 0});
        id = addCapsule(Role::Impactor, def.name, centre - h, centre + h, def.radius, s);
    }
    Collider& c = colliders[id];
    c.movable = true;
    c.vel = d * speed;
    Impactor im;
    im.def = &def;
    im.collider = id;
    im.vel = c.vel;
    im.speed0 = speed;
    impactors.push_back(im);
    // Nothing it strikes outruns it by much.
    speedCap_ = std::max(speedCap_, 1.3 * speed);
    ragdoll.setSpeedCap(speedCap_);
    log(fmt("%s sets off at %.1f m/s (%.0f km/h): %.0f kg carrying %.1f kJ", def.name.c_str(), speed, speed * 3.6,
            def.mass, 0.5 * def.mass * speed * speed * 1e-3));
    return static_cast<int>(impactors.size()) - 1;
}

int World::drop(const Vec3& at, const Vec3& velocity, double diameter, const std::string& kind) {
    const int src = spray.addSource(kind, at, time_);
    Droplet d;
    d.x = at;
    d.v = velocity;
    d.d = diameter;
    d.t = time_;
    d.source = src;
    spray.emit(d);
    initialVolume_ += kPi / 6 * diameter * diameter * diameter;
    return src;
}

const ShardPattern* World::shardsOf(int body) const {
    if (body < 0 || body >= static_cast<int>(bodyPattern_.size()) || bodyPattern_[body] < 0) return nullptr;
    return &patterns_[bodyPattern_[body]];
}

// ---- Running -------------------------------------------------------------------

double World::nextStep() const {
    double dt = 0.002;
    if (fluid.size() > 0) dt = std::min(dt, fluid.maxStep());
    for (const auto& p : projectiles)
        if (p.state == Projectile::State::Flying || p.state == Projectile::State::AtPane ||
            (p.state == Projectile::State::Waiting && p.launch < time_ + 0.002))
            dt = std::min(dt, 5e-5);
    double newest = -1;
    for (const auto& s : spray.sources) newest = std::max(newest, s.time);
    if (newest >= 0 && time_ - newest < 0.02 && !spray.drops.empty()) dt = std::min(dt, 2e-4);
    for (const auto& im : impactors) {
        const double sp = length(im.vel);
        if (im.stopped || sp < 0.01) continue;
        // Near the dummy, a centimetre a step: a thin limb must not slip through it.
        const Collider& c = colliders[im.collider];
        for (int id : dummyColliders_)
            if (c.sdf(colliders[id].pos) < colliders[id].boundRadius + 0.5) {
                dt = std::min(dt, 0.01 / sp);
                break;
            }
    }
    return std::max(dt, 1e-6);
}

int World::advance(double t, int maxSteps) {
    int n = 0;
    while (time_ < t && n < maxSteps) {
        step(nextStep());
        ++n;
    }
    return n;
}

void World::syncColliders() {
    for (auto& c : colliders)
        if (c.body >= 0) {
            c.pos = bodies[c.body].x;
            c.rot = bodies[c.body].q;
        }
}

void World::step(double dt) {
    const double tEnd = time_ + dt;
    stepProjectiles(tEnd);

    std::vector<Strike> strikes;
    for (auto& b : bodies)
        if (!b.broken) stepBody(b, colliders, dt, strikeSpeed, &strikes);
    for (const auto& s : strikes) {
        RigidBody& b = bodies[s.body];
        if (!b.breakable || b.broken || b.shard >= 0 || s.sample < 0) continue;
        const Collider& on = colliders[s.collider];
        const SurfaceDef* surf = on.surface ? on.surface : lib_.surface("ceramic_tile");
        if (!surf) continue;
        const double meff = b.effectiveMass(s.point, s.normal);
        const WallSample& ws = b.samples[s.sample];
        const double cn = std::abs(dot(b.q.unrotate(s.normal), ws.wallNormal));
        const ImpactAssessment a = assessImpact(b.geom->def(), ws, *surf, meff, s.speed, cn);
        const double u = rng_.uniform();
        char buf[400];
        std::snprintf(buf, sizeof buf,
                      "%s struck %s at %.2f m/s on its %s: effective mass %.0f g, peak force %.0f N, wall stress %.0f "
                      "MPa, P(break) %.2f -> %s",
                      b.name.c_str(), on.name.c_str(), s.speed, a.edge ? "edge" : "wall", meff * 1e3, a.force,
                      a.stress * 1e-6, a.probability, u < a.probability ? "breaks" : "holds");
        log(buf);
        if (u < a.probability) breakVessel(s.body, s, a);
    }
    // Loose pieces of glass fly on between the host's updates.
    for (auto& bp : broken_)
        for (auto& pc : bp.pieces)
            if (pc.posed && pc.state == 2) {
                pc.vel += Vec3{0, -kGravity, 0} * dt;
                pc.pos += pc.vel * dt;
                pc.moved = true;
                if (pc.collider >= 0) {
                    colliders[pc.collider].pos = pc.pos;
                    colliders[pc.collider].vel = pc.vel;
                }
            }
    if (!impactors.empty()) stepImpactors(dt);
    if (ragdoll.active()) {
        std::vector<Pusher> push;
        for (const auto& im : impactors)
            push.push_back({im.collider, im.stopped ? 1e12 : im.def->mass, im.stopped ? Vec3{} : im.vel, false});
        if (!ragdoll.released())
            for (const auto& im : impactors) {
                const Collider& c = colliders[im.collider];
                for (int id : dummyColliders_) {
                    const Collider& dc = colliders[id];
                    if (std::min({c.sdf(dc.a), c.sdf(dc.b), c.sdf(dc.pos)}) > dc.radius) continue;
                    ragdoll.release(time_);
                    log("the stand lets go: the dummy is a ragdoll");
                    break;
                }
                if (ragdoll.released()) break;
            }
        ragdoll.step(dt, colliders, push.empty() ? nullptr : &push);
        for (std::size_t k = 0; k < impactors.size(); ++k) {
            Impactor& im = impactors[k];
            if (im.stopped) continue;
            im.vel = push[k].vel;
            colliders[im.collider].vel = im.vel;
            if (im.contactAt >= 0) continue;
            const auto& H = ragdoll.hits();
            for (int s = 0; s < static_cast<int>(H.size()); ++s)
                if (H[s].collider == im.collider) {
                    im.contactAt = time_;
                    log(fmt("%s strikes the dummy's %s at %.1f m/s", im.def->name.c_str(), ragdoll.segmentName(s).c_str(),
                            length(im.vel)));
                    break;
                }
        }
        ragdoll.pose(colliders);
        for (auto& w : wounds) ragdoll.place(w.at, w.pos, w.normal);
        ruptures(dt);
    }
    syncColliders();

    std::vector<int> toSpray;
    fluid.step(dt, colliders, bodies, &toSpray);
    std::vector<Droplet> rejoin;
    spray.step(tEnd, colliders, &fluid, &rejoin, rng_);
    if (!toSpray.empty()) {
        const double d = std::cbrt(6 * fluid.particleMass() / (kPi * fluid.restDensity()));
        // Fast relative to the liquid it left: wine carried along in a falling
        // glass is not a splash.
        Vec3 mean;
        for (const auto& u : fluid.v) mean += u;
        mean = mean / static_cast<double>(std::max<std::size_t>(fluid.size(), 1));
        for (int i : toSpray) {
            Droplet dr;
            dr.x = fluid.x[i];
            dr.v = fluid.v[i];
            dr.d = d;
            dr.t = tEnd;
            dr.fromBulk = true;
            // Liquid thrown off fast is a sheet or a ligament, not a ball: it
            // tears into drops of a few millimetres (the size of the ligaments
            // a splash throws), log-normally spread. Slow liquid leaves as a
            // blob that joins the bulk again where it lands.
            const double sp = length(dr.v - mean);
            const double vol = kPi / 6 * d * d * d;
            if (sp > kFragmentSpeed && static_cast<int>(spray.drops.size()) + 40 < spray.maxDrops) {
                double left = vol;
                const Vec3 u = normalize(dr.v);
                while (left > 0) {
                    double dd = clamp(rng_.lognormal(kFragmentMedian, 0.5), 0.6e-3, 0.8 * d);
                    double v = kPi / 6 * dd * dd * dd;
                    if (v > left || left - v < kPi / 6 * 0.2e-9) {
                        v = left;
                        dd = std::cbrt(6 * v / kPi);
                    }
                    Droplet f = dr;
                    f.fromBulk = false;
                    f.d = dd;
                    f.x = dr.x + rng_.cone(u, kPi) * (0.4 * d);
                    f.v = dr.v * (1 + 0.08 * rng_.normal()) + rng_.cone(u, kPi) * (0.1 * sp);
                    spray.drops.push_back(f);
                    left -= v;
                }
                continue;
            }
            if (static_cast<int>(spray.drops.size()) < spray.maxDrops) spray.drops.push_back(dr);
            else rejoin.push_back(dr);  // no room in the air: it stays liquid
        }
        fluid.remove(toSpray);
    }
    for (const auto& r : rejoin) fluid.add(r.x, r.v);
    // Drops that land in the bulk liquid join it a particle at a time.
    if (fluid.liquid()) {
        const double vp = fluid.particleMass() / fluid.restDensity();
        while (spray.toBulk >= vp) {
            const double h = fluid.settings().spacing;
            fluid.add(spray.toBulkAt + Vec3{rng_.uniform(-0.4, 0.4), rng_.uniform(0, 0.4), rng_.uniform(-0.4, 0.4)} * h, spray.toBulkV);
            spray.toBulk -= vp;
        }
    }
    // New puddle liquid, scattered a little in the surface so it never sits
    // exactly on another particle.
    for (const auto& [at, n] : spray.pooled) {
        const Vec3 u = perpendicular(n), w = cross(n, u);
        const double h = fluid.settings().spacing;
        fluid.add(at + n * (0.5 * h) + u * rng_.uniform(-0.5 * h, 0.5 * h) + w * rng_.uniform(-0.5 * h, 0.5 * h), {});
    }
    spray.pooled.clear();
    carryStains();
    carryPieceStains();
    carryImpactorStains();
    bleed(tEnd);
    time_ = tEnd;
    ++steps_;
    updateWet();
}

bool World::quiet() const {
    for (const auto& p : projectiles)
        if (p.state != Projectile::State::Stopped) return false;
    if (!spray.drops.empty()) return false;
    for (const auto& b : bodies)
        if (!b.broken && !b.sleeping) return false;
    if (fluid.size() > 0 && fluid.maxSpeed() > 0.03) return false;
    for (const auto& w : wounds)
        if (dummy && dummy->bleedRate > 0 && time_ - w.time < kBleedFor) return false;
    if (ragdoll.moving()) return false;
    for (const auto& im : impactors)
        if (!im.stopped && length2(im.vel) > 1e-4) return false;
    return true;
}

DummyEnergy World::dummyEnergy() const {
    DummyEnergy e;
    e.inTissue = tissue_;
    e.burstVolume = burstVolume_;
    for (const auto& b : bruises_) e.bursts += b.bursts;
    if (ragdoll.active()) {
        e.book = ragdoll.energy();
        e.kinetic = ragdoll.released() ? ragdoll.kinetic() : 0.0;
        e.potential = ragdoll.potential();
    }
    return e;
}

void World::stepImpactors(double dt) {
    for (auto& im : impactors) {
        Collider& c = colliders[im.collider];
        if (im.stopped) {
            c.vel = {};
            continue;
        }
        if (im.contactAt >= 0 && im.def->brake > 0 && time_ - im.contactAt >= im.def->brakeDelay) {
            const double sp = length(im.vel), dv = im.def->brake * dt;
            im.vel = sp > dv ? im.vel * ((sp - dv) / sp) : Vec3{};
            if (sp <= dv) {
                im.stopped = true;
                c.vel = {};
                log(fmt("%s comes to a stop", im.def->name.c_str()));
                continue;
            }
        }
        const Vec3 move = im.vel * dt;
        // It stops a body's thickness short of a wall in its way, as if it had
        // hit the wall (what a dummy pinned there absorbs would depend on how it
        // crushes, which is not modelled). It is held up by a road or a chain.
        for (const auto& o : colliders) {
            if (o.kind != ShapeKind::Plane || !o.enabled || o.role == Role::Floor) continue;
            Vec3 n;
            const double d = o.sdf(c.pos + move, &n);
            double ext = c.radius + 0.5 * std::abs(dot(n, c.b - c.a));
            if (c.kind == ShapeKind::Box) {
                ext = 0;
                for (int k = 0; k < 3; ++k) {
                    Vec3 e;
                    e[k] = 1;
                    ext += std::abs(dot(n, c.rot.rotate(e))) * c.half[k];
                }
            }
            if (d - ext < 0.3 && dot(im.vel, n) < 0) {
                im.stopped = true;
                log(fmt("%s stops against the %s", im.def->name.c_str(), o.name.c_str()));
                break;
            }
        }
        if (im.stopped) {
            im.vel = c.vel = {};
            continue;
        }
        c.pos += move;
        c.a += move;
        c.b += move;
        c.vel = im.vel;
    }
}

void World::carryImpactorStains() {
    if (impactors.empty()) return;
    for (; impactorSeen_ < spray.stains.size(); ++impactorSeen_) {
        const Stain& st = spray.stains[impactorSeen_];
        if (st.collider >= 0 && colliders[st.collider].role == Role::Impactor)
            onImpactor_.push_back({impactorSeen_, st.pos - colliders[st.collider].pos});
    }
    for (const auto& [k, off] : onImpactor_) {
        Stain& st = spray.stains[k];
        const Collider& c = colliders[st.collider];
        const Vec3 p = c.pos + off;
        if (length2(p - st.pos) < 1e-12) continue;
        st.pos = p;
        markMoved(k);
    }
}

void World::ruptures(double dt) {
    if (!dummy || dummy->ruptureEnergy <= 0 || dummy->burst.volumePerJoule <= 0 || !ragdoll.released()) return;
    const auto& H = ragdoll.hits();
    if (bruises_.size() != H.size()) bruises_.assign(H.size(), {});
    const double boost = std::max(spatterBoost, 0.1), fade = std::exp(-dt / kBlowTime);
    for (int s = 0; s < static_cast<int>(H.size()); ++s) {
        Bruise& b = bruises_[s];
        const SegmentHit& h = H[s];
        // What it absorbs in one blow: older contacts fade (a body lying on the
        // floor is not slowly bursting), unless it is already bursting.
        if (b.since < 0) {
            b.energy *= fade;
            b.best *= fade;
            b.speed *= fade;
        }
        b.energy += h.energy;
        b.speed = std::max(b.speed, h.speed);
        if (h.best > b.best) {
            b.best = h.best;
            b.at = ragdoll.attach(s, h.point, -h.normal);
        }
        // A bag of liquid bursts on the pressure of a fast blow, not on a slow
        // squeeze: it must be struck fast enough, and take enough. More than
        // life, it bursts sooner: at 11 under a third of the energy, at 55% of the speed.
        const double threshold = dummy->ruptureEnergy * ragdoll.segmentMass(s) / std::sqrt(boost);
        const double fast = dummy->ruptureSpeed / std::pow(boost, 0.25);
        if (b.since < 0 && b.bursts < kMaxBursts && b.energy >= threshold && b.speed >= fast) b.since = time_;
        if (b.since < 0 || time_ - b.since < 0.5 * kBlowTime) continue;
        // The skin splits where the blow went in, and the blow squeezes the
        // liquid out both ways round the segment.
        Vec3 at, nrm;
        ragdoll.place(b.at, at, nrm);
        const Vec3 e = normalize(ragdoll.b(s) - ragdoll.a(s), {0, 1, 0});
        Vec3 side = cross(nrm, e);
        side = length2(side) > 1e-6 ? normalize(side) : perpendicular(nrm);
        const double perJoule = dummy->burst.volumePerJoule * boost;
        const double E = dummy->maxBurst > 0 ? std::min(b.energy, dummy->maxBurst / perJoule) : b.energy;
        log(fmt("the %s bursts: %.0f J absorbed in the blow at %.1f m/s (it bursts above %.0f J)",
                ragdoll.segmentName(s).c_str(), b.energy, b.speed, threshold));
        const double before = initialVolume_;
        const Vec3 carrier = ragdoll.velocityOf(b.at);
        for (int k = -1; k <= 1; k += 2)
            spatter(at, nrm, normalize(side * k + nrm * 0.6), dummy->burst, 0.5 * E, true, time_, "burst",
                    std::max(b.speed, 2.0), carrier);
        burstVolume_ += initialVolume_ - before;
        Wound wd{at, nrm, false, time_, 0, 0, static_cast<int>(spray.sources.size()) - 1, 0, b.at};
        wd.rupture = true;
        wounds.push_back(wd);
        ++b.bursts;
        b.since = -1;
        b.energy = b.best = b.speed = 0;
    }
}

VolumeBook World::volumes() const {
    VolumeBook v;
    v.initial = initialVolume_;
    v.bulk = fluid.volume() + spray.toBulk;
    v.air = spray.airVolume();
    v.stained = spray.stainedVolume;
    v.absorbed = spray.absorbedVolume;
    v.lost = spray.lostVolume;
    v.bled = bledVolume_;
    return v;
}

void World::updateWet() {
    if (floor_ < 0 || fluid.size() == 0) return;
    const Collider& f = colliders[floor_];
    const double reach = 0.8 * fluid.settings().spacing;
    for (const auto& p : fluid.x) {
        if (f.sdf(p) > reach) continue;
        const int i = static_cast<int>(std::floor(p.x / wetCell_)) + wetN_ / 2;
        const int k = static_cast<int>(std::floor(p.z / wetCell_)) + wetN_ / 2;
        if (i < 0 || k < 0 || i >= wetN_ || k >= wetN_) continue;
        unsigned char& c = wet_[static_cast<std::size_t>(k) * wetN_ + i];
        if (!c) {
            c = 1;
            ++wetCells_;
        }
    }
}

double World::coveredArea() const {
    if (floor_ < 0) return 0;
    const Collider& f = colliders[floor_];
    const double reach = 0.8 * fluid.settings().spacing;
    std::vector<unsigned char> now(wet_.size(), 0);
    int n = 0;
    for (const auto& p : fluid.x) {
        if (f.sdf(p) > reach) continue;
        const int i = static_cast<int>(std::floor(p.x / wetCell_)) + wetN_ / 2;
        const int k = static_cast<int>(std::floor(p.z / wetCell_)) + wetN_ / 2;
        if (i < 0 || k < 0 || i >= wetN_ || k >= wetN_) continue;
        unsigned char& c = now[static_cast<std::size_t>(k) * wetN_ + i];
        if (!c) {
            c = 1;
            ++n;
        }
    }
    return n * wetCell_ * wetCell_;
}

// ---- Breaking glass ----------------------------------------------------------------

void World::breakVessel(int bi, const Strike& s, const ImpactAssessment& a) {
    const RigidBody B = bodies[bi];  // copy: the body list grows below
    const VesselGeometry& g = *B.geom;
    const WallSample ws = B.samples[s.sample];
    const double energy = 0.5 * a.effectiveMass * a.speed * a.speed;
    const int count = std::clamp(static_cast<int>(std::round(6 + 10 * std::sqrt(energy / 0.1))), 4, 64);
    const double spread = 0.12 * g.length();
    patterns_.push_back(makeShards(g, ws.s, ws.theta, count, spread, rng_));
    const int pi = static_cast<int>(patterns_.size()) - 1;
    const ShardPattern& P = patterns_[pi];
    const double rho = g.def().glassDensity;
    std::vector<std::vector<int>> cells(P.count);
    for (int c = 0; c < g.cells(); ++c) cells[P.cellShard[c]].push_back(c);
    int made = 0;
    for (int k = 0; k < P.count; ++k) {
        double M = 0;
        Vec3 com;
        for (int c : cells[k]) {
            const double m = rho * g.cellVolume(c);
            M += m;
            com += g.cellPoint(c) * m;
        }
        if (M <= 0) continue;
        com = com / M;
        Mat3 I;
        for (int c : cells[k]) {
            const double m = rho * g.cellVolume(c);
            const Vec3 r = g.cellPoint(c) - com;
            I = I + (Mat3::identity() * length2(r) - outer(r, r)) * m;
            const double cell = std::sqrt(std::max(g.cellVolume(c), 1e-15) / std::max(g.def().profile[0].t, 1e-4));
            I = I + Mat3::identity() * (m * cell * cell / 6);
        }
        RigidBody sb;
        sb.id = static_cast<int>(bodies.size());
        sb.name = B.name + " shard " + std::to_string(k + 1);
        sb.mass = M;
        sb.inertia = I;
        sb.geom = &g;
        sb.shard = k;
        sb.comLocal = com;
        sb.q = B.q;
        sb.x = B.fromVessel(com);
        sb.w = B.w;
        const Vec3 away = normalize(sb.x - s.point, s.normal);
        sb.v = B.velocityAt(sb.x) + away * (0.12 * a.speed) + rng_.cone({0, 1, 0}, kPi) * (0.04 * a.speed);
        for (const auto& w : g.samples())
            if (P.shardAt(g, w.s, w.theta) == k) sb.samples.push_back(w);
        const int stride = std::max(1, static_cast<int>(cells[k].size()) / 24);
        for (std::size_t q = 0; q < cells[k].size(); q += stride) {
            WallSample w;
            w.p = g.cellPoint(cells[k][q]);
            sb.samples.push_back(w);
        }
        sb.restitution = 0.6;
        sb.born = time_;
        Collider c;
        c.name = sb.name;
        c.role = Role::Shard;
        c.kind = ShapeKind::Shell;
        c.shell = &g;
        c.shards = &P;
        c.shard = k;
        c.comLocal = com;
        c.body = sb.id;
        c.pos = sb.x;
        c.rot = sb.q;
        c.stains = false;
        c.surface = lib_.surface("glassware");
        sb.collider = addCollider(c);
        colliders[sb.collider].boundRadius = 0;
        for (int cc : cells[k]) colliders[sb.collider].boundRadius = std::max(colliders[sb.collider].boundRadius, length(g.cellPoint(cc) - com) + 0.01);
        bodies.push_back(sb);
        bodyPattern_.push_back(pi);
        ++made;
    }
    bodies[bi].broken = true;
    bodies[bi].brokenAt = time_;
    bodies[bi].samples.clear();
    bodyPattern_[bi] = pi;
    colliders[B.collider].enabled = false;
    ++breaks_;
    log(fmt("  shattered into %d pieces (%.2f J at the strike)", made, energy));
}

// ---- Rounds ------------------------------------------------------------------------

PaneResult World::estimatePane(const PaneStrike& s, const Collider& pane, const ProjectileDef& p) {
    PaneResult r;
    const double tmm = 2 * pane.half.z * 1e3;
    const double absorb = 18.0 * tmm;  // J per mm of glass (estimated, see docs/physics.md)
    const double ke = 0.5 * p.mass * s.speed * s.speed;
    if (ke <= absorb) {
        r.penetrated = false;
        r.speedOut = 0;
        return r;
    }
    r.speedOut = std::sqrt(2 * (ke - absorb) / p.mass);
    r.holeRadius = 1.5 * p.diameter;
    return r;
}

void World::punchHole(int collider, double x, double y, double radius) {
    if (collider < 0 || collider >= static_cast<int>(colliders.size())) return;
    colliders[collider].holes.push_back({x, y, radius});
}

void World::resolveStrike(int index, const PaneResult& r) {
    if (index < 0 || index >= static_cast<int>(pending_.size())) return;
    const PaneStrike s = pending_[index];
    pending_.erase(pending_.begin() + index);
    Projectile& p = projectiles[s.projectile];
    Collider& pane = colliders[s.collider];
    if (r.holeRadius > 0) {
        const Vec3 l = pane.toLocal(s.point);
        punchHole(s.collider, l.x, l.y, r.holeRadius);
    }
    if (r.shattered) pane.enabled = false;
    if (!r.penetrated) {
        p.state = Projectile::State::Stopped;
        p.fate = "stopped by the " + pane.name;
        events.push_back({s.time, fmt("round stopped by the pane at %.0f m/s", s.speed)});
        return;
    }
    // Through the glass: on from the far face at the speed it kept.
    double tin, tout;
    if (rayBox(pane, s.point - p.dir * 1e-4, p.dir, tin, tout)) p.pos = s.point + p.dir * (tout + 1e-3);
    p.speed = r.speedOut;
    p.time = s.time;
    p.path.push_back(p.pos);
    p.state = Projectile::State::Flying;
    events.push_back({s.time, fmt("round through the pane: %.0f -> %.0f m/s", s.speed, r.speedOut)});
}

void World::spatter(const Vec3& at, const Vec3& normal, const Vec3& axis, const SpatterModel& m, double energy, bool boosted,
                    double time, const std::string& kind, double speed, const Vec3& carrier,
                    double atLeast) {
    const double boost = boosted ? std::max(spatterBoost, 0.1) : 1.0;
    const double volume = m.volumePerJoule * energy * boost;
    if (volume <= 0) return;
    const int src = spray.addSource(kind, at, time);
    const double vmax = std::min(m.maxSpeed, m.speedFactor * speed) * std::pow(boost, 0.15);
    const double sizeUp = std::sqrt(boost);
    // A wound may use a third of the drops the budget has left, and always at
    // least `atLeast` (48 for spatter: a full budget is overrun by that little),
    // so every pellet of a shot load and every later shot still throws its
    // spatter; past that it is sampled, each drop standing for several.
    const double meanVolume = kPi / 6 * std::pow(m.medianDiameter * sizeUp, 3) * std::exp(4.5 * m.sigmaLn * m.sigmaLn);
    const double expected = volume / meanVolume;
    const double room = std::max(0, spray.maxDrops - static_cast<int>(spray.drops.size()));
    const double allowed = std::max(atLeast, room / 3);
    const double weight = expected > allowed ? expected / allowed : 1.0;
    double made = 0;
    int n = 0;
    while (made < volume && n < 4 * allowed) {
        double d = clamp(rng_.lognormal(m.medianDiameter * sizeUp, m.sigmaLn), m.minDiameter, m.maxDiameter * sizeUp);
        const double v = weight * kPi / 6 * d * d * d;
        Vec3 dir = rng_.cone(axis, m.coneHalfAngle);
        if (dot(dir, normal) < 0.05) dir = normalize(dir - normal * (2 * dot(dir, normal) - 0.05));
        const double sp = std::min(vmax, vmax * std::sqrt(rng_.uniform()) * std::cbrt(m.medianDiameter * sizeUp / d));
        Droplet dr;
        dr.d = d;
        dr.t = time;
        dr.source = src;
        dr.weight = weight;
        dr.v = dir * sp + carrier;  // they leave with the body they come from
        const Vec3 u = perpendicular(normal), w = cross(normal, u);
        const double jr = 0.004 * std::sqrt(rng_.uniform()), jp = 2 * kPi * rng_.uniform();
        dr.x = at + normal * (0.5 * d + 4e-3) + u * (jr * std::cos(jp)) + w * (jr * std::sin(jp));
        // Never inside the body it leaves (its capsules overlap round the wound).
        for (int pass = 0; pass < 3; ++pass)
            for (const auto& c : colliders) {
                if (c.role != Role::Dummy || !c.enabled) continue;
                Vec3 cn;
                const double cd = c.sdf(dr.x, &cn);
                if (cd < 0.5 * d + 2e-3) dr.x += cn * (0.5 * d + 2e-3 - cd);
            }
        spray.emit(dr, true);
        made += v;
        ++n;
    }
    initialVolume_ += made;
    if (weight > 1.01)
        log(fmt("  %d drops of %s, %.2f ml (median %.2f mm, up to %.0f m/s; sampled, each for %.1f)", n, kind.c_str(),
                made * 1e6, m.medianDiameter * sizeUp * 1e3, vmax, weight));
    else
        log(fmt("  %d drops of %s, %.2f ml (median %.2f mm, up to %.0f m/s)", n, kind.c_str(), made * 1e6,
                m.medianDiameter * sizeUp * 1e3, vmax));
}

void World::stepProjectiles(double tEnd) {
    for (auto& p : projectiles) {
        if (p.state == Projectile::State::Waiting) {
            if (p.launch >= tEnd) continue;
            p.state = Projectile::State::Flying;
            p.time = p.launch;
            if (p.def->pellets <= 1 || p.id == 0 || projectiles[p.id - 1].def != p.def || projectiles[p.id - 1].launch > p.launch + 0.01)
                events.push_back({p.launch, p.def->pellets > 1 ? fmt("shot fired: %d pellets of %.1f g at %.0f m/s", p.def->pellets, p.def->mass * 1e3, p.speed)
                                                                : fmt("round fired: %.1f g at %.0f m/s", p.def->mass * 1e3, p.speed)});
        }
        int guard = 0;
        while (p.state == Projectile::State::Flying && p.time < tEnd && ++guard < 64) {
            const double rem = tEnd - p.time;
            if (p.inside) {
                // Poncelet: m dv/dt = -(a + b v^2), so v = sqrt(a/b) tan(theta0 - w t)
                // with theta = atan(v sqrt(b/a)), w = sqrt(a b)/m, and the path
                // x = (m/b) ln(cos(theta0 - w t) / cos theta0).
                const double m = p.def->mass, A = p.a, B = p.b;
                const double q = std::sqrt(B / A), om = std::sqrt(A * B) / m;
                const double th0 = std::atan(p.speed * q);
                const double vEnd2 = ((A + B * p.speed * p.speed) * std::exp(-2 * B * p.exitDistance / m) - A) / B;
                const double vEnd = std::sqrt(std::max(vEnd2, 0.0));
                const double tExit = (th0 - std::atan(vEnd * q)) / om;
                if (tExit > rem) {
                    const double th1 = th0 - om * rem;
                    const double s = (m / B) * std::log(std::cos(th1) / std::cos(th0));
                    p.pos += p.dir * s;
                    p.exitDistance -= s;
                    p.speed = std::tan(th1) / q;
                    p.time = tEnd;
                    break;
                }
                p.pos += p.dir * p.exitDistance;
                p.speed = vEnd;
                p.time += tExit;
                p.inside = false;
                p.path.push_back(p.pos);
                const double deposited = 0.5 * p.def->mass * (p.entrySpeed * p.entrySpeed - p.speed * p.speed);
                if (p.speed < 25.5) {
                    p.state = Projectile::State::Stopped;
                    p.fate = "lodged in the dummy";
                    events.push_back({p.time, fmt("round lodged after depositing %.0f J", deposited)});
                    break;
                }
                // Exit wound: forward spatter.
                Vec3 n;
                double best = 1e30;
                int exitCollider = -1;
                // The capsule whose surface it is (the least signed distance: the
                // capsules overlap at the joints, and the round leaves through the
                // outside of their union, not an inner boundary).
                for (const auto& c : colliders)
                    if (c.role == Role::Dummy && c.enabled) {
                        Vec3 nn;
                        const double d = c.sdf(p.pos, &nn);
                        if (d < best) {
                            best = d;
                            n = nn;
                            exitCollider = c.id;
                        }
                    }
                events.push_back({p.time, fmt("exit wound: %.0f m/s left, %.0f J deposited", p.speed, deposited)});
                Attachment at;
                if (ragdoll.active()) at = ragdoll.attach(ragdoll.segmentOf(exitCollider), p.pos, n);
                const double tsave = time_;
                time_ = p.time;
                spatter(p.pos, n, normalize(p.dir + n * 0.25), dummy->forward, deposited, true, p.time, "forward spatter",
                        p.speed, ragdoll.velocityOf(at));
                time_ = tsave;
                Wound wd{p.pos, n, true, p.time, 0, 0, static_cast<int>(spray.sources.size()) - 1, 0, at};
                wounds.push_back(wd);
                continue;
            }
            // Next thing on the path.
            double hit = 1e30;
            int what = -1;
            for (const auto& c : colliders) {
                if (!c.enabled || c.body >= 0) continue;
                double tin = 1e30, tout = 0;
                bool ok = false;
                if (c.kind == ShapeKind::Capsule) ok = rayCapsule(p.pos, p.dir, c.a, c.b, c.radius, tin, tout);
                else if (c.kind == ShapeKind::Box) ok = rayBox(c, p.pos, p.dir, tin, tout);
                else if (c.kind == ShapeKind::Plane) {
                    Vec3 n;
                    const double d = c.sdf(p.pos, &n), dn = dot(p.dir, n);
                    if (dn < -1e-9 && d > 0) {
                        tin = d / -dn;
                        ok = true;
                    }
                }
                if (!ok || tin < 1e-7 || tin >= hit) continue;
                if (c.role == Role::Pane && c.inHole(c.toLocal(p.pos + p.dir * tin))) continue;
                hit = tin;
                what = c.id;
            }
            const double reach = p.speed * rem;
            if (what < 0 || hit > reach) {
                p.pos += p.dir * reach;
                p.time = tEnd;
                if (length(p.pos) > 200) {
                    p.state = Projectile::State::Stopped;
                    p.fate = "out of the scene";
                    p.path.push_back(p.pos);
                }
                break;
            }
            p.pos += p.dir * hit;
            p.time += hit / p.speed;
            p.path.push_back(p.pos);
            Collider& c = colliders[what];
            if (c.role == Role::Dummy && dummy) {
                // The chord through every dummy segment the path is in.
                double end = 0;
                std::vector<std::pair<double, double>> iv;
                for (const auto& o : colliders)
                    if (o.role == Role::Dummy && o.enabled) {
                        double a0, a1;
                        if (rayCapsule(p.pos, p.dir, o.a, o.b, o.radius, a0, a1) && a1 > 0) iv.push_back({a0, a1});
                    }
                std::sort(iv.begin(), iv.end());
                for (const auto& [a0, a1] : iv)
                    if (a0 <= end + 1e-6) end = std::max(end, a1);
                const double area = kPi / 4 * p.def->diameter * p.def->diameter;
                p.a = dummy->tissueStrength * area;
                p.b = 0.5 * dummy->tissueDensity * p.def->cdTissue * area;
                p.entrySpeed = p.speed;
                const double m = p.def->mass;
                double chord = end;
                const double vOut2 = ((p.a + p.b * p.speed * p.speed) * std::exp(-2 * p.b * chord / m) - p.a) / p.b;
                const double stop = 25.0;
                if (vOut2 < stop * stop)  // it stops inside: the path to 25 m/s
                    chord = m / (2 * p.b) * std::log((p.a + p.b * p.speed * p.speed) / (p.a + p.b * stop * stop));
                p.exitDistance = chord;
                p.inside = true;
                const double vOut = vOut2 < stop * stop ? 0.0 : std::sqrt(vOut2);
                const double deposited = 0.5 * m * (p.speed * p.speed - vOut * vOut);
                tissue_ += deposited;
                Vec3 n;
                c.sdf(p.pos, &n);
                events.push_back({p.time, fmt("entry wound in the %s at %.0f m/s; %.0f mm of simulant to cross",
                                              c.name.c_str(), p.speed, end * 1e3)});
                Wound wd{p.pos, n, false, p.time, 0, 0, -1, 0, {}};
                if (ragdoll.active()) {
                    // The momentum the round leaves in the dummy, and the stand lets go.
                    const int seg = ragdoll.segmentOf(c.id);
                    double J = m * (p.speed - vOut);
                    if (knockback > 0) {
                        // A fraction of the energy deposited as the body's motion,
                        // shared out over a shot load's pellets so the whole load
                        // carries that fraction of its whole energy.
                        const double loads = std::max(1, p.def->pellets);
                        J = std::max(J, std::sqrt(2 * ragdoll.mass() * knockback * deposited / loads));
                    }
                    ragdoll.impulse(seg, p.pos, p.dir * J);
                    wd.at = ragdoll.attach(seg, p.pos, n);
                    if (releaseOnHit && !ragdoll.released()) {
                        ragdoll.release(p.time);
                        events.push_back({p.time, "the stand lets go: the dummy is a ragdoll"});
                    }
                }
                // After the push, so the spatter leaves with the body it comes from.
                const double tsave = time_;
                time_ = p.time;
                spatter(p.pos, n, normalize(-p.dir + n * 0.35), dummy->back, deposited, true, p.time, "back spatter", p.speed,
                        ragdoll.velocityOf(wd.at));
                time_ = tsave;
                wd.source = static_cast<int>(spray.sources.size()) - 1;
                wounds.push_back(wd);
                continue;
            }
            if (c.role == Role::Pane) {
                PaneStrike s;
                s.projectile = p.id;
                s.collider = c.id;
                s.point = p.pos;
                const Vec3 l = c.toLocal(p.pos), ld = c.rot.unrotate(p.dir);
                s.x = l.x + c.half.x;
                s.y = l.y + c.half.y;
                s.speed = p.speed;
                s.fromBack = ld.z > 0;
                const double nz = std::abs(ld.z);
                s.yaw = std::atan2(ld.x, nz);
                s.pitch = std::atan2(ld.y, nz);
                s.time = p.time;
                s.glazing = c.glazing;
                events.push_back({p.time, fmt("round reaches the pane at %.0f m/s", p.speed)});
                if (paneResolver) {
                    pending_.push_back(s);
                    resolveStrike(static_cast<int>(pending_.size()) - 1, paneResolver(s));
                } else if (deferPanes) {
                    p.state = Projectile::State::AtPane;
                    pending_.push_back(s);
                } else {
                    pending_.push_back(s);
                    resolveStrike(static_cast<int>(pending_.size()) - 1, estimatePane(s, c, *p.def));
                }
                continue;
            }
            p.state = Projectile::State::Stopped;
            p.fate = "struck the " + c.name;
            events.push_back({p.time, "round struck the " + c.name});
        }
    }
}

void World::carryStains() {
    if (!ragdoll.active()) return;
    for (; stainsSeen_ < spray.stains.size(); ++stainsSeen_) {
        const Stain& st = spray.stains[stainsSeen_];
        if (st.collider < 0 || colliders[st.collider].role != Role::Dummy) continue;
        carried_.push_back({stainsSeen_, ragdoll.attach(ragdoll.segmentOf(st.collider), st.pos, st.normal)});
    }
    if (!ragdoll.moving()) return;
    for (const auto& [k, a] : carried_) {
        Stain& st = spray.stains[k];
        const Vec3 old = st.normal;
        ragdoll.place(a, st.pos, st.normal);
        st.carried = true;
        markMoved(k);
        // The tail turns with the surface.
        const Vec3 ax = cross(old, st.normal);
        const double sn = length(ax);
        if (sn > 1e-9) st.dir = normalize(Quat::axisAngle(ax, std::asin(std::min(1.0, sn))).rotate(st.dir));
    }
}

// ---- Panes in pieces ---------------------------------------------------------------

PanePieces& World::piecesOf(int pane) {
    for (auto& b : broken_)
        if (b.pane == pane) return b;
    broken_.emplace_back();
    broken_.back().pane = pane;
    return broken_.back();
}

const PanePieces* World::panePieces(int pane) const {
    for (const auto& b : broken_)
        if (b.pane == pane) return &b;
    return nullptr;
}

void World::setPaneMap(int pane, int nx, int ny, double cell, const std::vector<int>& front, const std::vector<int>& back,
                       bool laminated) {
    if (pane < 0 || pane >= static_cast<int>(colliders.size()) || colliders[pane].kind != ShapeKind::Box) return;
    PanePieces& P = piecesOf(pane);
    Collider& c = colliders[pane];
    P.map.nx = nx;
    P.map.ny = ny;
    P.map.cell = cell;
    P.map.x0 = -c.half.x;
    P.map.y0 = -c.half.y;
    P.map.front = front;
    P.map.back = back;
    P.map.laminated = laminated;
    int top = 0;
    for (int k : front) top = std::max(top, k + 1);
    for (int k : back) top = std::max(top, k + 1);
    if (static_cast<int>(P.map.held.size()) < top) P.map.held.resize(top, 1);
    c.map = &P.map;
    piecesDirty_ = true;
}

void World::addPanePiece(int pane, int piece, const std::vector<double>& outline, double thickness, bool chip,
                         const std::vector<std::vector<double>>& holes) {
    if (pane < 0 || pane >= static_cast<int>(colliders.size()) || piece < 0) return;
    PanePieces& P = piecesOf(pane);
    if (static_cast<int>(P.pieces.size()) <= piece) P.pieces.resize(piece + 1);
    if (static_cast<int>(P.map.held.size()) <= piece) P.map.held.resize(piece + 1, 1);
    PanePiece& pc = P.pieces[piece];
    pc.outline = outline;
    pc.holes = holes;
    pc.thickness = thickness;
    pc.chip = chip;
    pc.radius = 0;
    for (std::size_t i = 0; i + 1 < outline.size(); i += 2) pc.radius = std::max(pc.radius, std::hypot(outline[i], outline[i + 1]));
    piecesDirty_ = true;
}

void World::setPanePiecePose(int pane, int piece, const Vec3& pos, const Quat& rot, const Vec3& vel, int state) {
    if (pane < 0 || pane >= static_cast<int>(colliders.size()) || piece < 0) return;
    PanePieces& P = piecesOf(pane);
    if (static_cast<int>(P.pieces.size()) <= piece) P.pieces.resize(piece + 1);
    if (static_cast<int>(P.map.held.size()) <= piece) P.map.held.resize(piece + 1, 1);
    const Collider pc0 = colliders[pane];
    // Glass Dynamics' pane frame to the pane's own (centred, +z the front face): z is reflected.
    const Vec3 l{pos.x - pc0.half.x, pos.y - pc0.half.y, pc0.half.z - pos.z};
    const Vec3 wp = pc0.pos + pc0.rot.rotate(l);
    const Quat wq = (pc0.rot * Quat{rot.w, -rot.x, -rot.y, rot.z}).normalized();
    const Vec3 wv = pc0.rot.rotate({vel.x, vel.y, -vel.z});
    PanePiece& pc = P.pieces[piece];
    const bool first = !pc.posed;
    if (first || length(wp - pc.pos) > 1e-7 || std::abs(dot({wq.x, wq.y, wq.z}, {pc.rot.x, pc.rot.y, pc.rot.z}) + wq.w * pc.rot.w) < 1 - 1e-10)
        pc.moved = true;
    pc.pos = wp;
    pc.rot = wq;
    pc.vel = wv;
    pc.posed = true;
    const int was = pc.state;
    pc.state = state;
    P.map.held[piece] = state == 0 || state == 1;
    if (first) piecesDirty_ = true;
    // A loose piece is a slab drops can land on (not the crumbs: chips and
    // pieces under 3 mm across take nothing worth keeping).
    if ((state == 2 || state == 3) && pc.collider < 0 && !pc.chip && pc.radius >= 3e-3 && pc.outline.size() >= 6 &&
        P.loose < maxLoosePieces) {
        Collider c;
        c.name = pc0.name + " piece";
        c.role = Role::Shard;
        c.kind = ShapeKind::Slab;
        c.outline = pc.outline;
        c.outlineHoles = pc.holes;
        c.half = {0, 0, 0.5 * pc.thickness};
        c.surface = pc0.surface;
        c.glazing = pc0.glazing;
        c.pane = pane;
        c.piece = piece;
        c.movable = true;
        c.pos = wp;
        c.rot = wq;
        c.vel = wv;
        pc.collider = addCollider(c);
        ++P.loose;
    } else if (pc.collider >= 0) {
        Collider& c = colliders[pc.collider];
        c.pos = wp;
        c.rot = wq;
        c.vel = state == 2 ? wv : Vec3{};
        c.enabled = state != -1;
    }
    if (state == -1 && was != -1) piecesDirty_ = true;
}

void World::markMoved(std::size_t k) {
    if (movedFlag_.size() < spray.stains.size()) movedFlag_.resize(spray.stains.size(), 0);
    if (!movedFlag_[k]) {
        movedFlag_[k] = 1;
        moved_.push_back(k);
    }
}

std::vector<std::size_t> World::takeMovedStains() {
    std::vector<std::size_t> out;
    out.swap(moved_);
    for (std::size_t k : out) movedFlag_[k] = 0;
    return out;
}

namespace {
bool insidePolygon(const std::vector<double>& poly, double px, double py) {
    bool in = false;
    const std::size_t n = poly.size() / 2;
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const double ax = poly[2 * j], ay = poly[2 * j + 1], bx = poly[2 * i], by = poly[2 * i + 1];
        if ((ay > py) != (by > py) && px < ax + (py - ay) * (bx - ax) / (by - ay)) in = !in;
    }
    return in;
}
}  // namespace

void World::carryPieceStains() {
    if (broken_.empty()) return;
    auto attach = [&](PieceStain& ps, const PanePiece& pc) {
        const Stain& st = spray.stains[ps.stain];
        ps.p = pc.rot.unrotate(st.pos - pc.pos);
        ps.n = pc.rot.unrotate(st.normal);
        ps.d = pc.rot.unrotate(st.dir);
    };
    // New stains on a broken pane or on its pieces.
    for (std::size_t k = pieceSeen_; k < spray.stains.size(); pieceSeen_ = ++k) {
        const Stain& st = spray.stains[k];
        if (st.collider < 0) continue;
        const Collider& c = colliders[st.collider];
        if (c.kind == ShapeKind::Slab) pieceStains_.push_back({k, c.pane, c.piece, {}, {}, {}});
        else if (c.role == Role::Pane) pieceStains_.push_back({k, c.id, -1, {}, {}, {}});
        else continue;
        PieceStain& ps = pieceStains_.back();
        if (ps.piece >= 0) attach(ps, piecesOf(ps.pane).pieces[ps.piece]);
        else piecesDirty_ = true;
    }
    for (auto& ps : pieceStains_) {
        Stain& st = spray.stains[ps.stain];
        if (st.gone) continue;
        PanePieces* P = nullptr;
        for (auto& b : broken_)
            if (b.pane == ps.pane) P = &b;
        if (!P) continue;
        // A stain on the pane: which piece it landed on, once there are pieces.
        if (ps.piece < 0) {
            if (!piecesDirty_ || P->map.nx <= 0) continue;
            const Collider& pane = colliders[ps.pane];
            const Vec3 l = pane.toLocal(st.pos);
            const int k = P->map.pieceAt(l.x, l.y, dot(pane.rot.unrotate(st.normal), {0, 0, 1}) > 0);
            if (k < 0 || k >= static_cast<int>(P->pieces.size()) || !P->pieces[k].posed) continue;
            ps.piece = k;
            attach(ps, P->pieces[k]);
            continue;
        }
        PanePiece* pc = &P->pieces[ps.piece];
        if (pc->state == -1) {
            // The piece broke again: the stain goes with the new piece under it,
            // or with the dice it crumbled into.
            int found = -1;
            for (int j = ps.piece + 1; j < static_cast<int>(P->pieces.size()) && found < 0; ++j) {
                const PanePiece& q = P->pieces[j];
                if (!q.posed || q.state == -1 || q.outline.size() < 6) continue;
                const Vec3 l = q.rot.unrotate(st.pos - q.pos);
                if (std::abs(l.z) <= 0.5 * q.thickness + 2e-3 && insidePolygon(q.outline, l.x, l.y)) found = j;
            }
            if (found < 0) {
                st.gone = true;
                markMoved(ps.stain);
                continue;
            }
            ps.piece = found;
            pc = &P->pieces[found];
            attach(ps, *pc);
        }
        if (!pc->moved) continue;
        st.pos = pc->pos + pc->rot.rotate(ps.p);
        st.normal = pc->rot.rotate(ps.n);
        st.dir = pc->rot.rotate(ps.d);
        if (pc->state != 0 && pc->state != 1) st.carried = true;
        markMoved(ps.stain);
    }
    for (auto& b : broken_)
        for (auto& pc : b.pieces) pc.moved = false;
    piecesDirty_ = false;
}

void World::bleed(double tEnd) {
    if (!dummy || dummy->bleedRate <= 0 || !spray.liquid()) return;
    const LiquidDef& L = *spray.liquid();
    // Tate's law with the Harkins-Brown factor (~0.6) for a 1.5 mm wound edge.
    const double vDrip = 0.6 * 2 * kPi * 0.0015 * L.surfaceTension / (L.density * kGravity);
    const double dDrip = std::cbrt(6 * vDrip / kPi);
    const double boost = std::max(spatterBoost, 0.1);
    for (auto& w : wounds) {
        const double age = time_ - w.time;
        if (age < 0 || age > kBleedFor) continue;
        // The cinema's cut artery: a few pulses a heartbeat apart, each a jet
        // that arcs out across the room (not as in life: a dummy has no heart).
        if (boost > 5 && w.spurts < 3 && age >= 0.12 + 0.8 * w.spurts) {
            const double k = (boost - 5) / 6;
            SpatterModel jet;
            jet.volumePerJoule = 10e-6 * std::pow(k, 1.5) * std::pow(0.8, w.spurts);
            jet.medianDiameter = 1.6e-3;
            jet.sigmaLn = 0.5;
            jet.minDiameter = 0.4e-3;
            jet.maxDiameter = 4.5e-3;
            jet.speedFactor = jet.maxSpeed = (2.5 + 4.5 * k) * std::pow(0.85, w.spurts);
            jet.coneHalfAngle = 0.12;
            const double tsave = time_;
            // A jet is a few hundred big drops at least, or it reads as a sprinkle
            // (they land within a second, and the budget recovers).
            spatter(w.pos + w.normal * 0.01, w.normal, normalize(w.normal + Vec3{0, 0.35, 0}), jet, 1.0, false, time_,
                    "arterial spurt", 1.0, ragdoll.velocityOf(w.at), 400);
            time_ = tsave;
            ++w.spurts;
        }
        const double dt = tEnd - time_;
        // Wounds ooze for a second or two, not a minute: what they drip stays
        // under them and leaves the spatter as it landed.
        w.owed += dummy->bleedRate * std::sqrt(boost) * std::exp(-age / kBleedDecay) * dt;
        while (w.owed >= vDrip) {
            w.owed -= vDrip;
            if (w.bled == 0) w.source = spray.addSource("drip", w.pos, tEnd);
            w.bled += vDrip;
            bledVolume_ += vDrip;
            Droplet d;
            d.d = dDrip;
            d.t = tEnd;
            d.source = w.source;
            // It runs out of the wound and drips off the skin just below it.
            d.x = w.pos + w.normal * (0.5 * dDrip + 0.015) + Vec3{0, -0.02, 0};
            d.v = w.normal * 0.05 + Vec3{0, -0.05, 0} + ragdoll.velocityOf(w.at);
            spray.emit(d);
            initialVolume_ += vDrip;
        }
    }
}

}  // namespace lqd
