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

std::string fmt(const char* f, ...) {
    char buf[512];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof buf, f, ap);
    va_end(ap);
    return buf;
}

// Ray o + t d (|d| = 1) against a sphere: entry and exit.
bool raySphere(const Vec3& o, const Vec3& d, const Vec3& c, double r, double& t0, double& t1) {
    const Vec3 oc = o - c;
    const double b = dot(oc, d), cc = length2(oc) - r * r, disc = b * b - cc;
    if (disc < 0) return false;
    const double s = std::sqrt(disc);
    t0 = -b - s;
    t1 = -b + s;
    return true;
}

// Ray against a capsule (convex: the hull of the cylinder and end spheres).
bool rayCapsule(const Vec3& o, const Vec3& d, const Vec3& a, const Vec3& b, double r, double& tin, double& tout) {
    double lo = 1e30, hi = -1e30, s0, s1;
    if (raySphere(o, d, a, r, s0, s1)) {
        lo = std::min(lo, s0);
        hi = std::max(hi, s1);
    }
    if (raySphere(o, d, b, r, s0, s1)) {
        lo = std::min(lo, s0);
        hi = std::max(hi, s1);
    }
    const Vec3 ax = b - a;
    const double L = length(ax);
    if (L > 1e-9) {
        const Vec3 u = ax / L;
        const Vec3 oa = o - a;
        const Vec3 dp = d - u * dot(d, u), op = oa - u * dot(oa, u);
        const double A = length2(dp), B = 2 * dot(dp, op), C = length2(op) - r * r;
        double c0 = -1e30, c1 = 1e30;
        bool ok = true;
        if (A > 1e-12) {
            const double disc = B * B - 4 * A * C;
            if (disc < 0) ok = false;
            else {
                const double s = std::sqrt(disc);
                c0 = (-B - s) / (2 * A);
                c1 = (-B + s) / (2 * A);
            }
        } else if (C > 0) ok = false;
        if (ok) {
            // Slab between the end planes.
            const double du = dot(d, u), ou = dot(oa, u);
            double p0 = -1e30, p1 = 1e30;
            if (std::abs(du) > 1e-12) {
                p0 = (0 - ou) / du;
                p1 = (L - ou) / du;
                if (p0 > p1) std::swap(p0, p1);
            } else if (ou < 0 || ou > L) ok = false;
            const double e0 = std::max(c0, p0), e1 = std::min(c1, p1);
            if (ok && e0 <= e1) {
                lo = std::min(lo, e0);
                hi = std::max(hi, e1);
            }
        }
    }
    if (lo > hi) return false;
    tin = lo;
    tout = hi;
    return true;
}

// Ray against an oriented box (slabs in the box frame).
bool rayBox(const Collider& c, const Vec3& o, const Vec3& d, double& tin, double& tout) {
    const Vec3 lo = c.toLocal(o), ld = c.rot.unrotate(d);
    double t0 = -1e30, t1 = 1e30;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(ld[i]) < 1e-12) {
            if (std::abs(lo[i]) > c.half[i]) return false;
            continue;
        }
        double a = (-c.half[i] - lo[i]) / ld[i], b = (c.half[i] - lo[i]) / ld[i];
        if (a > b) std::swap(a, b);
        t0 = std::max(t0, a);
        t1 = std::min(t1, b);
    }
    if (t0 > t1) return false;
    tin = t0;
    tout = t1;
    return true;
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
    if (ragdoll.active()) {
        ragdoll.step(dt, colliders);
        ragdoll.pose(colliders);
        for (auto& w : wounds) ragdoll.place(w.at, w.pos, w.normal);
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
        if (dummy && dummy->bleedRate > 0 && time_ - w.time < 20) return false;
    if (ragdoll.moving()) return false;
    return true;
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
        log(fmt("round stopped by the pane at %.0f m/s", s.speed));
        return;
    }
    // Through the glass: on from the far face at the speed it kept.
    double tin, tout;
    if (rayBox(pane, s.point - p.dir * 1e-4, p.dir, tin, tout)) p.pos = s.point + p.dir * (tout + 1e-3);
    p.speed = r.speedOut;
    p.time = s.time;
    p.path.push_back(p.pos);
    p.state = Projectile::State::Flying;
    log(fmt("round through the pane: %.0f -> %.0f m/s", s.speed, r.speedOut));
}

void World::spatter(const Vec3& at, const Vec3& normal, const Vec3& axis, const SpatterModel& m, double energy,
                    double time, const std::string& kind, double speed) {
    const double volume = m.volumePerJoule * energy;
    if (volume <= 0) return;
    const int src = spray.addSource(kind, at, time);
    const double vmax = std::min(m.maxSpeed, m.speedFactor * speed);
    double made = 0;
    int n = 0;
    while (made < volume && n < spray.maxDrops) {
        double d = clamp(rng_.lognormal(m.medianDiameter, m.sigmaLn), m.minDiameter, m.maxDiameter);
        const double v = kPi / 6 * d * d * d;
        Vec3 dir = rng_.cone(axis, m.coneHalfAngle);
        if (dot(dir, normal) < 0.05) dir = normalize(dir - normal * (2 * dot(dir, normal) - 0.05));
        const double sp = std::min(vmax, vmax * std::sqrt(rng_.uniform()) * std::cbrt(m.medianDiameter / d));
        Droplet dr;
        dr.d = d;
        dr.t = time;
        dr.source = src;
        dr.v = dir * sp;
        const Vec3 u = perpendicular(normal), w = cross(normal, u);
        const double jr = 0.004 * std::sqrt(rng_.uniform()), jp = 2 * kPi * rng_.uniform();
        dr.x = at + normal * (0.5 * d + 1e-3) + u * (jr * std::cos(jp)) + w * (jr * std::sin(jp));
        const std::size_t before = spray.drops.size();
        spray.emit(dr);
        if (spray.drops.size() == before) break;
        made += v;
        ++n;
    }
    initialVolume_ += made;
    log(fmt("  %d drops of %s, %.2f ml (median %.2f mm, up to %.0f m/s)", n, kind.c_str(), made * 1e6,
            m.medianDiameter * 1e3, vmax));
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
                for (const auto& c : colliders)
                    if (c.role == Role::Dummy && c.enabled) {
                        Vec3 nn;
                        const double d = std::abs(c.sdf(p.pos, &nn));
                        if (d < best) {
                            best = d;
                            n = nn;
                            exitCollider = c.id;
                        }
                    }
                events.push_back({p.time, fmt("exit wound: %.0f m/s left, %.0f J deposited", p.speed, deposited)});
                const double tsave = time_;
                time_ = p.time;
                spatter(p.pos, n, normalize(p.dir + n * 0.25), dummy->forward, deposited, p.time, "forward spatter",
                        p.speed);
                time_ = tsave;
                Wound wd{p.pos, n, true, p.time, 0, 0, static_cast<int>(spray.sources.size()) - 1, {}};
                if (ragdoll.active()) wd.at = ragdoll.attach(ragdoll.segmentOf(exitCollider), p.pos, n);
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
                Vec3 n;
                c.sdf(p.pos, &n);
                events.push_back({p.time, fmt("entry wound in the %s at %.0f m/s; %.0f mm of simulant to cross",
                                              c.name.c_str(), p.speed, end * 1e3)});
                const double tsave = time_;
                time_ = p.time;
                spatter(p.pos, n, normalize(-p.dir + n * 0.35), dummy->back, deposited, p.time, "back spatter", p.speed);
                time_ = tsave;
                Wound wd{p.pos, n, false, p.time, 0, 0, static_cast<int>(spray.sources.size()) - 1, {}};
                if (ragdoll.active()) {
                    // The momentum the round leaves in the dummy, and the stand lets go.
                    const int seg = ragdoll.segmentOf(c.id);
                    ragdoll.impulse(seg, p.pos, p.dir * (m * (p.speed - vOut)));
                    wd.at = ragdoll.attach(seg, p.pos, n);
                    if (releaseOnHit && !ragdoll.released()) {
                        ragdoll.release(p.time);
                        events.push_back({p.time, "the stand lets go: the dummy is a ragdoll"});
                    }
                }
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
        // The tail turns with the surface.
        const Vec3 ax = cross(old, st.normal);
        const double sn = length(ax);
        if (sn > 1e-9) st.dir = normalize(Quat::axisAngle(ax, std::asin(std::min(1.0, sn))).rotate(st.dir));
    }
}

void World::bleed(double tEnd) {
    if (!dummy || dummy->bleedRate <= 0 || !spray.liquid()) return;
    const LiquidDef& L = *spray.liquid();
    // Tate's law with the Harkins-Brown factor (~0.6) for a 1.5 mm wound edge.
    const double vDrip = 0.6 * 2 * kPi * 0.0015 * L.surfaceTension / (L.density * kGravity);
    const double dDrip = std::cbrt(6 * vDrip / kPi);
    for (auto& w : wounds) {
        const double age = time_ - w.time;
        if (age < 0 || age > 20) continue;
        const double dt = tEnd - time_;
        w.owed += dummy->bleedRate * std::exp(-age / 8.0) * dt;
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
