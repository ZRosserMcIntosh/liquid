// Liquid Dynamics — C ABI for the WebAssembly build (web/build_wasi.sh).
// SPDX-License-Identifier: MIT
//
// One world at a time. Text goes in and out through one shared buffer
// (lqd_text), numbers through float/double arrays the loader views in place.
// Every array starts with its element count where the count is not implied.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "lqd/bpa.hpp"
#include "lqd/export.hpp"
#include "lqd/library.hpp"
#include "lqd/mesher.hpp"
#include "lqd/scenarios.hpp"
#include "lqd/world.hpp"

using namespace lqd;

#define LQD_EXPORT extern "C" __attribute__((visibility("default")))

namespace {

struct Sim {
    Library lib;
    std::vector<std::string> errors;
    std::unique_ptr<World> world;
    SceneConfig config;
    std::string text;
    std::vector<float> f;
    std::vector<double> state;
    std::unique_ptr<Recorder> rec;
    std::vector<std::uint8_t> blob;
    std::vector<std::uint32_t> in;  // the host's numbers going in (lqd_input)
};

Sim& sim() {
    static Sim s;
    return s;
}

void putText(const std::string& s) { sim().text = s; }

void jsonString(std::string& o, const std::string& s) {
    o += '"';
    for (char c : s) {
        if (c == '"' || c == '\\') o += '\\';
        if (c == '\n') {
            o += "\\n";
            continue;
        }
        o += c;
    }
    o += '"';
}
void jsonNum(std::string& o, double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.9g", std::isfinite(v) ? v : 0.0);
    o += b;
}
void jsonVec(std::string& o, const Vec3& v) {
    o += "[";
    jsonNum(o, v.x);
    o += ",";
    jsonNum(o, v.y);
    o += ",";
    jsonNum(o, v.z);
    o += "]";
}
void provenance(std::string& o, const Provenance& p) {
    o += "{";
    bool first = true;
    for (const auto& [k, v] : p) {
        if (!first) o += ",";
        first = false;
        jsonString(o, k);
        o += ":";
        jsonString(o, v);
    }
    o += "}";
}
float* floats() { return sim().f.data(); }

}  // namespace

// ---- Text and library ----------------------------------------------------------------

// A buffer of n bytes for the host to write into (and the text output buffer).
LQD_EXPORT char* lqd_text(int n) {
    Sim& S = sim();
    if (n > 0) S.text.assign(static_cast<std::size_t>(n), '\0');
    return S.text.data();
}
LQD_EXPORT int lqd_text_length() { return static_cast<int>(std::strlen(sim().text.c_str())); }

LQD_EXPORT int lqd_load(int n) {
    Sim& S = sim();
    const std::string doc(S.text.data(), static_cast<std::size_t>(n));
    for (auto& e : S.lib.load(doc)) S.errors.push_back(e);
    return static_cast<int>(S.errors.size());
}
LQD_EXPORT int lqd_resolve() {
    Sim& S = sim();
    for (auto& e : S.lib.resolve()) S.errors.push_back(e);
    return static_cast<int>(S.errors.size());
}
LQD_EXPORT int lqd_error_text(int i) {
    Sim& S = sim();
    putText(i >= 0 && i < static_cast<int>(S.errors.size()) ? S.errors[i] : "");
    return lqd_text_length();
}

LQD_EXPORT int lqd_catalog() {
    const Library& L = sim().lib;
    std::string o = "{\"liquids\":[";
    for (std::size_t i = 0; i < L.liquids.size(); ++i) {
        const auto& d = L.liquids[i];
        o += i ? ",{" : "{";
        o += "\"id\":"; jsonString(o, d.id);
        o += ",\"name\":"; jsonString(o, d.name);
        o += ",\"category\":"; jsonString(o, d.category);
        o += ",\"description\":"; jsonString(o, d.description);
        o += ",\"density\":"; jsonNum(o, d.density);
        o += ",\"surfaceTension\":"; jsonNum(o, d.surfaceTension);
        o += ",\"viscosity\":"; jsonNum(o, d.viscosity);
        o += ",\"impactViscosity\":"; jsonNum(o, d.impactViscosity());
        o += ",\"shearThinning\":"; o += d.rheology == Rheology::Carreau ? "true" : "false";
        o += ",\"color\":"; jsonVec(o, d.color);
        o += ",\"stainColor\":"; jsonVec(o, d.stainColor);
        o += ",\"opacity\":"; jsonNum(o, d.opacity);
        o += ",\"provenance\":"; provenance(o, d.provenance);
        o += "}";
    }
    o += "],\"surfaces\":[";
    for (std::size_t i = 0; i < L.surfaces.size(); ++i) {
        const auto& d = L.surfaces[i];
        o += i ? ",{" : "{";
        o += "\"id\":"; jsonString(o, d.id);
        o += ",\"name\":"; jsonString(o, d.name);
        o += ",\"description\":"; jsonString(o, d.description);
        o += ",\"youngs\":"; jsonNum(o, d.youngs);
        o += ",\"layer\":"; jsonNum(o, d.layer);
        o += ",\"absorbency\":"; jsonNum(o, d.absorbency);
        o += ",\"splashFactor\":"; jsonNum(o, d.splashFactor);
        o += ",\"color\":"; jsonVec(o, d.color);
        o += ",\"provenance\":"; provenance(o, d.provenance);
        o += "}";
    }
    o += "],\"vessels\":[";
    for (std::size_t i = 0; i < L.vessels.size(); ++i) {
        const auto& d = L.vessels[i];
        const VesselGeometry g(d, 32);
        o += i ? ",{" : "{";
        o += "\"id\":"; jsonString(o, d.id);
        o += ",\"name\":"; jsonString(o, d.name);
        o += ",\"category\":"; jsonString(o, d.category);
        o += ",\"description\":"; jsonString(o, d.description);
        o += ",\"mass\":"; jsonNum(o, g.mass());
        o += ",\"capacity\":"; jsonNum(o, g.cavityCapacity());
        o += ",\"height\":"; jsonNum(o, g.height());
        o += ",\"fill\":"; jsonNum(o, d.fill);
        o += ",\"strength\":"; jsonNum(o, d.strength);
        o += ",\"liquid\":"; jsonString(o, d.liquid);
        o += ",\"provenance\":"; provenance(o, d.provenance);
        o += "}";
    }
    o += "],\"projectiles\":[";
    for (std::size_t i = 0; i < L.projectiles.size(); ++i) {
        const auto& d = L.projectiles[i];
        o += i ? ",{" : "{";
        o += "\"id\":"; jsonString(o, d.id);
        o += ",\"name\":"; jsonString(o, d.name);
        o += ",\"description\":"; jsonString(o, d.description);
        o += ",\"mass\":"; jsonNum(o, d.mass);
        o += ",\"diameter\":"; jsonNum(o, d.diameter);
        o += ",\"speed\":"; jsonNum(o, d.speed);
        o += ",\"glassImpactor\":"; jsonString(o, d.glassImpactor);
        o += ",\"pellets\":"; jsonNum(o, d.pellets);
        o += ",\"spread\":"; jsonNum(o, d.spread);
        o += ",\"provenance\":"; provenance(o, d.provenance);
        o += "}";
    }
    o += "],\"dummies\":[";
    for (std::size_t i = 0; i < L.dummies.size(); ++i) {
        const auto& d = L.dummies[i];
        o += i ? ",{" : "{";
        o += "\"id\":"; jsonString(o, d.id);
        o += ",\"name\":"; jsonString(o, d.name);
        o += ",\"description\":"; jsonString(o, d.description);
        o += ",\"liquid\":"; jsonString(o, d.liquid);
        o += ",\"provenance\":"; provenance(o, d.provenance);
        o += "}";
    }
    o += "],\"tiers\":[";
    for (int t = 0; t <= 4; ++t) {
        const Budget b = budgetForTier(t);
        o += t ? ",{" : "{";
        o += "\"name\":"; jsonString(o, b.name);
        o += ",\"spacing\":"; jsonNum(o, b.spacing);
        o += ",\"drops\":"; jsonNum(o, b.maxDrops);
        o += "}";
    }
    o += "]}";
    putText(o);
    return lqd_text_length();
}

// ---- Scene -------------------------------------------------------------------------

// Builds a scene from the JSON config in the text buffer (n bytes). Returns 1,
// or 0 with the error in the text buffer.
LQD_EXPORT int lqd_build(int n) {
    Sim& S = sim();
    const std::string doc(S.text.data(), static_cast<std::size_t>(n));
    SceneConfig c;
    std::string err;
    if (!parseSceneConfig(doc, c, &err)) {
        putText(err);
        return 0;
    }
    auto w = std::make_unique<World>(S.lib, budgetForTier(c.tier), c.seed);
    if (!buildScene(*w, c, &err)) {
        putText(err);
        return 0;
    }
    S.world = std::move(w);
    S.config = c;
    S.rec.reset();
    putText("");
    return 1;
}

// The scene's shapes for the renderer (JSON in the text buffer).
LQD_EXPORT int lqd_scene() {
    Sim& S = sim();
    std::string o = "{\"colliders\":[";
    if (S.world) {
        bool first = true;
        for (const auto& c : S.world->colliders) {
            if (!first) o += ",";
            first = false;
            o += "{\"id\":"; jsonNum(o, c.id);
            o += ",\"name\":"; jsonString(o, c.name);
            o += ",\"role\":"; jsonString(o, roleName(c.role));
            o += ",\"kind\":"; jsonString(o, c.kind == ShapeKind::Plane ? "plane" : c.kind == ShapeKind::Box ? "box" : c.kind == ShapeKind::Capsule ? "capsule" : c.kind == ShapeKind::Slab ? "slab" : "shell");
            o += ",\"pos\":"; jsonVec(o, c.pos);
            o += ",\"rot\":[";
            jsonNum(o, c.rot.x); o += ","; jsonNum(o, c.rot.y); o += ","; jsonNum(o, c.rot.z); o += ","; jsonNum(o, c.rot.w);
            o += "],\"half\":"; jsonVec(o, c.half);
            o += ",\"a\":"; jsonVec(o, c.a);
            o += ",\"b\":"; jsonVec(o, c.b);
            o += ",\"radius\":"; jsonNum(o, c.radius);
            o += ",\"color\":"; jsonVec(o, c.color);
            o += ",\"surface\":"; jsonString(o, c.surface ? c.surface->id : "");
            o += ",\"glazing\":"; jsonString(o, c.glazing);
            o += ",\"body\":"; jsonNum(o, c.body);
            o += ",\"enabled\":"; o += c.enabled ? "true" : "false";
            o += "}";
        }
        const LiquidDef* L = S.world->fluid.liquid() ? S.world->fluid.liquid() : S.world->spray.liquid();
        o += "],\"liquid\":";
        jsonString(o, L ? L->id : "");
        o += ",\"spacing\":"; jsonNum(o, S.world->fluid.settings().spacing);
        o += ",\"scenario\":"; jsonString(o, S.config.scenario);
        o += ",\"tier\":"; jsonString(o, S.world->budget().name);
        const ShotLine sl = shotLine(S.config);
        o += ",\"muzzle\":"; jsonVec(o, sl.muzzle);
        o += ",\"aim\":"; jsonVec(o, sl.aim);
        o += "}";
    } else {
        o += "]}";
    }
    putText(o);
    return lqd_text_length();
}

LQD_EXPORT int lqd_advance(double t, int maxSteps) { return sim().world ? sim().world->advance(t, maxSteps) : 0; }
LQD_EXPORT double lqd_time() { return sim().world ? sim().world->time() : 0; }
LQD_EXPORT double lqd_next_step() { return sim().world ? sim().world->nextStep() : 0; }

// [0] time, [1] quiet, [2] particles, [3] drops, [4] stains, [5] bodies,
// [6] events, [7] steps, [8] spacing, [9] wetted m^2, [10] covered m^2,
// [11] volume poured m^3, [12] bulk, [13] in the air, [14] stained,
// [15] absorbed, [16] lost, [17] breaks, [18] fastest liquid m/s,
// [19] projectiles, [20] pending pane strikes, [21] sources, [22] next step s,
// [23] wounds, [24] breakups, [25] splashes, [26] colliders, [27] bled m^3
LQD_EXPORT double* lqd_state() {
    Sim& S = sim();
    S.state.assign(32, 0);
    if (!S.world) return S.state.data();
    const World& w = *S.world;
    const VolumeBook v = w.volumes();
    double* s = S.state.data();
    s[0] = w.time();
    s[1] = w.quiet() ? 1 : 0;
    s[2] = static_cast<double>(w.fluid.size());
    s[3] = static_cast<double>(w.spray.drops.size());
    s[4] = static_cast<double>(w.spray.stains.size());
    s[5] = static_cast<double>(w.bodies.size());
    s[6] = static_cast<double>(w.events.size());
    s[7] = w.steps();
    s[8] = w.fluid.settings().spacing;
    s[9] = w.wettedArea();
    s[10] = w.coveredArea();
    s[11] = v.initial;
    s[12] = v.bulk;
    s[13] = v.air;
    s[14] = v.stained;
    s[15] = v.absorbed;
    s[16] = v.lost;
    s[17] = w.breaks();
    s[18] = w.fluid.size() ? w.fluid.maxSpeed() : 0;
    s[19] = static_cast<double>(w.projectiles.size());
    s[20] = static_cast<double>(w.pendingStrikes().size());
    s[21] = static_cast<double>(w.spray.sources.size());
    s[22] = w.nextStep();
    s[23] = static_cast<double>(w.wounds.size());
    s[24] = w.spray.breakups;
    s[25] = w.spray.splashes;
    s[26] = static_cast<double>(w.colliders.size());
    s[27] = v.bled;
    return s;
}

// Bulk particles: [n, then x y z speed per particle].
LQD_EXPORT float* lqd_particles() {
    Sim& S = sim();
    S.f.clear();
    if (!S.world) {
        S.f.push_back(0);
        return floats();
    }
    const Fluid& F = S.world->fluid;
    S.f.reserve(1 + 4 * F.size());
    S.f.push_back(static_cast<float>(F.size()));
    for (std::size_t i = 0; i < F.size(); ++i)
        S.f.insert(S.f.end(), {float(F.x[i].x), float(F.x[i].y), float(F.x[i].z), float(length(F.v[i]))});
    return floats();
}

// Droplets in the air: [n, then 8 per drop: x y z, diameter, vx vy vz, blob of the bulk].
LQD_EXPORT float* lqd_drops() {
    Sim& S = sim();
    S.f.clear();
    if (!S.world) {
        S.f.push_back(0);
        return floats();
    }
    const auto& D = S.world->spray.drops;
    S.f.reserve(1 + 8 * D.size());
    S.f.push_back(static_cast<float>(D.size()));
    for (const auto& d : D)
        S.f.insert(S.f.end(), {float(d.x.x), float(d.x.y), float(d.x.z), float(d.d), float(d.v.x), float(d.v.y),
                               float(d.v.z), d.fromBulk ? 1.0f : 0.0f});
    return floats();
}

// The ragdoll: [n, then 7 per capsule: a xyz, b xyz, radius]; [0] = 0 without one.
LQD_EXPORT float* lqd_dummy() {
    Sim& S = sim();
    S.f.assign(1, 0);
    if (!S.world || !S.world->ragdoll.active()) return floats();
    const Ragdoll& R = S.world->ragdoll;
    for (int k = 0; k < R.segments(); ++k) {
        const Vec3 a = R.a(k), b = R.b(k);
        S.f.insert(S.f.end(), {float(a.x), float(a.y), float(a.z), float(b.x), float(b.y), float(b.z), float(R.radius(k))});
        S.f[0] += 1;
    }
    return floats();
}

// Stains on the dummy, where they are now: [n, then 10 per stain: index, pos xyz, normal xyz, tail xyz].
LQD_EXPORT float* lqd_dummy_stains() {
    Sim& S = sim();
    S.f.assign(1, 0);
    if (!S.world) return floats();
    for (const auto& [k, a] : S.world->carried()) {
        const Stain& st = S.world->spray.stains[k];
        S.f.insert(S.f.end(), {float(k), float(st.pos.x), float(st.pos.y), float(st.pos.z), float(st.normal.x),
                               float(st.normal.y), float(st.normal.z), float(st.dir.x), float(st.dir.y), float(st.dir.z)});
        S.f[0] += 1;
    }
    return floats();
}

// Stains from index `from`: [n, then 16 per stain: pos xyz, normal xyz, tail
// xyz, width, length, impact angle, time, collider, satellite, source].
LQD_EXPORT float* lqd_stains(int from) {
    Sim& S = sim();
    S.f.clear();
    S.f.push_back(0);
    if (!S.world) return floats();
    const auto& St = S.world->spray.stains;
    for (std::size_t i = static_cast<std::size_t>(std::max(from, 0)); i < St.size(); ++i) {
        const Stain& s = St[i];
        S.f.insert(S.f.end(), {float(s.pos.x), float(s.pos.y), float(s.pos.z), float(s.normal.x), float(s.normal.y),
                               float(s.normal.z), float(s.dir.x), float(s.dir.y), float(s.dir.z), float(s.width),
                               float(s.length), float(s.alpha), float(s.time), float(s.collider),
                               s.satellite ? 1.0f : 0.0f, float(s.source)});
        S.f[0] += 1;
    }
    return floats();
}

// Bodies: [n, then 12 per body: x y z, qx qy qz qw, shown, shard, mass, born, broken at].
LQD_EXPORT float* lqd_bodies() {
    Sim& S = sim();
    S.f.clear();
    S.f.push_back(0);
    if (!S.world) return floats();
    for (const auto& b : S.world->bodies) {
        S.f.insert(S.f.end(), {float(b.x.x), float(b.x.y), float(b.x.z), float(b.q.x), float(b.q.y), float(b.q.z),
                               float(b.q.w), b.broken ? 0.0f : 1.0f, float(b.shard), float(b.mass), float(b.born),
                               float(b.brokenAt)});
        S.f[0] += 1;
    }
    return floats();
}

namespace {
float* packMesh(const std::vector<float>& pos, const std::vector<float>& nrm, const std::vector<std::uint32_t>& tri) {
    Sim& S = sim();
    S.f.clear();
    S.f.push_back(static_cast<float>(pos.size() / 3));
    S.f.push_back(static_cast<float>(tri.size() / 3));
    S.f.insert(S.f.end(), pos.begin(), pos.end());
    S.f.insert(S.f.end(), nrm.begin(), nrm.end());
    for (auto t : tri) S.f.push_back(static_cast<float>(t));
    return floats();
}
}  // namespace

// A body's mesh in its own frame: [vertices, triangles, positions, normals, indices].
LQD_EXPORT float* lqd_body_mesh(int body) {
    Sim& S = sim();
    if (!S.world) return packMesh({}, {}, {});
    const VesselMesh m = bodyMesh(*S.world, body);
    std::vector<float> pos, nrm;
    for (std::size_t i = 0; i < m.pos.size(); ++i) {
        pos.insert(pos.end(), {float(m.pos[i].x), float(m.pos[i].y), float(m.pos[i].z)});
        nrm.insert(nrm.end(), {float(m.nrm[i].x), float(m.nrm[i].y), float(m.nrm[i].z)});
    }
    return packMesh(pos, nrm, m.tri);
}

// The liquid's surface now: [vertices, triangles, positions, normals, indices].
LQD_EXPORT float* lqd_surface(double cellFactor) {
    Sim& S = sim();
    if (!S.world || S.world->fluid.size() == 0) return packMesh({}, {}, {});
    MesherSettings ms;
    const double h = S.world->fluid.settings().spacing;
    if (cellFactor > 0) ms.cell = cellFactor * h;
    const SurfaceMesh m = meshLiquid(S.world->fluid.x, h, ms);
    return packMesh(m.pos, m.nrm, m.tri);
}

// Event log lines from index `from` ("time<TAB>text", newline separated).
LQD_EXPORT int lqd_events(int from) {
    Sim& S = sim();
    std::string o;
    if (S.world)
        for (std::size_t i = static_cast<std::size_t>(std::max(from, 0)); i < S.world->events.size(); ++i) {
            char b[32];
            std::snprintf(b, sizeof b, "%.6f\t", S.world->events[i].time);
            o += b + S.world->events[i].text + "\n";
        }
    putText(o);
    return lqd_text_length();
}

// Rounds: [n, then 10 per round: x y z, dir xyz, speed, state (0 waiting,
// 1 flying, 2 at a pane, 3 stopped), launch, path corners].
LQD_EXPORT float* lqd_projectiles() {
    Sim& S = sim();
    S.f.clear();
    S.f.push_back(0);
    if (!S.world) return floats();
    for (const auto& p : S.world->projectiles) {
        S.f.insert(S.f.end(), {float(p.pos.x), float(p.pos.y), float(p.pos.z), float(p.dir.x), float(p.dir.y),
                               float(p.dir.z), float(p.speed), float(static_cast<int>(p.state)), float(p.launch),
                               float(p.path.size())});
        S.f[0] += 1;
    }
    return floats();
}
// The corners of round k's path: [n, then x y z].
LQD_EXPORT float* lqd_projectile_path(int k) {
    Sim& S = sim();
    S.f.clear();
    S.f.push_back(0);
    if (!S.world || k < 0 || k >= static_cast<int>(S.world->projectiles.size())) return floats();
    for (const auto& c : S.world->projectiles[k].path) {
        S.f.insert(S.f.end(), {float(c.x), float(c.y), float(c.z)});
        S.f[0] += 1;
    }
    return floats();
}
LQD_EXPORT int lqd_shoot(int projectile, double ox, double oy, double oz, double dx, double dy, double dz, double speed,
                         double delay) {
    Sim& S = sim();
    if (!S.world || projectile < 0 || projectile >= static_cast<int>(S.lib.projectiles.size())) return -1;
    const ProjectileDef& p = S.lib.projectiles[projectile];
    return S.world->shoot(p, {ox, oy, oz}, {dx, dy, dz}, speed > 0 ? speed : p.speed, delay);
}

// Wounds: [n, then 8 per wound: x y z, normal xyz, exit, time].
LQD_EXPORT float* lqd_wounds() {
    Sim& S = sim();
    S.f.clear();
    S.f.push_back(0);
    if (!S.world) return floats();
    for (const auto& w : S.world->wounds) {
        S.f.insert(S.f.end(), {float(w.pos.x), float(w.pos.y), float(w.pos.z), float(w.normal.x), float(w.normal.y),
                               float(w.normal.z), w.exit ? 1.0f : 0.0f, float(w.time)});
        S.f[0] += 1;
    }
    return floats();
}

// Spatter sources as JSON (text buffer).
LQD_EXPORT int lqd_sources() {
    Sim& S = sim();
    std::string o = "[";
    if (S.world)
        for (std::size_t i = 0; i < S.world->spray.sources.size(); ++i) {
            const Source& s = S.world->spray.sources[i];
            if (i) o += ",";
            o += "{\"kind\":"; jsonString(o, s.kind);
            o += ",\"origin\":"; jsonVec(o, s.origin);
            o += ",\"time\":"; jsonNum(o, s.time);
            o += ",\"volume\":"; jsonNum(o, s.volume);
            o += ",\"drops\":"; jsonNum(o, s.drops);
            o += "}";
        }
    o += "]";
    putText(o);
    return lqd_text_length();
}

// ---- Panes (the Glass Dynamics hand-off) ---------------------------------------------

// Pending pane strikes as JSON (text buffer).
LQD_EXPORT int lqd_strikes() {
    Sim& S = sim();
    std::string o = "[";
    if (S.world)
        for (std::size_t i = 0; i < S.world->pendingStrikes().size(); ++i) {
            const PaneStrike& s = S.world->pendingStrikes()[i];
            if (i) o += ",";
            o += "{\"projectile\":"; jsonNum(o, s.projectile);
            o += ",\"collider\":"; jsonNum(o, s.collider);
            o += ",\"point\":"; jsonVec(o, s.point);
            o += ",\"x\":"; jsonNum(o, s.x);
            o += ",\"y\":"; jsonNum(o, s.y);
            o += ",\"speed\":"; jsonNum(o, s.speed);
            o += ",\"yaw\":"; jsonNum(o, s.yaw);
            o += ",\"pitch\":"; jsonNum(o, s.pitch);
            o += ",\"fromBack\":"; o += s.fromBack ? "true" : "false";
            o += ",\"time\":"; jsonNum(o, s.time);
            o += ",\"glazing\":"; jsonString(o, s.glazing);
            o += ",\"impactor\":";
            const Projectile& p = S.world->projectiles[s.projectile];
            jsonString(o, p.def ? p.def->glassImpactor : "");
            o += "}";
        }
    o += "]";
    putText(o);
    return lqd_text_length();
}
LQD_EXPORT int lqd_resolve_strike(int index, int penetrated, double speedOut, double holeRadius, int shattered) {
    Sim& S = sim();
    if (!S.world) return 0;
    PaneResult r;
    r.penetrated = penetrated != 0;
    r.speedOut = speedOut;
    r.holeRadius = holeRadius;
    r.shattered = shattered != 0;
    S.world->resolveStrike(index, r);
    return 1;
}
// Resolves a pending strike with the built-in estimate (no Glass Dynamics on the page).
LQD_EXPORT int lqd_estimate_strike(int index) {
    Sim& S = sim();
    if (!S.world || index < 0 || index >= static_cast<int>(S.world->pendingStrikes().size())) return 0;
    const PaneStrike s = S.world->pendingStrikes()[index];
    const PaneResult r = World::estimatePane(s, S.world->colliders[s.collider], *S.world->projectiles[s.projectile].def);
    S.world->resolveStrike(index, r);
    return 1;
}
// ---- Panes in pieces (Glass Dynamics) ---------------------------------------------------

// A buffer of n 32-bit words for the host to write numbers into (int32 or float32).
LQD_EXPORT std::uint32_t* lqd_input(int n) {
    Sim& S = sim();
    if (n > 0 && static_cast<int>(S.in.size()) < n) S.in.resize(static_cast<std::size_t>(n));
    return S.in.data();
}
// The pane's fragment map: int32 input, nx * ny piece indices of the front ply
// then (if `back`) of the back ply, rows from the bottom.
LQD_EXPORT void lqd_pane_map(int pane, int nx, int ny, double cell, int back, int laminated) {
    Sim& S = sim();
    if (!S.world || nx <= 0 || ny <= 0) return;
    const std::size_t n = static_cast<std::size_t>(nx) * ny;
    const std::int32_t* p = reinterpret_cast<const std::int32_t*>(S.in.data());
    std::vector<int> f(p, p + n), b;
    if (back) b.assign(p + n, p + 2 * n);
    S.world->setPaneMap(pane, nx, ny, cell, f, b, laminated != 0);
}
// A piece: float32 input, its outline as x, y pairs about its centroid (n
// points), then the number of holes and for each its point count and pairs.
LQD_EXPORT void lqd_pane_piece(int pane, int piece, int n, double thickness, int chip) {
    Sim& S = sim();
    if (!S.world || n < 0) return;
    const float* p = reinterpret_cast<const float*>(S.in.data());
    std::vector<std::vector<double>> holes;
    const float* q = p + 2 * n;
    for (int h = 0, nh = static_cast<int>(*q++); h < nh; ++h) {
        const int m = static_cast<int>(*q++);
        holes.emplace_back(q, q + 2 * m);
        q += 2 * m;
    }
    S.world->addPanePiece(pane, piece, std::vector<double>(p, p + 2 * n), thickness, chip != 0, holes);
}
// Every piece's pose: float32 input, 12 per piece as Glass Dynamics gives them
// (pos xyz, quat wxyz, vel xyz, state, born).
LQD_EXPORT void lqd_pane_poses(int pane, int n) {
    Sim& S = sim();
    if (!S.world) return;
    const float* p = reinterpret_cast<const float*>(S.in.data());
    for (int k = 0; k < n; ++k, p += 12)
        S.world->setPanePiecePose(pane, k, {p[0], p[1], p[2]}, Quat{p[3], p[4], p[5], p[6]}, {p[7], p[8], p[9]},
                                  static_cast<int>(std::lround(p[10])));
}
// Stains that moved since the last call: [n, then 11 per stain: index, pos xyz,
// normal xyz, tail xyz, gone (it went with glass that crumbled)].
LQD_EXPORT float* lqd_moved_stains() {
    Sim& S = sim();
    S.f.assign(1, 0);
    if (!S.world) return floats();
    for (std::size_t k : S.world->takeMovedStains()) {
        const Stain& st = S.world->spray.stains[k];
        S.f.insert(S.f.end(), {float(k), float(st.pos.x), float(st.pos.y), float(st.pos.z), float(st.normal.x),
                               float(st.normal.y), float(st.normal.z), float(st.dir.x), float(st.dir.y), float(st.dir.z),
                               st.gone ? 1.0f : 0.0f});
        S.f[0] += 1;
    }
    return floats();
}

LQD_EXPORT void lqd_collider_enable(int id, int on) {
    Sim& S = sim();
    if (S.world && id >= 0 && id < static_cast<int>(S.world->colliders.size())) S.world->colliders[id].enabled = on != 0;
}
LQD_EXPORT void lqd_hole(int id, double x, double y, double r) {
    if (sim().world) sim().world->punchHole(id, x, y, r);
}

// ---- Stain analysis -------------------------------------------------------------------

// Bloodstain pattern analysis of one surface; JSON in the text buffer.
LQD_EXPORT int lqd_bpa(int collider, int source, double readErrorMm, int satellites) {
    Sim& S = sim();
    std::string o = "{";
    if (S.world) {
        BpaOptions opt;
        opt.collider = collider;
        opt.source = source;
        opt.readError = readErrorMm * 1e-3;
        opt.satellites = satellites != 0;
        const BpaResult r = analyseStains(S.world->spray, S.world->colliders, opt);
        o += "\"ok\":"; o += r.ok ? "true" : "false";
        o += ",\"considered\":"; jsonNum(o, r.considered);
        o += ",\"used\":"; jsonNum(o, r.used);
        o += ",\"convergence\":"; jsonVec(o, r.convergence);
        o += ",\"convergenceRms\":"; jsonNum(o, r.convergenceRms);
        o += ",\"height\":"; jsonNum(o, r.height);
        o += ",\"heightLow\":"; jsonNum(o, r.heightLow);
        o += ",\"heightHigh\":"; jsonNum(o, r.heightHigh);
        o += ",\"origin\":"; jsonVec(o, r.origin);
        o += ",\"normal\":"; jsonVec(o, r.planeNormal);
        o += ",\"meanAngleError\":"; jsonNum(o, r.meanAngleError);
        o += ",\"truth\":"; o += r.truth ? "true" : "false";
        o += ",\"trueOrigin\":"; jsonVec(o, r.trueOrigin);
        o += ",\"convergenceError\":"; jsonNum(o, r.convergenceError);
        o += ",\"heightError\":"; jsonNum(o, r.heightError);
        o += ",\"summary\":"; jsonString(o, r.summary);
        o += ",\"lines\":[";
        for (std::size_t i = 0; i < r.lines.size(); ++i) {
            const BpaLine& l = r.lines[i];
            if (i) o += ",";
            o += "["; jsonNum(o, l.stain);
            o += ","; jsonNum(o, l.alpha);
            o += ","; jsonNum(o, l.alphaTrue);
            o += ","; jsonNum(o, l.distance);
            o += ","; jsonNum(o, l.height);
            o += "]";
        }
        o += "]";
    }
    o += "}";
    putText(o);
    return lqd_text_length();
}

// ---- Recording and export --------------------------------------------------------------

LQD_EXPORT void lqd_record(double fps, int maxFrames) {
    Sim& S = sim();
    S.rec = std::make_unique<Recorder>(fps, static_cast<std::size_t>(std::max(maxFrames, 1)));
}
// Captures a frame if one is due (call after each advance); returns frames recorded.
LQD_EXPORT int lqd_capture() {
    Sim& S = sim();
    if (!S.rec || !S.world) return 0;
    S.rec->capture(*S.world);
    return static_cast<int>(S.rec->frames.size());
}
LQD_EXPORT std::uint8_t* lqd_blob() { return sim().blob.data(); }
// Binary glTF of the recording (or of this moment if nothing was recorded); returns its length.
LQD_EXPORT int lqd_export_glb(int scene, int stride) {
    Sim& S = sim();
    if (!S.world) return 0;
    Recorder now(1);
    if (!S.rec || S.rec->frames.empty()) now.captureNow(*S.world);
    ExportOptions o;
    o.scene = scene != 0;
    o.stride = std::max(1, stride);
    S.blob = exportGlb(*S.world, S.rec && !S.rec->frames.empty() ? *S.rec : now, o);
    return static_cast<int>(S.blob.size());
}
// PLY point cache of this moment; returns its length.
LQD_EXPORT int lqd_export_ply() {
    Sim& S = sim();
    if (!S.world) return 0;
    Recorder now(1);
    now.captureNow(*S.world);
    S.blob = exportPly(now.frames[0], S.world->fluid.settings().spacing);
    return static_cast<int>(S.blob.size());
}
LQD_EXPORT int lqd_export_stains() {
    Sim& S = sim();
    putText(S.world ? exportStainsJson(*S.world) : "{}");
    return lqd_text_length();
}

// Milliseconds per step of a 2,000-particle splash on this device (for the
// page's automatic tier).
LQD_EXPORT double lqd_bench() {
    Sim& S = sim();
    World w(S.lib, budgetForTier(2), 1);
    const LiquidDef* water = S.lib.liquid("water");
    if (!water) return 0;
    w.addPlane(Role::Floor, "floor", {0, 0, 0}, {0, 1, 0}, S.lib.surface("ceramic_tile"));
    w.addLiquidBox(*water, {0, 0.02, 0}, {0.06, 0.08, 0.06}, {0, -1, 0});
    w.advance(0.01);
    const auto t0 = std::chrono::steady_clock::now();
    const int n = 12;
    for (int i = 0; i < n; ++i) w.step(w.nextStep());
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() / n;
}
