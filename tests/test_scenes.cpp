// Liquid Dynamics — scene tests: the liquid at rest and in a dam break, the
// glass, the bottle, the shot, drips and impact-angle cards, determinism.
// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>

#include "lqd/bpa.hpp"
#include "lqd/export.hpp"
#include "lqd/json.hpp"
#include "lqd/scenarios.hpp"
#include "lqd_test.hpp"

using namespace lqd;

namespace {
const Library& lib() {
    static Library l = [] {
        Library x;
        loadDefaultLibrary(x);
        return x;
    }();
    return l;
}

struct Run {
    std::unique_ptr<World> w;
    bool ok = false;
};
Run run(const std::string& json, double until, double chunk = 0.02) {
    Run r;
    SceneConfig c;
    std::string err;
    if (!parseSceneConfig(json, c, &err)) {
        std::cerr << "  " << err << "\n";
        ++lqdtest::failures();
        return r;
    }
    r.w = std::make_unique<World>(lib(), budgetForTier(c.tier), c.seed);
    if (!buildScene(*r.w, c, &err)) {
        std::cerr << "  " << err << "\n";
        ++lqdtest::failures();
        return r;
    }
    for (double t = chunk; t <= until + 1e-9; t += chunk) r.w->advance(t);
    r.ok = true;
    return r;
}
bool logged(const World& w, const std::string& what) {
    for (const auto& e : w.events)
        if (e.text.find(what) != std::string::npos) return true;
    return false;
}
void conserved(const World& w) {
    const VolumeBook v = w.volumes();
    CHECK_NEAR(v.total(), v.initial, 1e-9 * std::max(v.initial, 1e-6) + 1e-12);
}
}  // namespace

TEST_CASE("liquid: water at rest in a tank stays at rest density and settles") {
    World w(lib(), budgetForTier(1), 1);
    const SurfaceDef* s = lib().surface("glass_pane");
    w.addPlane(Role::Floor, "floor", {0, 0, 0}, {0, 1, 0}, s);
    w.addPlane(Role::Wall, "x0", {0, 0, 0}, {1, 0, 0}, s);
    w.addPlane(Role::Wall, "x1", {0.1, 0, 0}, {-1, 0, 0}, s);
    w.addPlane(Role::Wall, "z0", {0, 0, 0}, {0, 0, 1}, s);
    w.addPlane(Role::Wall, "z1", {0, 0, 0.1}, {0, 0, -1}, s);
    w.addLiquidBox(*lib().liquid("water"), {0, 0, 0}, {0.1, 0.1, 0.1}, {});
    w.advance(lqdtest::kOptimisedBuild ? 1.5 : 0.4);
    CHECK(w.fluid.maxCompression() < 0.03);
    CHECK(w.fluid.maxSpeed() < 0.1);
    std::cerr << "  tank: compression " << w.fluid.maxCompression() << ", speed " << w.fluid.maxSpeed() << "\n";
    // The free surface is level to within a particle: the top of every 2 cm
    // column inside the walls.
    double hi = 0, lo = 1;
    for (int cx = 1; cx < 4; ++cx)
        for (int cz = 1; cz < 4; ++cz) {
            double top = 0;
            for (const auto& p : w.fluid.x)
                if (p.x >= 0.02 * cx && p.x < 0.02 * (cx + 1) && p.z >= 0.02 * cz && p.z < 0.02 * (cz + 1)) top = std::max(top, p.y);
            hi = std::max(hi, top);
            lo = std::min(lo, top);
        }
    std::cerr << "  tank: free surface from " << lo << " to " << hi << "\n";
    CHECK(hi - lo < w.fluid.settings().spacing);
    // And the poured 10 cm stays 10 cm: the top particles' centres sit within
    // two spacings of it.
    CHECK_RANGE(hi, 0.1 - 2 * w.fluid.settings().spacing, 0.1);
    conserved(w);
}

TEST_CASE("liquid: dam-break front against Martin & Moyce (1952)") {
    // Surge front position Z = x/a against T = t sqrt(2g/a) for a column of width a
    // and height 2a. Martin & Moyce measured Z = 2.0 at T = 1.8 and Z = 3.0 at
    // T = 2.5 (their figure 3, the n^2 = 2 column); the shallow-water (Ritter)
    // front 1 + 2 T / sqrt 2 runs ahead of both.
    const Run r = run(R"({"scenario":"dam_break","tier":1})", 0.0);
    World& w = *r.w;
    const double a = 0.146, Ts = std::sqrt(a / (2 * kGravity));
    auto front = [&] {
        double x = 0;
        for (const auto& p : w.fluid.x)
            if (p.y < 0.03) x = std::max(x, p.x);
        return x / a;
    };
    w.advance(1.8 * Ts);
    const double z18 = front();
    w.advance(2.5 * Ts);
    const double z25 = front();
    std::cerr << "  dam break: Z(1.8) = " << z18 << ", Z(2.5) = " << z25 << "\n";
    CHECK_RANGE(z18, 1.7, 2.6);
    CHECK_RANGE(z25, 2.6, 3.6);
    CHECK(z25 < 1 + 2 * 2.5 / std::sqrt(2.0));
}

TEST_CASE("wine glass: knocked off the table it breaks on tile, holds on carpet; volume is conserved") {
    const Run tile = run(R"({"scenario":"wine_glass","tier":0,"floor":"ceramic_tile"})", lqdtest::kOptimisedBuild ? 1.6 : 0.9);
    CHECK(logged(*tile.w, "struck floor"));
    CHECK(tile.w->breaks() == 1);
    CHECK(logged(*tile.w, "shattered into"));
    CHECK(tile.w->wettedArea() > 0.02);
    CHECK(!tile.w->spray.stains.empty());
    conserved(*tile.w);
    const Run carpet = run(R"({"scenario":"wine_glass","tier":0,"floor":"carpet"})", lqdtest::kOptimisedBuild ? 1.6 : 0.9);
    CHECK(logged(*carpet.w, "struck floor"));
    CHECK(carpet.w->breaks() == 0);
    conserved(*carpet.w);
    // The glass's knock on the table itself never breaks it.
    for (const auto& e : tile.w->events)
        if (e.text.find("struck table") != std::string::npos) CHECK(e.text.find("holds") != std::string::npos);
}

TEST_CASE("bottle: dropped from a metre onto concrete it bursts, with the wine conserved") {
    const Run r = run(R"({"scenario":"bottle_drop","tier":0,"floor":"concrete"})", lqdtest::kOptimisedBuild ? 0.8 : 0.55);
    const World& w = *r.w;
    bool fast = false;
    for (const auto& e : w.events)
        if (e.text.find("struck floor at 4.") != std::string::npos) fast = true;  // sqrt(2 g h) = 4.4 m/s
    CHECK(fast);
    CHECK(w.breaks() == 1);
    int shards = 0;
    for (const auto& b : w.bodies) shards += b.shard >= 0 ? 1 : 0;
    CHECK(shards >= 20);
    conserved(w);
}

TEST_CASE("gunshot: a 9 mm round through the dummy, back and forward spatter, stains on the wall") {
    const Run r = run(R"({"scenario":"gunshot","tier":1})", lqdtest::kOptimisedBuild ? 2.5 : 1.2);
    const World& w = *r.w;
    CHECK(logged(w, "entry wound"));
    CHECK(logged(w, "exit wound"));
    const Projectile& p = w.projectiles[0];
    CHECK(p.state == Projectile::State::Stopped);
    CHECK(p.fate.find("back wall") != std::string::npos);
    // 9 mm FMJ keeps most of its speed through a 0.2 m torso (it penetrates 50-70 cm of gelatin).
    double back = 0, fwd = 0;
    for (const auto& s : w.spray.sources) {
        if (s.kind == "back spatter") back += s.volume;
        if (s.kind == "forward spatter") fwd += s.volume;
    }
    CHECK(fwd > 2 * back);
    CHECK(back > 0);
    int wall = 0, floor = 0;
    for (const auto& s : w.spray.stains) {
        const Role role = w.colliders[s.collider].role;
        wall += role == Role::Wall ? 1 : 0;
        floor += role == Role::Floor ? 1 : 0;
    }
    CHECK(wall > 10);
    CHECK(floor > 50);
    conserved(w);
}

TEST_CASE("gunshot: rounds through a pane lose speed; a pane behind the dummy keeps a spatter record") {
    const Run front = run(R"({"scenario":"gunshot","tier":1,"pane":"front"})", 0.05, 0.01);
    CHECK(logged(*front.w, "through the pane"));
    CHECK(logged(*front.w, "entry wound in the"));
    const Run behind = run(R"({"scenario":"gunshot","tier":1,"pane":"behind"})", lqdtest::kOptimisedBuild ? 1.5 : 0.8);
    const World& w = *behind.w;
    int pane = -1;
    for (const auto& c : w.colliders)
        if (c.role == Role::Pane) pane = c.id;
    CHECK(pane >= 0);
    CHECK(!w.colliders[pane].holes.empty());
    BpaOptions o;
    o.collider = pane;
    const BpaResult r = analyseStains(w.spray, w.colliders, o);
    std::cerr << "  pane behind: " << r.summary << "\n";
    CHECK(r.ok);
    // 0.7 m from the exit wound, the stains point back to it to within 0.15 m.
    CHECK(r.convergenceError < 0.15);
    CHECK(std::abs(r.heightError) < 0.2);
}

TEST_CASE("gunshot: a host can resolve the pane itself (Glass Dynamics hand-off)") {
    SceneConfig c;
    std::string err;
    CHECK(parseSceneConfig(R"({"scenario":"gunshot","tier":1,"pane":"front","defer_panes":true})", c, &err));
    World w(lib(), budgetForTier(1), 1);
    CHECK(buildScene(w, c, &err));
    w.advance(0.03);
    CHECK(w.pendingStrikes().size() == 1);
    const PaneStrike s = w.pendingStrikes()[0];
    CHECK_NEAR(s.speed, 360, 1);
    CHECK(!s.fromBack);
    CHECK_RANGE(s.x, 0.3, 0.7);
    PaneResult r;
    r.speedOut = 100;  // a laminated pane takes most of it
    r.holeRadius = 0.01;
    w.resolveStrike(0, r);
    w.advance(0.05);
    CHECK(logged(w, "through the pane: 360 -> 100"));
    // At 100 m/s the round stops in the dummy.
    CHECK(w.projectiles[0].fate == "lodged in the dummy");
}

TEST_CASE("drip: passive drops leave round stains that grow with the fall") {
    const Run r = run(R"({"scenario":"drip","tier":1,"height_m":2.1})", 1.2);
    const World& w = *r.w;
    std::vector<std::pair<double, double>> hd;  // fall height, stain width
    for (const auto& s : w.spray.stains) {
        if (s.satellite) continue;
        CHECK(s.width / s.length > 0.99);
        hd.push_back({w.spray.sources[s.source].origin.y, s.width});
    }
    CHECK(hd.size() == 7);
    std::sort(hd.begin(), hd.end());
    for (std::size_t i = 1; i < hd.size(); ++i) CHECK(hd[i].second >= hd[i - 1].second);
    // A 50 ul drop leaves a stain of 10-25 mm from 0.3-2 m.
    CHECK_RANGE(hd.front().second * 1e3, 8, 25);
    CHECK_RANGE(hd.back().second * 1e3, hd.front().second * 1e3, 30);
}

TEST_CASE("angles: the impact-angle card reads back its angles from width over length") {
    const Run r = run(R"({"scenario":"angles","tier":1})", 0.3);
    const World& w = *r.w;
    int n = 0;
    for (const auto& s : w.spray.stains) {
        if (s.satellite || w.colliders[s.collider].name != "card") continue;
        ++n;
        const double est = std::asin(std::min(1.0, s.width / s.length));
        CHECK_NEAR(est, s.alpha, 0.02);
    }
    CHECK(n == 9);
}

TEST_CASE("determinism: same seed, same stains; the result does not depend on the host's frame rate") {
    const Run a = run(R"({"scenario":"gunshot","tier":1,"seed":5})", 0.6, 1.0 / 60);
    const Run b = run(R"({"scenario":"gunshot","tier":1,"seed":5})", 0.6, 1.0 / 240);
    CHECK(a.w->spray.stains.size() == b.w->spray.stains.size());
    bool same = a.w->spray.stains.size() == b.w->spray.stains.size();
    for (std::size_t i = 0; same && i < a.w->spray.stains.size(); ++i)
        same = length(a.w->spray.stains[i].pos - b.w->spray.stains[i].pos) == 0;
    CHECK(same);
    const Run c = run(R"({"scenario":"gunshot","tier":1,"seed":6})", 0.6);
    CHECK(c.w->spray.stains.size() != a.w->spray.stains.size() ||
          length(c.w->spray.stains[0].pos - a.w->spray.stains[0].pos) > 0);
}

TEST_CASE("export: the glTF file is valid and holds the liquid, the glass and the stains") {
    SceneConfig c;
    std::string err;
    parseSceneConfig(R"({"scenario":"wine_glass","tier":0})", c, &err);
    World w(lib(), budgetForTier(0), 1);
    buildScene(w, c, &err);
    Recorder rec(30);
    for (double t = 0; t < 0.9; t += 1.0 / 30) {
        w.advance(t);
        rec.capture(w);
    }
    const std::vector<std::uint8_t> glb = exportGlb(w, rec);
    CHECK(glb.size() > 1000);
    auto u32 = [&](std::size_t o) {
        return static_cast<std::uint32_t>(glb[o]) | static_cast<std::uint32_t>(glb[o + 1]) << 8 |
               static_cast<std::uint32_t>(glb[o + 2]) << 16 | static_cast<std::uint32_t>(glb[o + 3]) << 24;
    };
    CHECK(u32(0) == 0x46546C67);
    CHECK(u32(4) == 2);
    CHECK(u32(8) == glb.size());
    const std::uint32_t jlen = u32(12);
    CHECK(u32(16) == 0x4E4F534A);
    const std::string js(glb.begin() + 20, glb.begin() + 20 + jlen);
    const json::ParseResult pr = json::parse(js);
    CHECK(pr.ok());
    if (pr.ok()) {
        const json::Value* nodes = pr.value.find("nodes");
        const json::Value* anim = pr.value.find("animations");
        CHECK(nodes && nodes->items().size() > rec.frames.size());
        CHECK(anim != nullptr);
        bool liquid = false, glass = false, stains = false;
        for (const auto& n : nodes->items()) {
            const std::string name = n.find("name")->asString();
            liquid = liquid || name == "liquid";
            glass = glass || name.find("Bordeaux wine glass") == 0;
            stains = stains || name == "stains";
        }
        CHECK(liquid);
        CHECK(glass);
        CHECK(stains == !w.spray.stains.empty());
        const std::uint32_t blen = u32(20 + jlen);
        CHECK(pr.value.find("buffers")->items()[0].find("byteLength")->asNumber() == blen);
    }
    const auto ply = exportPly(rec.frames.back(), rec.spacing);
    CHECK(std::string(ply.begin(), ply.begin() + 3) == "ply");
    const json::ParseResult sj = json::parse(exportStainsJson(w));
    CHECK(sj.ok());
}

TEST_CASE("ragdoll: shot, the dummy comes off its stand and falls, carrying its wounds") {
    const Run r = run(R"({"scenario":"gunshot","tier":1})", lqdtest::kOptimisedBuild ? 2.5 : 1.2);
    const World& w = *r.w;
    CHECK(w.ragdoll.active());
    CHECK_NEAR(w.ragdoll.mass(), 78.2, 0.01);  // Hybrid III 50th percentile
    CHECK(w.ragdoll.released());
    CHECK(logged(w, "the stand lets go"));
    // The head is near the floor, and the wounds went down with the chest.
    int head = -1;
    for (int s = 0; s < w.ragdoll.segments(); ++s)
        if (w.colliders[w.dummyCollider(s)].name == "head") head = s;
    CHECK(head >= 0);
    if (head >= 0) CHECK(w.ragdoll.a(head).y < 0.5);
    for (const auto& wd : w.wounds) CHECK(wd.pos.y < 0.9);
    // Nothing went through the floor.
    for (int s = 0; s < w.ragdoll.segments(); ++s) CHECK(std::min(w.ragdoll.a(s).y, w.ragdoll.b(s).y) > -0.01);
    conserved(w);
}

TEST_CASE("buckshot: nine pellets in a pattern, each its own wound") {
    const Run r = run(R"({"scenario":"gunshot","tier":1,"projectile":"buckshot_00","ragdoll":false})", 0.1, 0.01);
    const World& w = *r.w;
    CHECK(w.projectiles.size() == 9);
    int entries = 0;
    for (const auto& wd : w.wounds) entries += wd.exit ? 0 : 1;
    CHECK(entries >= 7);
    // A cylinder-bore pattern: within ~10 cm of the aim at 3 m.
    double spread = 0;
    for (const auto& p : w.projectiles) spread = std::max(spread, length(cross(p.dir, w.projectiles[0].dir)) * 3.0);
    CHECK(spread > 0.002);
    CHECK(spread < 0.12);
}

TEST_CASE("spill: most of the wine pools, the rest stays where it lands") {
    const Run r = run(R"({"scenario":"wine_glass","tier":1})", lqdtest::kOptimisedBuild ? 2.0 : 1.2);
    const World& w = *r.w;
    const VolumeBook v = w.volumes();
    CHECK(v.bulk > 0.7 * v.initial);
    CHECK(v.stained > 0.02 * v.initial);
    conserved(w);
}
