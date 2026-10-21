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
    // The worst single particle, so it moves with the compiler's rounding
    // (2.6% with gcc, 3.1% with MSVC); the solver's tolerance is a few percent.
    CHECK(w.fluid.maxCompression() < 0.035);
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
    const Run r = run(R"({"scenario":"gunshot","tier":1,"ragdoll":true})", lqdtest::kOptimisedBuild ? 2.5 : 1.2);
    const World& w = *r.w;
    CHECK(w.ragdoll.active());
    CHECK_NEAR(w.ragdoll.mass(), 78.2, 0.01);  // Hybrid III 50th percentile
    CHECK(w.ragdoll.released());
    CHECK(logged(w, "the stand lets go"));
    // A 9 mm round's momentum barely moves a body: it crumples where it stood,
    // the wounds going down with the chest.
    int head = -1, pelvis = -1;
    for (int s = 0; s < w.ragdoll.segments(); ++s) {
        if (w.colliders[w.dummyCollider(s)].name == "head") head = s;
        if (w.colliders[w.dummyCollider(s)].name == "pelvis") pelvis = s;
    }
    CHECK(head >= 0 && pelvis >= 0);
    if (head >= 0 && pelvis >= 0) {
        CHECK(w.ragdoll.a(head).y < 1.0);
        CHECK(w.ragdoll.a(pelvis).y < 0.45);
        CHECK(std::abs(w.ragdoll.a(pelvis).z) < 0.5);
    }
    for (const auto& wd : w.wounds) CHECK(wd.pos.y < 1.1);
    // Nothing went through the floor or flew off.
    for (int s = 0; s < w.ragdoll.segments(); ++s) {
        CHECK(std::min(w.ragdoll.a(s).y, w.ragdoll.b(s).y) > -0.01);
        CHECK(std::max(w.ragdoll.a(s).y, w.ragdoll.b(s).y) < 1.8);
        CHECK(std::hypot(w.ragdoll.a(s).x, w.ragdoll.a(s).z) < 2.0);
    }
    conserved(w);
}

TEST_CASE("ragdoll: with knockback, a load of buckshot puts the dummy on its back, away from the gun") {
    const Run r = run(R"({"scenario":"gunshot","tier":1,"ragdoll":true,"knockback":0.05,"projectile":"buckshot_00"})",
                      lqdtest::kOptimisedBuild ? 2.0 : 1.2);
    const World& w = *r.w;
    int head = -1;
    for (int s = 0; s < w.ragdoll.segments(); ++s)
        if (w.colliders[w.dummyCollider(s)].name == "head") head = s;
    CHECK(head >= 0);
    if (head < 0) return;
    // The gun is at +z: down, and over backwards.
    CHECK(w.ragdoll.a(head).y < 0.5);
    CHECK(w.ragdoll.a(head).z < -0.3);
    for (int s = 0; s < w.ragdoll.segments(); ++s) CHECK(std::min(w.ragdoll.a(s).y, w.ragdoll.b(s).y) > -0.01);
}

TEST_CASE("ragdoll: the forward spatter flies on behind a dummy that falls, as behind one that stands") {
    auto behind = [](const char* json) {
        const Run r = run(json, 0.6, 0.02);
        int wall = 0, own = 0;
        for (const auto& st : r.w->spray.stains) {
            if (r.w->spray.sources[st.source].kind != "forward spatter") continue;
            wall += r.w->colliders[st.collider].name == "back wall" ? 1 : 0;
            // None may land on the body it left the moment it leaves.
            own += r.w->colliders[st.collider].role == Role::Dummy && st.time - r.w->spray.sources[st.source].time < 1e-3 ? 1 : 0;
        }
        return std::make_pair(wall, own);
    };
    const auto [fixedWall, fixedOwn] = behind(R"({"scenario":"gunshot","tier":1,"spatter_boost":3})");
    const auto [ragWall, ragOwn] = behind(R"({"scenario":"gunshot","tier":1,"spatter_boost":3,"ragdoll":true,"knockback":0.05})");
    std::cerr << "  back wall: " << fixedWall << " stains behind a standing dummy, " << ragWall << " behind a ragdoll\n";
    CHECK(fixedOwn == 0);
    CHECK(ragOwn == 0);
    CHECK(ragWall > 0.7 * fixedWall);
}

TEST_CASE("ragdoll: off by default, the dummy stands where it is") {
    const Run r = run(R"({"scenario":"gunshot","tier":1})", 0.3, 0.01);
    const World& w = *r.w;
    CHECK(!w.ragdoll.active());
    CHECK(!logged(w, "the stand lets go"));
    CHECK(logged(w, "exit wound"));
    for (const auto& c : w.colliders)
        if (c.name == "head") CHECK(c.a.y > 1.4);
}

TEST_CASE("spatter boost: more than life on request, three times the spatter in bigger drops") {
    auto spattered = [](const World& w) {
        double v = 0, d = 0;
        int n = 0;
        for (const auto& src : w.spray.sources)
            if (src.kind == "back spatter" || src.kind == "forward spatter") v += src.volume;
        for (const auto& st : w.spray.stains)
            if (!st.satellite) {
                d += st.diameter;
                ++n;
            }
        return std::pair<double, double>{v, n ? d / n : 0};
    };
    const Run a = run(R"({"scenario":"gunshot","tier":1})", 0.25, 0.01);
    const Run b = run(R"({"scenario":"gunshot","tier":1,"spatter_boost":3})", 0.25, 0.01);
    const auto [va, da] = spattered(*a.w);
    const auto [vb, db] = spattered(*b.w);
    std::cerr << "  boost 3: " << vb * 1e6 << " ml against " << va * 1e6 << " ml; mean stain drop " << db * 1e3 << " mm against " << da * 1e3 << " mm\n";
    CHECK_RANGE(vb / va, 2.5, 3.5);
    CHECK(db > da);
    conserved(*b.w);
    SceneConfig c;
    std::string err;
    CHECK(!parseSceneConfig(R"({"scenario":"gunshot","spatter_boost":12})", c, &err));
}

TEST_CASE("spatter boost: at 11 the wounds spurt across the room, then stop bleeding") {
    const Run r = run(R"({"scenario":"gunshot","tier":1,"spatter_boost":11})", lqdtest::kOptimisedBuild ? 4.0 : 2.5, 0.02);
    const World& w = *r.w;
    int spurts = 0;
    double far = 0;
    for (const auto& src : w.spray.sources)
        if (src.kind == "arterial spurt") ++spurts;
    for (const auto& st : w.spray.stains) {
        const auto& src = w.spray.sources[st.source];
        if (src.kind == "arterial spurt") far = std::max(far, std::hypot(st.pos.x - src.origin.x, st.pos.z - src.origin.z));
    }
    std::cerr << "  boost 11: " << spurts << " spurts, landing up to " << far << " m away\n";
    CHECK(spurts == 6);  // three from each wound
    CHECK(far > 1.0);
    // The wounds ooze for a second or two: each bleeds no more than its rate
    // times the 1.5 s it decays over (0.8 ml/s, times the square root of 11).
    for (const auto& wd : w.wounds) {
        CHECK(wd.bled > 0.5e-6);
        CHECK(wd.bled < 0.8e-6 * std::sqrt(11.0) * 1.5 * 1.05);
    }
    conserved(w);
}

TEST_CASE("settling: fine drops at their terminal speed land where they would have, at once") {
    auto rain = [](bool settle) {
        World w(lib(), budgetForTier(1), 7);
        w.addPlane(Role::Floor, "floor", {0, 0, 0}, {0, 1, 0}, lib().surface("ceramic_tile"));
        w.addBox(Role::Table, "table", {0.3, 0.7, 0}, {0.2, 0.02, 0.2}, Quat{}, lib().surface("oak_table"));
        w.setLiquid(*lib().liquid("stage_blood_green"));
        w.spray.settle = settle;
        for (int i = 0; i < 40; ++i) w.drop({0.03 * i - 0.3, 1.5, 0.01 * (i % 7)}, {0.4, 0, -0.2}, 3e-4, "mist");
        w.advance(4.0);
        return std::make_pair(w.spray.stains, w.spray.settled);
    };
    const auto [flown, none] = rain(false);
    const auto [settled, n] = rain(true);
    std::cerr << "  settling: " << n << " of 40 drops settled\n";
    CHECK(none == 0);
    CHECK(n > 30);
    CHECK(flown.size() == settled.size());
    double worst = 0;
    for (std::size_t i = 0; i < std::min(flown.size(), settled.size()); ++i) {
        // Same source order: match each drop's stain by its collider and nearest position.
        double best = 1e30;
        for (const auto& s : settled)
            if (s.collider == flown[i].collider) best = std::min(best, length(s.pos - flown[i].pos));
        worst = std::max(worst, best);
    }
    std::cerr << "  settling: worst miss " << worst * 1e3 << " mm\n";
    CHECK(worst < 3e-3);  // within 3 mm, on the table or the floor as before
}

TEST_CASE("buckshot: nine pellets in a pattern, each its own wound") {
    const Run r = run(R"({"scenario":"gunshot","tier":1,"projectile":"buckshot_00"})", 0.1, 0.01);
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
    // Every wound throws its spatter, however full the drop budget: past a
    // share of it the drops are sampled, each standing for several.
    int sprays = 0, sampled = 0;
    for (const auto& src : w.spray.sources)
        if (src.kind == "back spatter" || src.kind == "forward spatter") {
            ++sprays;
            CHECK(src.drops > 0);
        }
    for (const auto& d : w.spray.drops) sampled += d.weight > 1.01 ? 1 : 0;
    CHECK(sprays >= entries);
    CHECK(sampled > 0);
    std::cerr << "  buckshot: " << sprays << " sprays, " << w.spray.drops.size() << " drops in the air of a budget of "
              << w.spray.maxDrops << ", " << sampled << " of them sampled\n";
    // A full budget is overrun by at most a few dozen drops a wound.
    CHECK(static_cast<int>(w.spray.drops.size()) <= w.spray.maxDrops + 48 * sprays);
    conserved(w);
}

TEST_CASE("pane pieces: back spatter lands on the pieces of a broken pane and falls with them") {
    SceneConfig c;
    std::string err;
    CHECK(parseSceneConfig(R"({"scenario":"gunshot","tier":1,"pane":"front","pane_gap_m":0.5,"defer_panes":true})", c, &err));
    World w(lib(), budgetForTier(1), 1);
    CHECK(buildScene(w, c, &err));
    w.advance(0.03);
    CHECK(w.pendingStrikes().size() == 1);
    PaneResult res;
    res.speedOut = 330;
    res.holeRadius = 0.006;
    w.resolveStrike(0, res);
    int pane = -1;
    for (const auto& k : w.colliders)
        if (k.role == Role::Pane) pane = k.id;
    CHECK(pane >= 0);
    if (pane < 0) return;
    // The host's Glass Dynamics: the pane cracked down the middle into two
    // pieces, both still in the frame (Glass Dynamics' frame: from the lower
    // left corner, z into the glass from the front face).
    const Collider P = w.colliders[pane];
    const double W = 2 * P.half.x, H = 2 * P.half.y, T = 2 * P.half.z, cell = 0.02;
    const int nx = static_cast<int>(std::ceil(W / cell)), ny = static_cast<int>(std::ceil(H / cell));
    std::vector<int> label(static_cast<std::size_t>(nx) * ny);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) label[static_cast<std::size_t>(j) * nx + i] = (i + 0.5) * cell < W / 2 ? 0 : 1;
    w.setPaneMap(pane, nx, ny, cell, label, label, false);
    const std::vector<double> half{-W / 4, -H / 2, W / 4, -H / 2, W / 4, H / 2, -W / 4, H / 2};
    for (int k = 0; k < 2; ++k) {
        w.addPanePiece(pane, k, half, T, false);
        w.setPanePiecePose(pane, k, {W / 4 + k * W / 2, H / 2, T / 2}, Quat{}, {}, 0);
    }
    // Back spatter flies back to the pane, half a metre in front of the
    // dummy, onto the face towards it (a windscreen seen from inside).
    w.advance(0.35);
    std::vector<std::pair<std::size_t, Vec3>> onLeft, onRight;
    for (const auto& ps : w.pieceStains()) {
        const Stain& st = w.spray.stains[ps.stain];
        CHECK(ps.piece == 0 || ps.piece == 1);
        CHECK(st.normal.z < -0.9);  // the inside face
        (ps.piece == 0 ? onLeft : onRight).push_back({ps.stain, st.pos});
    }
    std::cerr << "  pane pieces: " << onLeft.size() << " stains on the left piece, " << onRight.size() << " on the right\n";
    CHECK(onLeft.size() + onRight.size() > 5);
    // The left piece falls out of the frame and drops 20 cm, tipping over;
    // the host updates it every 10 ms.
    const double t0 = w.time();
    for (int k = 1; k <= 20; ++k) {
        const double t = 0.01 * k, y = H / 2 - 0.5 * kGravity * t * t;
        const Quat tip = Quat::axisAngle({1, 0, 0}, 0.6 * t);
        w.setPanePiecePose(pane, 0, {W / 4, std::max(y, H / 2 - 0.2), T / 2 + 0.3 * t}, tip, {0, -kGravity * t, 0.3}, y > H / 2 - 0.2 ? 2 : 3);
        w.advance(t0 + t);
    }
    const PanePieces* bp = w.panePieces(pane);
    CHECK(bp && bp->pieces[0].collider >= 0);
    if (!bp || bp->pieces[0].collider < 0) return;
    const Collider& slab = w.colliders[bp->pieces[0].collider];
    CHECK(slab.kind == ShapeKind::Slab);
    double worst = 0, fell = 1;
    for (const auto& [k, at] : onLeft) {
        const Stain& st = w.spray.stains[k];
        CHECK(st.carried);
        worst = std::max(worst, std::abs(slab.sdf(st.pos)));
        fell = std::min(fell, at.y - st.pos.y);
    }
    if (!onLeft.empty()) {
        CHECK(worst < 1e-6);  // still on the glass
        CHECK(fell > 0.1);    // and down with it
    }
    for (const auto& [k, at] : onRight) CHECK(length(w.spray.stains[k].pos - at) < 1e-12);
    // Where the left piece was there is now no glass; the right piece is still there.
    CHECK(w.colliders[pane].inHole({-W / 4, 0, 0}));
    CHECK(!w.colliders[pane].inHole({W / 4 + 0.05, 0.2, 0}));
    CHECK(!w.takeMovedStains().empty() == !onLeft.empty());
    CHECK(w.takeMovedStains().empty());
    // It crumbles where it lies: its stains go with it.
    w.setPanePiecePose(pane, 0, {W / 4, H / 2 - 0.2, T / 2 + 0.06}, Quat{}, {}, -1);
    w.advance(w.time() + 0.01);
    for (const auto& [k, at] : onLeft) CHECK(w.spray.stains[k].gone);
    conserved(w);
}

TEST_CASE("pane pieces: the glTF holds each piece, keyed every frame, with its stains on it") {
    SceneConfig c;
    std::string err;
    CHECK(parseSceneConfig(R"({"scenario":"gunshot","tier":0,"pane":"front","defer_panes":true})", c, &err));
    World w(lib(), budgetForTier(0), 1);
    CHECK(buildScene(w, c, &err));
    Recorder rec(60);
    int pane = -1;
    for (const auto& k : w.colliders)
        if (k.role == Role::Pane) pane = k.id;
    const Collider P = w.colliders[pane];
    const double W = 2 * P.half.x, H = 2 * P.half.y, T = 2 * P.half.z;
    bool broke = false;
    for (double t = 1.0 / 60; t < 0.25; t += 1.0 / 60) {
        w.advance(t);
        if (!w.pendingStrikes().empty()) {
            w.resolveStrike(0, PaneResult{true, 330, 0.006, false});
            // One piece with the bullet hole through it, and the rest.
            std::vector<int> label(10 * 10, 1);
            for (int j = 3; j < 7; ++j)
                for (int i = 3; i < 7; ++i) label[j * 10 + i] = 0;
            w.setPaneMap(pane, 10, 10, std::max(W, H) / 10, label, label, false);
            w.addPanePiece(pane, 0, {-0.1, -0.1, 0.1, -0.1, 0.1, 0.1, -0.1, 0.1}, T, false, {{-0.01, -0.01, -0.01, 0.01, 0.01, 0.01, 0.01, -0.01}});
            w.addPanePiece(pane, 1, {-W / 2, -H / 2, W / 2, -H / 2, W / 2, H / 2, -W / 2, H / 2}, T, false);
            broke = true;
        }
        if (broke) {
            w.setPanePiecePose(pane, 0, {W / 2, H / 2 - 0.5 * (t > 0.1 ? t - 0.1 : 0), T / 2}, Quat{}, {}, t > 0.1 ? 2 : 0);
            w.setPanePiecePose(pane, 1, {W / 2, H / 2, T / 2}, Quat{}, {}, 0);
        }
        rec.capture(w);
    }
    CHECK(broke);
    const std::vector<std::uint8_t> glb = exportGlb(w, rec);
    const std::uint32_t jlen = static_cast<std::uint32_t>(glb[12]) | static_cast<std::uint32_t>(glb[13]) << 8 |
                               static_cast<std::uint32_t>(glb[14]) << 16 | static_cast<std::uint32_t>(glb[15]) << 24;
    const json::ParseResult pr = json::parse(std::string(glb.begin() + 20, glb.begin() + 20 + jlen));
    CHECK(pr.ok());
    if (!pr.ok()) return;
    bool piece0 = false, box = false;
    for (const auto& n : pr.value.find("nodes")->items()) {
        const std::string name = n.find("name")->asString();
        piece0 = piece0 || name == "pane piece 0";
        box = box || name == "pane";
    }
    CHECK(piece0);
    CHECK(!box);  // the intact pane is drawn as its pieces
}

TEST_CASE("spill: most of the wine pools, the rest stays where it lands") {
    const Run r = run(R"({"scenario":"wine_glass","tier":1})", lqdtest::kOptimisedBuild ? 2.0 : 1.2);
    const World& w = *r.w;
    const VolumeBook v = w.volumes();
    CHECK(v.bulk > 0.7 * v.initial);
    CHECK(v.stained > 0.02 * v.initial);
    conserved(w);
}
