// Liquid Dynamics — lqdemo: headless runs, exports for game engines, and the validation table.
// SPDX-License-Identifier: MIT
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "lqd/bpa.hpp"
#include "lqd/export.hpp"
#include "lqd/library.hpp"
#include "lqd/rigid.hpp"
#include "lqd/scenarios.hpp"
#include "lqd/spray.hpp"
#include "lqd/world.hpp"

using namespace lqd;

namespace {

bool build(const Library& lib, const std::string& json, std::unique_ptr<World>& w, SceneConfig& c) {
    std::string err;
    if (!parseSceneConfig(json, c, &err)) {
        std::cerr << err << "\n";
        return false;
    }
    w = std::make_unique<World>(lib, budgetForTier(c.tier), c.seed);
    if (!buildScene(*w, c, &err)) {
        std::cerr << err << "\n";
        return false;
    }
    return true;
}

void analyse(const World& w) {
    for (const auto& col : w.colliders) {
        int n = 0;
        for (const auto& st : w.spray.stains) n += st.collider == col.id ? 1 : 0;
        if (n < 5 || col.role == Role::Dummy || col.role == Role::Impactor) continue;
        BpaOptions o;
        o.collider = col.id;
        o.readError = 0.25e-3;
        const BpaResult r = analyseStains(w.spray, w.colliders, o);
        std::printf("stain analysis, %s (%d stains): %s\n", col.name.c_str(), n, r.summary.c_str());
    }
}

int run(const Library& lib, int argc, char** argv) {
    const std::string json = argc > 2 ? argv[2] : "{}";
    const double until = argc > 3 ? std::atof(argv[3]) : 2.0;
    std::unique_ptr<World> wp;
    SceneConfig c;
    if (!build(lib, json, wp, c)) return 1;
    World& w = *wp;
    const auto t0 = std::chrono::steady_clock::now();
    std::size_t shown = 0;
    for (double t = 0; t < until; t += 0.05) {
        w.advance(t + 0.05);
        for (; shown < w.events.size(); ++shown) std::printf("%8.4f s  %s\n", w.events[shown].time, w.events[shown].text.c_str());
        const VolumeBook v = w.volumes();
        std::printf("%8.3f s  steps %6d  bulk %5zu (%.1f ml)  drops %5zu  stains %5zu  wet %.3f m2  err %.3f  vmax %.2f  dt %.1e\n",
                    w.time(), w.steps(), w.fluid.size(), v.bulk * 1e6, w.spray.drops.size(), w.spray.stains.size(),
                    w.wettedArea(), w.fluid.maxCompression(), w.fluid.maxSpeed(), w.nextStep());
        if (w.quiet() && t > 0.2) break;
    }
    for (; shown < w.events.size(); ++shown) std::printf("%8.4f s  %s\n", w.events[shown].time, w.events[shown].text.c_str());
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const VolumeBook v = w.volumes();
    std::printf("volume: in %.3f ml = bulk %.3f + air %.3f + stains %.3f + absorbed %.3f + lost %.3f (%.3f)\n",
                v.initial * 1e6, v.bulk * 1e6, v.air * 1e6, v.stained * 1e6, v.absorbed * 1e6, v.lost * 1e6, v.total() * 1e6);
    if (w.ragdoll.active() || !w.impactors.empty()) {
        const DummyEnergy e = w.dummyEnergy();
        const EnergyBook& b = e.book;
        std::printf("dummy energy: started %.0f J; given %.0f J (pushes %.0f, impactors %.0f, rounds into tissue %.0f); "
                    "absorbed %.0f J (in the blows %.0f, afterwards %.0f, tissue %.0f); moving it now %.0f J kinetic + %.0f J "
                    "potential; peak kinetic %.0f J; %d bursts, %.1f ml\n",
                    b.start, e.given(), b.fromHits, b.fromImpactors, e.inTissue, e.absorbed(), b.inBlows, b.afterwards,
                    e.inTissue, e.kinetic, e.potential, b.peakKinetic, e.bursts, e.burstVolume * 1e6);
    }
    std::printf("%.2f s of simulation in %.2f s\n", w.time(), secs);
    analyse(w);
    return 0;
}

// lqdemo export '<scene json>' seconds fps out-dir [--ply] [--no-scene]
int exportScene(const Library& lib, int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "lqdemo export '<scene json>' <seconds> <fps> <out-dir> [--ply] [--no-scene]\n";
        return 1;
    }
    const double until = std::atof(argv[3]), fps = std::atof(argv[4]);
    const std::filesystem::path dir = argv[5];
    bool ply = false, scene = true;
    for (int i = 6; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--ply")) ply = true;
        if (!std::strcmp(argv[i], "--no-scene")) scene = false;
    }
    std::unique_ptr<World> wp;
    SceneConfig c;
    if (!build(lib, argv[2], wp, c)) return 1;
    World& w = *wp;
    std::filesystem::create_directories(dir);
    Recorder rec(fps, static_cast<std::size_t>(std::ceil(until * fps)) + 2);
    rec.capture(w);
    const auto t0 = std::chrono::steady_clock::now();
    for (int k = 1; k <= static_cast<int>(std::ceil(until * fps)); ++k) {
        w.advance(k / fps);
        rec.capture(w);
        if (ply) {
            char name[64];
            std::snprintf(name, sizeof name, "points_%05d.ply", k);
            const auto bytes = exportPly(rec.frames.back(), w.fluid.settings().spacing);
            std::ofstream(dir / name, std::ios::binary).write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        }
        if (k % static_cast<int>(std::max(1.0, fps / 4)) == 0)
            std::fprintf(stderr, "\r%.2f / %.2f s  (%zu particles, %zu drops, %zu stains)   ", w.time(), until, w.fluid.size(),
                         w.spray.drops.size(), w.spray.stains.size());
    }
    std::fprintf(stderr, "\n");
    ExportOptions o;
    o.scene = scene;
    const auto glb = exportGlb(w, rec, o);
    std::ofstream(dir / (c.scenario + ".glb"), std::ios::binary).write(reinterpret_cast<const char*>(glb.data()), static_cast<std::streamsize>(glb.size()));
    std::ofstream(dir / (c.scenario + "_stains.json")) << exportStainsJson(w);
    std::ofstream log(dir / (c.scenario + "_events.txt"));
    for (const auto& e : w.events) log << e.time << "\t" << e.text << "\n";
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("%zu frames at %g fps, %.1f MB glTF, %zu stains -> %s (%.1f s)\n", rec.frames.size(), fps, static_cast<double>(glb.size()) / 1e6,
                w.spray.stains.size(), dir.string().c_str(), secs);
    return 0;
}

// The comparisons with published data, as a table (docs/validation.md).
int validate(const Library& lib) {
    std::printf("| Case | Reference | Model | Source |\n|---|---|---|---|\n");
    // Drops in air.
    const double gk[3][2] = {{0.5e-3, 2.06}, {1e-3, 4.03}, {2e-3, 6.49}};
    for (const auto& r : gk)
        std::printf("| Terminal speed, %.1f mm water drop | %.2f m/s | %.2f m/s | Gunn & Kinzer (1949) |\n", r[0] * 1e3, r[1],
                    terminalSpeed(r[0], 998.2));
    // Blood.
    const LiquidDef& b = *lib.liquid("stage_blood_green");
    std::printf("| Blood viscosity at 1 / 100 / 1000 1/s | Carreau fit | %.1f / %.2f / %.2f mPa s | Cho & Kensey (1991) |\n",
                b.viscosityAt(1) * 1e3, b.viscosityAt(100) * 1e3, b.viscosityAt(1000) * 1e3);
    // Dam break.
    {
        std::unique_ptr<World> w;
        SceneConfig c;
        build(lib, R"({"scenario":"dam_break","tier":1})", w, c);
        const double a = 0.146, Ts = std::sqrt(a / (2 * kGravity));
        auto front = [&] {
            double x = 0;
            for (const auto& p : w->fluid.x)
                if (p.y < 0.03) x = std::max(x, p.x);
            return x / a;
        };
        w->advance(1.8 * Ts);
        const double z18 = front();
        w->advance(2.5 * Ts);
        const double z25 = front();
        std::printf("| Dam-break surge front Z at T = 1.8 / 2.5 | 2.0 / 3.0 | %.2f / %.2f | Martin & Moyce (1952) |\n", z18, z25);
    }
    // Penetration.
    {
        const DummyDef& d = *lib.dummy("crash_dummy");
        for (const char* id : {"pistol_9mm_fmj", "magnum_44", "rifle_762"}) {
            const ProjectileDef& p = *lib.projectile(id);
            const double A = kPi / 4 * p.diameter * p.diameter, aa = d.tissueStrength * A, bb = 0.5 * d.tissueDensity * p.cdTissue * A;
            std::printf("| Penetration in 10%% gelatin, %s | %s | %.2f m | %s |\n", p.name.c_str(),
                        std::strcmp(id, "pistol_9mm_fmj") == 0 ? "0.60-0.75 m" : "-",
                        p.mass / (2 * bb) * std::log(1 + bb * p.speed * p.speed / aa),
                        std::strcmp(id, "pistol_9mm_fmj") == 0 ? "FBI-protocol tests (calibration point)" : "model");
        }
    }
    // Vessels.
    for (const auto& v : lib.vessels) {
        const VesselGeometry g(v);
        std::printf("| %s: mass, brimful capacity | %s | %.0f g, %.0f ml | %s |\n", v.name.c_str(),
                    v.id == "bordeaux_bottle" ? "400-550 g, 765-775 ml" : "140-220 g, -", g.mass() * 1e3,
                    g.cavityCapacity() * 1e6, v.id == "bordeaux_bottle" ? "750 ml bottle norms" : "machine-made stemware");
    }
    // Strikes.
    {
        const VesselDef& v = *lib.vessel("bordeaux_glass");
        const VesselGeometry g(v);
        WallSample bowl;
        for (const auto& w : g.samples())
            if (!w.edge && std::abs(w.p.y - 0.165) < 0.004 && w.p.x > 0) {
                bowl = w;
                break;
            }
        const struct {
            const char* s;
            double v;
            const char* expect;
        } cases[] = {{"oak_table", 0.3, "never breaks"}, {"oak_table", 1.0, "sometimes"}, {"ceramic_tile", 3.8, "breaks"},
                     {"carpet", 3.8, "rarely"}};
        for (const auto& k : cases) {
            const ImpactAssessment a = assessImpact(v, bowl, *lib.surface(k.s), 0.08, k.v, 1);
            std::printf("| Wine glass bowl at %.1f m/s on %s | %s | P(break) %.3f (%.0f MPa) | calibration point |\n", k.v,
                        lib.surface(k.s)->name.c_str(), k.expect, a.probability, a.stress * 1e-6);
        }
    }
    // Stain analysis behind the dummy.
    {
        std::unique_ptr<World> w;
        SceneConfig c;
        build(lib, R"({"scenario":"gunshot","tier":2,"pane":"behind"})", w, c);
        w->advance(1.5);
        int pane = -1, floor = -1;
        for (const auto& col : w->colliders) {
            if (col.role == Role::Pane) pane = col.id;
            if (col.role == Role::Floor) floor = col.id;
        }
        BpaOptions o;
        o.readError = 0.25e-3;
        o.collider = pane;
        const BpaResult rp = analyseStains(w->spray, w->colliders, o);
        o.collider = floor;
        const BpaResult rf = analyseStains(w->spray, w->colliders, o);
        std::printf("| Area of origin from a pane 0.7 m behind the exit wound | %.2f m | %.2f m (convergence off by %.0f mm) | tangent method |\n",
                    rp.height - rp.heightError, rp.height, rp.convergenceError * 1e3);
        std::printf("| Area of origin from the floor below a 1.3 m wound | %.2f m | %.2f m | tangent method: overestimates (Carter 2001) |\n",
                    rf.height - rf.heightError, rf.height);
    }
    // Blunt impacts and falls.
    {
        std::unique_ptr<World> w;
        SceneConfig c;
        build(lib, R"({"scenario":"impact","tier":0,"impactor":"car_50"})", w, c);
        const double M = w->impactors[0].def->mass, m = w->ragdoll.mass(), v = w->impactors[0].speed0;
        w->advance(0.1);
        const double first = w->dummyEnergy().book.fromImpactors;
        w->advance(4.0);
        double z = 0;
        for (int s = 0; s < w->ragdoll.segments(); ++s)
            if (w->ragdoll.segmentName(s) == "pelvis") z = -w->ragdoll.a(s).z;
        std::printf("| Energy a car at 50 km/h loses to the dummy in the first 0.1 s | %.1f kJ (inelastic, M m v^2 / (M + m)) | %.1f kJ | momentum conservation |\n",
                    M * m * v * v / (M + m) * 1e-3, first * 1e-3);
        std::printf("| Throw distance, car at 50 km/h braking at 0.75 g | 10-20 m | %.1f m | pedestrian throw distances (Wood 1988; Searle 1993) |\n", z);
    }
    {
        std::unique_ptr<World> w;
        SceneConfig c;
        build(lib, R"({"scenario":"fall","tier":0,"fall_m":10})", w, c);
        w->advance(2.0);
        const DummyEnergy e = w->dummyEnergy();
        std::printf("| Speed at landing, a fall from 10 m (centre of mass 10.5 m above where it stops) | %.1f m/s | %.1f m/s | free fall, less air drag (vt 55 m/s) |\n",
                    std::sqrt(2 * kGravity * 10.5), std::sqrt(2 * e.book.peakKinetic / w->ragdoll.mass()));
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Library lib;
    std::vector<std::string> errors;
    if (!loadDefaultLibrary(lib, &errors)) {
        for (const auto& e : errors) std::cerr << e << "\n";
        return 1;
    }
    const std::string cmd = argc > 1 ? argv[1] : "help";
    if (cmd == "run") return run(lib, argc, argv);
    if (cmd == "export") return exportScene(lib, argc, argv);
    if (cmd == "validate") return validate(lib);
    if (cmd == "list") {
        for (const auto& d : lib.liquids) std::printf("liquid      %-20s %s\n", d.id.c_str(), d.name.c_str());
        for (const auto& d : lib.surfaces) std::printf("surface     %-20s %s\n", d.id.c_str(), d.name.c_str());
        for (const auto& d : lib.vessels) std::printf("vessel      %-20s %s\n", d.id.c_str(), d.name.c_str());
        for (const auto& d : lib.projectiles) std::printf("projectile  %-20s %s\n", d.id.c_str(), d.name.c_str());
        for (const auto& d : lib.dummies) std::printf("dummy       %-20s %s\n", d.id.c_str(), d.name.c_str());
        for (const auto& d : lib.impactors) std::printf("impactor    %-20s %s\n", d.id.c_str(), d.name.c_str());
        return 0;
    }
    std::printf(
        "lqdemo run '<scene json>' [seconds]                       run a scene, print its log and the stain analysis\n"
        "lqdemo export '<scene json>' <s> <fps> <dir> [--ply] [--no-scene]   glTF + stains JSON (+ PLY per frame)\n"
        "lqdemo validate                                           the comparisons with published data\n"
        "lqdemo list                                               the library\n"
        "scene json keys: scenario (wine_glass, bottle_drop, gunshot, drip, angles, dam_break, impact, fall), tier 0-4,\n"
        "  seed, liquid, vessel, floor, table, wall, fill_ml, knock_ms, height_m, tilt_deg, projectile, dummy, target,\n"
        "  distance_m, shooter_height_m, speed_ms, pane (none/front/behind/both), glazing, pane_width_m, pane_height_m,\n"
        "  pane_thickness_mm, pane_gap_m, defer_panes, spatter_boost, knockback, fire, ragdoll, impactor, wall_m,\n"
        "  dummy_yaw_deg, fall_m\n");
    return 0;
}
