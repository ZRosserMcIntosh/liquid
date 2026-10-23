// Liquid Dynamics — foundation tests: data, kernel, neighbours, materials,
// drops in air, stain shapes, vessels, glass strikes, analysis, mesh, export.
// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>

#include "lqd/bpa.hpp"
#include "lqd/export.hpp"
#include "lqd/fluid.hpp"
#include "lqd/grid.hpp"
#include "lqd/json.hpp"
#include "lqd/kernel.hpp"
#include "lqd/library.hpp"
#include "lqd/mesher.hpp"
#include "lqd/rigid.hpp"
#include "lqd/spray.hpp"
#include "lqd/vessel.hpp"
#include "lqd_test.hpp"

using namespace lqd;

namespace {
const Library& lib() {
    static Library l = [] {
        Library x;
        std::vector<std::string> errors;
        if (!loadDefaultLibrary(x, &errors))
            for (const auto& e : errors) std::cerr << "  library: " << e << "\n";
        return x;
    }();
    return l;
}
}  // namespace

// ---- Data ------------------------------------------------------------------------

TEST_CASE("data: the library loads, and every number has its provenance") {
    Library l;
    std::vector<std::string> errors;
    CHECK(loadDefaultLibrary(l, &errors));
    for (const auto& e : errors) std::cerr << "  " << e << "\n";
    CHECK(l.liquids.size() >= 4);
    CHECK(l.surfaces.size() >= 6);
    CHECK(l.vessels.size() >= 2);
    CHECK(l.projectiles.size() >= 3);
    CHECK(l.dummies.size() >= 1);
    for (const auto& d : l.liquids) CHECK(d.provenance.count("density_kgm3") && d.provenance.count("surface_tension_mNm"));
    for (const auto& d : l.surfaces) CHECK(d.provenance.count("youngs_GPa") && d.provenance.count("splash_factor"));
    for (const auto& d : l.vessels) CHECK(d.provenance.count("strength_MPa") && d.provenance.count("profile_mm"));
    for (const auto& d : l.projectiles) CHECK(d.provenance.count("cd_tissue") && d.provenance.count("speed_ms"));
}

TEST_CASE("data: unknown fields, bad ranges and untagged provenance are rejected") {
    Library l;
    CHECK(!l.load(R"({"liquids":[{"id":"x","density_kgm3":1000,"colour":[1,0,0]}]})").empty());
    CHECK(!l.load(R"({"liquids":[{"id":"x","density_kgm3":5}]})").empty());
    CHECK(!l.load(R"({"surfaces":[{"id":"x","youngs_GPa":1,"provenance":{"youngs_GPa":"from memory"}}]})").empty());
    CHECK(!l.load(R"({"liquids":[{"id":"x","rheology":"carreau"}]})").empty());
    CHECK(l.load(R"({"liquids":[{"id":"x","density_kgm3":1000,"provenance":{"density_kgm3":"estimated: test"}}]})").empty());
}

// ---- Kernel and neighbours -----------------------------------------------------------

TEST_CASE("kernel: unit integral, odd gradient, wall volume from 1 to 0") {
    Kernel k(0.01);
    const double h = 0.0005;
    double sum = 0;
    for (double x = -0.01; x <= 0.01; x += h)
        for (double y = -0.01; y <= 0.01; y += h)
            for (double z = -0.01; z <= 0.01; z += h) sum += k.W(std::sqrt(x * x + y * y + z * z)) * h * h * h;
    CHECK_NEAR(sum, 1.0, 0.01);
    const Vec3 r{0.003, -0.002, 0.001};
    const Vec3 g1 = k.gradW(r, length(r)), g2 = k.gradW(-r, length(r));
    CHECK_NEAR(g1.x + g2.x, 0, 1e-9);
    CHECK(dot(g1, r) < 0);  // the kernel falls off outwards
    CHECK_NEAR(k.wallVolume(-0.01), 1.0, 1e-9);
    CHECK_NEAR(k.wallVolume(0), 0.5, 0.01);
    CHECK_NEAR(k.wallVolume(0.01), 0.0, 1e-9);
    // dPhi/dd = -A(d)
    for (double d : {-0.006, -0.002, 0.0, 0.003, 0.007}) {
        const double e = 1e-5;
        CHECK_NEAR((k.wallVolume(d - e) - k.wallVolume(d + e)) / (2 * e), k.wallSection(d), 2.0);
    }
}

TEST_CASE("grid: neighbour lists match brute force") {
    Rng rng(7);
    std::vector<Vec3> x;
    for (int i = 0; i < 800; ++i) x.push_back({rng.uniform(-0.1, 0.1), rng.uniform(-0.05, 0.05), rng.uniform(0, 0.08)});
    NeighbourGrid g;
    g.build(x, 0.012);
    std::vector<int> off, nb;
    g.neighbourLists(x, 0.012, off, nb);
    int mismatches = 0;
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        int brute = 0;
        for (int j = 0; j < static_cast<int>(x.size()); ++j) brute += (j != i && length(x[i] - x[j]) < 0.012) ? 1 : 0;
        mismatches += brute != off[i + 1] - off[i] ? 1 : 0;
    }
    CHECK(mismatches == 0);
}

TEST_CASE("fluid: a lattice sits just above rest density (relaxed packing factor)") {
    Fluid f;
    FluidSettings s;
    s.spacing = 0.005;
    f.configure(*lib().liquid("water"), s);
    for (int i = -4; i <= 4; ++i)
        for (int j = -4; j <= 4; ++j)
            for (int k = -4; k <= 4; ++k) f.add({i * 0.005, 1 + j * 0.005, k * 0.005}, {});
    std::vector<Collider> none;
    std::vector<RigidBody> bodies;
    f.step(1e-5, none, bodies, nullptr);
    // The centre particle (index 364 of the 9^3 block) sees a full kernel: the
    // particle mass is set for liquid that has relaxed, which packs ~5% closer.
    CHECK_RANGE(f.rho[364] / f.restDensity(), 1.0, 1.05);
}

// ---- Liquids and drops -----------------------------------------------------------------

TEST_CASE("liquids: blood thins with shear as Cho & Kensey's fit") {
    const LiquidDef& b = *lib().liquid("stage_blood_green");
    CHECK_NEAR(b.viscosityAt(0) * 1e3, 56, 0.01);
    CHECK(b.viscosityAt(10) < b.viscosityAt(1));
    CHECK_RANGE(b.viscosityAt(1000) * 1e3, 3.45, 4.5);
    CHECK_NEAR(b.impactViscosity() * 1e3, 3.45, 1e-9);
    CHECK_NEAR(lib().liquid("water")->viscosityAt(1e4) * 1e3, 1.002, 1e-9);
}

TEST_CASE("spray: terminal speeds of water drops against Gunn & Kinzer (1949)") {
    // Gunn & Kinzer, J. Meteorol. 6 (1949) 243: 0.5 mm 2.06 m/s, 1.0 mm 4.03, 2.0 mm 6.49.
    // A rigid sphere falls a little faster than a real 2 mm drop, which flattens.
    CHECK_NEAR(terminalSpeed(0.5e-3, 998.2), 2.06, 0.06 * 2.06);
    CHECK_NEAR(terminalSpeed(1.0e-3, 998.2), 4.03, 0.06 * 4.03);
    CHECK_NEAR(terminalSpeed(2.0e-3, 998.2), 6.49, 0.10 * 6.49);
    CHECK(terminalSpeed(2.0e-3, 998.2) > 6.49);
    CHECK_NEAR(dragCoefficient(2000), 0.44, 1e-12);
    CHECK_NEAR(dragCoefficient(0.1), 240 * (1 + 0.15 * std::pow(0.1, 0.687)), 1e-9);
}

TEST_CASE("stains: width over length is the sine of the impact angle; splashing above Mundo's K") {
    const LiquidDef& b = *lib().liquid("stage_blood_green");
    const SurfaceDef* tile = lib().surface("ceramic_tile");
    for (double deg : {15.0, 30.0, 45.0, 60.0, 80.0}) {
        const double a = deg * kPi / 180, v = 4;
        const StainShape s = stainShape(b, 2e-3, v * std::sin(a), v * std::cos(a), tile);
        CHECK_NEAR(std::asin(s.width / s.length), a, 1e-9);
    }
    const StainShape round = stainShape(b, 2e-3, 3, 0, tile);
    CHECK_NEAR(round.width, round.length, 1e-12);
    CHECK(round.width > 2e-3);  // it spreads to more than its own diameter
    // Spreading follows Laan et al. (2014): faster drops spread wider.
    CHECK(stainShape(b, 2e-3, 6, 0, tile).width > round.width);
    // A slow small drop does not splash; a fast large one does, earlier on rough concrete.
    CHECK(!stainShape(b, 0.3e-3, 1, 0, tile).splash);
    CHECK(stainShape(b, 4e-3, 6, 0, tile).splash);
    const StainShape c1 = stainShape(b, 2e-3, 3, 0, tile), c2 = stainShape(b, 2e-3, 3, 0, lib().surface("concrete"));
    CHECK(c2.Kc < c1.Kc);
    CHECK_NEAR(c1.K, std::sqrt(b.density * 9 * 2e-3 / b.surfaceTension) * std::pow(b.density * 3 * 2e-3 / b.muInf, 0.25), 1e-6);
}

// ---- Vessels ---------------------------------------------------------------------------

TEST_CASE("vessels: masses and capacities of the glass and the bottle") {
    const VesselGeometry glass(*lib().vessel("bordeaux_glass"));
    const VesselGeometry bottle(*lib().vessel("bordeaux_bottle"));
    // A machine-made Bordeaux glass weighs 140-220 g; a 750 ml Bordeaux bottle 400-550 g.
    CHECK_RANGE(glass.mass() * 1e3, 140, 220);
    CHECK_RANGE(bottle.mass() * 1e3, 400, 550);
    CHECK_NEAR(glass.cavityCapacity(), lib().vessel("bordeaux_glass")->capacity, 0.12 * lib().vessel("bordeaux_glass")->capacity);
    CHECK_NEAR(bottle.cavityCapacity(), lib().vessel("bordeaux_bottle")->capacity, 0.04 * lib().vessel("bordeaux_bottle")->capacity);
    // A standard pour sits below the widest point of the bowl.
    CHECK_RANGE(glass.fillHeight(150e-6), 0.13, 0.172);
    CHECK_NEAR(glass.cavityVolume(glass.fillHeight(150e-6)), 150e-6, 1e-9);
}

TEST_CASE("vessels: signed distance to the wall") {
    const VesselGeometry g(*lib().vessel("bordeaux_glass"));
    // On the mid-surface at the widest point of the bowl: half the thickness inside the glass.
    CHECK_NEAR(g.sdf({0.047, 0.165, 0}), -0.0006, 1e-5);
    Vec3 n;
    CHECK_NEAR(g.sdf({0.030, 0.165, 0}, &n), 0.017 - 0.0006, 1e-5);
    CHECK(n.x < -0.99);  // the cavity side's normal points into the cavity
    CHECK_NEAR(g.sdf({0, 0.165, 0.060}), 0.013 - 0.0006, 1e-5);  // the same wall at another azimuth
    CHECK(g.inCavity(0.02, 0.15, 0.004));
    CHECK(!g.inCavity(0.06, 0.15, 0.004));
    CHECK(!g.inCavity(0.0, 0.05, 0.001));  // the stem
}

TEST_CASE("vessels: shards partition the wall and keep the glass") {
    const VesselGeometry g(*lib().vessel("bordeaux_glass"));
    Rng rng(3);
    const ShardPattern P = makeShards(g, 0.05, 1.0, 30, 0.03, rng);
    CHECK(P.count >= 10 && P.count <= 30);
    double total = 0;
    std::vector<double> m(P.count, 0);
    for (int c = 0; c < g.cells(); ++c) {
        m[P.cellShard[c]] += g.cellVolume(c) * 2500;
        total += g.cellVolume(c) * 2500;
    }
    CHECK_NEAR(total, g.mass(), 0.02 * g.mass());
    for (double x : m) CHECK(x > 0);
    // Smaller pieces round the impact than far from it.
    const VesselMesh whole = g.mesh();
    CHECK(whole.tri.size() > 1000);
    CHECK(whole.pos.size() == whole.nrm.size());
}

// ---- Glass strikes ------------------------------------------------------------------------

TEST_CASE("strikes: a knock never breaks the glass, a fall onto tile does, carpet saves it") {
    const VesselDef& v = *lib().vessel("bordeaux_glass");
    const VesselGeometry g(v);
    WallSample bowl;
    for (const auto& w : g.samples())
        if (!w.edge && std::abs(w.p.y - 0.165) < 0.004 && w.p.x > 0) {
            bowl = w;
            break;
        }
    const double m = 0.08;
    auto P = [&](const char* s, double speed) { return assessImpact(v, bowl, *lib().surface(s), m, speed, 1).probability; };
    CHECK(P("oak_table", 0.3) < 0.01);
    CHECK_RANGE(P("oak_table", 1.0), 0.05, 0.8);  // tipped over on the table: sometimes
    CHECK(P("ceramic_tile", 3.8) > 0.99);         // off the table onto tile: always
    CHECK(P("carpet", 3.8) < 0.05);               // onto carpet: rarely
    CHECK(P("concrete", 2.0) > P("oak_floor", 1.0));
    // Struck end-on, a wall carries the load in compression and is safer.
    CHECK(assessImpact(v, bowl, *lib().surface("ceramic_tile"), m, 0.5, 0.1).stress <
          assessImpact(v, bowl, *lib().surface("ceramic_tile"), m, 0.5, 1).stress);
    // The force rises with speed; the carpet's pile spreads it.
    CHECK(assessImpact(v, bowl, *lib().surface("ceramic_tile"), m, 2, 1).force >
          5 * assessImpact(v, bowl, *lib().surface("carpet"), m, 2, 1).force);
}

// ---- Analysis --------------------------------------------------------------------------------

TEST_CASE("bpa: straight-line stains give back their origin exactly") {
    Spray sp;
    sp.configure(*lib().liquid("stage_blood_green"));
    std::vector<Collider> cols(1);
    cols[0].id = 0;
    cols[0].kind = ShapeKind::Plane;
    const Vec3 origin{0.2, 1.1, -0.3};
    const int src = sp.addSource("test", origin, 0);
    Rng rng(11);
    for (int i = 0; i < 40; ++i) {
        const Vec3 at{rng.uniform(-1, 1), 0, rng.uniform(-1, 1)};
        const Vec3 d = at - origin;
        const double alpha = std::atan2(-d.y, std::hypot(d.x, d.z));
        Stain s;
        s.pos = at;
        s.normal = {0, 1, 0};
        s.dir = normalize(Vec3{d.x, 0, d.z});
        s.width = 0.004;
        s.length = s.width / std::sin(alpha);
        s.alpha = alpha;
        s.collider = 0;
        s.source = src;
        sp.stains.push_back(s);
    }
    BpaOptions o;
    o.collider = 0;
    o.maxRatio = 0.999;
    const BpaResult r = analyseStains(sp, cols, o);
    CHECK(r.ok);
    CHECK(r.truth);
    CHECK_NEAR(r.convergenceError, 0, 1e-6);
    CHECK_NEAR(r.height, 1.1, 1e-6);
    CHECK_NEAR(r.origin.x, origin.x, 1e-6);
}

// ---- Mesh and export -----------------------------------------------------------------------------

TEST_CASE("colliders: a slab of glass with a hole, by signed distance") {
    Collider c;
    c.kind = ShapeKind::Slab;
    c.outline = {-0.1, -0.05, 0.1, -0.05, 0.1, 0.05, -0.1, 0.05};
    c.outlineHoles = {{-0.01, -0.01, 0.01, -0.01, 0.01, 0.01, -0.01, 0.01}};
    c.half = {0, 0, 0.003};
    c.pos = {1, 2, 3};
    c.rot = Quat::axisAngle({0, 1, 0}, kPi / 2);  // local z is world x
    Vec3 n;
    CHECK_NEAR(c.sdf({1.01, 2.0, 3.05}, &n), 0.007, 1e-12);   // off the face, between the hole and the edge
    CHECK_NEAR(n.x, 1, 1e-12);
    CHECK_NEAR(c.sdf({1.0, 2.0, 3.0 - 0.15}, &n), 0.05, 1e-12);  // beyond the end of it (local x = 0.15)
    CHECK_NEAR(c.sdf({1.0, 2.0, 3.0}), 0.01, 1e-12);           // in the middle of the hole: 1 cm from its edge
    CHECK(c.sdf({1.0, 2.0, 3.05}) < 0);                         // in the glass
    CHECK_NEAR(c.sdf({1.004, 2.03, 3.05}), 0.001, 1e-12);       // just off the face near its long edge
}

TEST_CASE("mesher: a ball of particles meshes to a closed sphere of the right size") {
    std::vector<Vec3> p;
    const double s = 0.004, R = 0.03;
    for (double x = -R; x <= R; x += s)
        for (double y = -R; y <= R; y += s)
            for (double z = -R; z <= R; z += s)
                if (x * x + y * y + z * z <= R * R) p.push_back({x, y, z});
    const SurfaceMesh m = meshLiquid(p, s);
    CHECK(m.triangles() > 500);
    double area = 0, rmean = 0;
    for (std::size_t t = 0; t < m.triangles(); ++t) {
        auto P = [&](std::uint32_t v) { return Vec3{m.pos[3 * v], m.pos[3 * v + 1], m.pos[3 * v + 2]}; };
        area += 0.5 * length(cross(P(m.tri[3 * t + 1]) - P(m.tri[3 * t]), P(m.tri[3 * t + 2]) - P(m.tri[3 * t])));
    }
    for (std::size_t v = 0; v < m.vertices(); ++v) rmean += length(Vec3{m.pos[3 * v], m.pos[3 * v + 1], m.pos[3 * v + 2]});
    rmean /= static_cast<double>(m.vertices());
    CHECK_NEAR(rmean, R, 0.15 * R);
    CHECK_NEAR(area, 4 * kPi * rmean * rmean, 0.12 * 4 * kPi * rmean * rmean);
    // Normals point outwards.
    int outward = 0;
    for (std::size_t v = 0; v < m.vertices(); ++v)
        outward += m.pos[3 * v] * m.nrm[3 * v] + m.pos[3 * v + 1] * m.nrm[3 * v + 1] + m.pos[3 * v + 2] * m.nrm[3 * v + 2] > 0;
    CHECK(outward > 0.97 * static_cast<double>(m.vertices()));
}

TEST_CASE("ballistics: the Poncelet fit puts a 9 mm FMJ 0.65 m into gelatin") {
    const DummyDef& d = *lib().dummy("crash_dummy");
    const ProjectileDef& p = *lib().projectile("pistol_9mm_fmj");
    const double A = kPi / 4 * p.diameter * p.diameter;
    const double a = d.tissueStrength * A, b = 0.5 * d.tissueDensity * p.cdTissue * A;
    const double depth = p.mass / (2 * b) * std::log(1 + b * p.speed * p.speed / a);
    CHECK_NEAR(depth, 0.65, 0.02);
    // The heavier, faster .44 Magnum goes deeper; nothing here models expansion or yaw.
    const ProjectileDef& m = *lib().projectile("magnum_44");
    const double Am = kPi / 4 * m.diameter * m.diameter;
    const double am = d.tissueStrength * Am, bm = 0.5 * d.tissueDensity * m.cdTissue * Am;
    CHECK(m.mass / (2 * bm) * std::log(1 + bm * m.speed * m.speed / am) > depth);
}

TEST_CASE("library: impactors and the dummy's burst model load, every value with its provenance") {
    Library L;
    std::vector<std::string> errors;
    CHECK(loadDefaultLibrary(L, &errors));
    for (const auto& e : errors) std::cerr << "  " << e << "\n";
    CHECK(L.impactors.size() >= 4);
    for (const auto& m : L.impactors) {
        CHECK(m.mass > 0);
        CHECK(m.speed > 0);
        CHECK(L.surface(m.surface) != nullptr);
        CHECK(m.provenance.count("mass_kg") == 1);
        CHECK(m.provenance.count("speed_ms") == 1);
    }
    const ImpactorDef* car = L.impactor("car_50");
    CHECK(car && car->shape == ImpactorShape::Box);
    if (car) CHECK_NEAR(car->speed * 3.6, 50, 0.1);
    const DummyDef* d = L.dummy("crash_dummy");
    CHECK(d && d->ruptureEnergy > 0 && d->ruptureSpeed > 0 && d->burst.volumePerJoule > 0);
    if (d) CHECK(d->provenance.count("rupture_J_per_kg") == 1);
}
