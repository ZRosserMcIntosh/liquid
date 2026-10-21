// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/scenarios.hpp"

#include <cmath>

#include "lqd/json.hpp"

namespace lqd {

bool parseSceneConfig(const std::string& text, SceneConfig& c, std::string* error) {
    const json::ParseResult pr = json::parse(text);
    auto fail = [&](const std::string& m) {
        if (error) *error = m;
        return false;
    };
    if (!pr.ok()) return fail("JSON " + pr.error);
    if (!pr.value.isObject()) return fail("the scene config must be a JSON object");
    for (const auto& [k, v] : pr.value.members()) {
        auto str = [&](std::string& o) {
            if (!v.isString()) return fail(k + " must be a string");
            o = v.asString();
            return true;
        };
        auto num = [&](double& o, double scale) {
            if (!v.isNumber()) return fail(k + " must be a number");
            o = v.asNumber() * scale;
            return true;
        };
        bool ok = true;
        double tmp = 0;
        if (k == "scenario") ok = str(c.scenario);
        else if (k == "tier") { ok = num(tmp, 1); c.tier = static_cast<int>(tmp); }
        else if (k == "seed") { ok = num(tmp, 1); c.seed = static_cast<std::uint64_t>(std::max(0.0, tmp)); }
        else if (k == "liquid") ok = str(c.liquid);
        else if (k == "vessel") ok = str(c.vessel);
        else if (k == "floor") ok = str(c.floor);
        else if (k == "table") ok = str(c.table);
        else if (k == "wall") ok = str(c.wall);
        else if (k == "fill_ml") ok = num(c.fill, 1e-6);
        else if (k == "knock_ms") ok = num(c.knock, 1);
        else if (k == "height_m") ok = num(c.height, 1);
        else if (k == "tilt_deg") ok = num(c.tilt, 1);
        else if (k == "projectile") ok = str(c.projectile);
        else if (k == "dummy") ok = str(c.dummy);
        else if (k == "target") ok = str(c.target);
        else if (k == "distance_m") ok = num(c.distance, 1);
        else if (k == "shooter_height_m") ok = num(c.shooterHeight, 1);
        else if (k == "speed_ms") ok = num(c.speed, 1);
        else if (k == "pane") ok = str(c.pane);
        else if (k == "glazing") ok = str(c.glazing);
        else if (k == "pane_width_m") ok = num(c.paneWidth, 1);
        else if (k == "pane_height_m") ok = num(c.paneHeight, 1);
        else if (k == "pane_thickness_mm") ok = num(c.paneThickness, 1e-3);
        else if (k == "pane_gap_m") ok = num(c.paneGap, 1);
        else if (k == "ragdoll") {
            if (!v.isBool()) return fail("ragdoll must be true or false");
            c.ragdoll = v.asBool();
        } else if (k == "defer_panes") {
            if (!v.isBool()) return fail("defer_panes must be true or false");
            c.deferPanes = v.asBool();
        } else return fail("unknown scene key \"" + k + "\"");
        if (!ok) return false;
    }
    return true;
}

namespace {

Vec3 targetPoint(const std::string& t) {
    if (t == "abdomen") return {0, 1.05, 0};
    if (t == "shoulder") return {0.15, 1.43, 0};
    if (t == "thigh") return {0.10, 0.72, 0};
    return {0, 1.30, 0};  // chest
}

void room(World& w, const SurfaceDef* floor, const SurfaceDef* wall, double backZ) {
    w.addPlane(Role::Floor, "floor", {0, 0, 0}, {0, 1, 0}, floor);
    w.addPlane(Role::Wall, "back wall", {0, 0, backZ}, {0, 0, 1}, wall);
}

}  // namespace

ShotLine shotLine(const SceneConfig& c) {
    ShotLine s;
    s.aim = targetPoint(c.target);
    s.muzzle = {s.aim.x * 0.3, c.shooterHeight, c.distance};
    return s;
}

bool buildScene(World& w, const SceneConfig& c, std::string* error) {
    const Library& lib = w.library();
    auto fail = [&](const std::string& m) {
        if (error) *error = m;
        return false;
    };
    const SurfaceDef* floor = lib.surface(c.floor);
    const SurfaceDef* wall = lib.surface(c.wall);
    const SurfaceDef* table = lib.surface(c.table);
    if (!floor) return fail("unknown floor surface " + c.floor);
    if (!wall) return fail("unknown wall surface " + c.wall);

    if (c.scenario == "wine_glass" || c.scenario == "bottle_drop") {
        const bool glass = c.scenario == "wine_glass";
        const VesselDef* v = lib.vessel(!c.vessel.empty() ? c.vessel : glass ? "bordeaux_glass" : "bordeaux_bottle");
        if (!v) return fail("unknown vessel " + c.vessel);
        const LiquidDef* l = lib.liquid(!c.liquid.empty() ? c.liquid : v->liquid);
        if (!l) return fail("unknown liquid " + c.liquid);
        room(w, floor, wall, -1.2);
        if (glass) {
            if (!table) return fail("unknown table surface " + c.table);
            const double top = 0.75;
            w.addBox(Role::Table, "table", {0, top - 0.02, 0}, {0.6, 0.02, 0.4}, Quat{}, table);
            const double footR = v->profile.back().r + 0.5 * v->profile.back().t;
            const int b = w.addVessel(*v, {0.6 - footR - 0.025, top + 1e-4, 0.05}, Quat{});
            w.fillVessel(b, *l, c.fill > 0 ? c.fill : v->fill);
            // The knock: a sideways blow at the bowl's widest point.
            RigidBody& B = w.bodies[b];
            double hBowl = 0, rMax = 0;
            for (int k = 0; k <= v->cavityEnd; ++k)
                if (v->profile[k].r > rMax) {
                    rMax = v->profile[k].r;
                    hBowl = v->profile[k].y;
                }
            const Vec3 at = B.fromVessel({-rMax, hBowl, 0});
            B.applyImpulse({B.mass * c.knock, 0, 0}, at);
            w.log("knocked towards the table edge at the bowl, " + std::to_string(c.knock).substr(0, 4) + " m/s");
        } else {
            const Quat q = Quat::axisAngle({0, 0, 1}, c.tilt * kPi / 180);
            const int b = w.addVessel(*v, {0, c.height, 0}, q);
            w.fillVessel(b, *l, c.fill > 0 ? c.fill : v->fill);
            w.log("let go " + std::to_string(c.height).substr(0, 4) + " m above the floor");
        }
        return true;
    }

    if (c.scenario == "gunshot") {
        const DummyDef* d = lib.dummy(c.dummy);
        const ProjectileDef* p = lib.projectile(c.projectile);
        if (!d) return fail("unknown dummy " + c.dummy);
        if (!p) return fail("unknown projectile " + c.projectile);
        const LiquidDef* l = lib.liquid(!c.liquid.empty() ? c.liquid : d->liquid);
        if (!l) return fail("unknown liquid " + c.liquid);
        room(w, floor, wall, -2.5);
        w.setLiquid(*l);
        w.addDummy(*d, {0, 0, 0}, 0, c.ragdoll);
        const ShotLine s = shotLine(c);
        if (c.pane != "none" && c.pane != "front" && c.pane != "behind" && c.pane != "both")
            return fail("pane must be none, front, behind or both");
        for (int side = 0; side < 2; ++side) {
            const bool front = side == 0;
            if (c.pane == "none" || (front && c.pane == "behind") || (!front && c.pane == "front")) continue;
            const double gap = c.paneGap > 0 ? std::min(c.paneGap, front ? 0.9 * c.distance : 2.0) : -1;
            const double z = front ? (gap > 0 ? gap : std::min(1.0, 0.6 * c.distance)) : (gap > 0 ? -gap : -0.8);
            // The pane straddles the path at its own depth.
            const Vec3 dir = normalize(s.aim - s.muzzle);
            const double t = (z - s.muzzle.z) / dir.z;
            const Vec3 at = s.muzzle + dir * t;
            const double cy = std::max(at.y, 0.5 * c.paneHeight + 0.05);
            const int id = w.addBox(Role::Pane, c.pane == "both" ? (front ? "front pane" : "back pane") : "pane",
                                    {at.x, cy, z}, {0.5 * c.paneWidth, 0.5 * c.paneHeight, 0.5 * c.paneThickness},
                                    Quat{}, lib.surface("glass_pane"));
            w.colliders[id].glazing = c.glazing;
        }
        w.deferPanes = c.deferPanes;
        w.shoot(*p, s.muzzle, s.aim - s.muzzle, c.speed > 0 ? c.speed : p->speed, 0.01);
        return true;
    }

    if (c.scenario == "drip") {
        const LiquidDef* l = lib.liquid(!c.liquid.empty() ? c.liquid : "stage_blood_green");
        if (!l) return fail("unknown liquid " + c.liquid);
        room(w, floor, wall, -1.5);
        w.setLiquid(*l);
        // A 50 microlitre drop (the drop size a slow wound drips).
        const double d = std::cbrt(6 * 50e-9 / kPi);
        const int n = 7;
        for (int i = 0; i < n; ++i) {
            const double h = 0.3 + (c.height - 0.3) * i / (n - 1);
            w.drop({-0.6 + 1.2 * i / (n - 1), h, 0}, {0, 0, 0}, d, "drip");
        }
        return true;
    }

    if (c.scenario == "angles") {
        const LiquidDef* l = lib.liquid(!c.liquid.empty() ? c.liquid : "stage_blood_green");
        if (!l) return fail("unknown liquid " + c.liquid);
        w.addPlane(Role::Floor, "floor", {0, 0, 0}, {0, 1, 0}, floor);
        w.addPlane(Role::Wall, "card", {0, 0, 0}, {0, 0, 1}, wall);
        w.setLiquid(*l);
        const double d = 0.0035;
        const double sp = 4.0;
        for (int i = 0; i < 9; ++i) {
            const double a = (10 + 10 * i) * kPi / 180;
            const Vec3 v{0, -sp * std::cos(a), -sp * std::sin(a)};  // down the card at angle a
            const Vec3 hit{-0.8 + 0.2 * i, 1.0, 0};
            w.drop(hit - v * 0.02 + Vec3{0, 0, d}, v, d, "card");
        }
        return true;
    }

    if (c.scenario == "dam_break") {
        const LiquidDef* l = lib.liquid(!c.liquid.empty() ? c.liquid : "water");
        if (!l) return fail("unknown liquid " + c.liquid);
        const SurfaceDef* s = lib.surface("glass_pane");
        // The 0.146 m x 0.292 m column of Koshizuka & Oka (1996), after Martin & Moyce (1952),
        // a few particles deep between glass walls (quasi two-dimensional).
        const double a = 0.146, W = 0.584, depth = 4 * w.budget().spacing;
        w.addPlane(Role::Floor, "floor", {0, 0, 0}, {0, 1, 0}, s);
        w.addPlane(Role::Wall, "left", {0, 0, 0}, {1, 0, 0}, s);
        w.addPlane(Role::Wall, "right", {W, 0, 0}, {-1, 0, 0}, s);
        w.addPlane(Role::Wall, "front", {0, 0, depth}, {0, 0, -1}, s);
        w.addPlane(Role::Wall, "back", {0, 0, 0}, {0, 0, 1}, s);
        w.addLiquidBox(*l, {0, 0, 0}, {a, 2 * a, depth}, {0, 0, 0});
        return true;
    }

    return fail("unknown scenario " + c.scenario);
}

}  // namespace lqd
