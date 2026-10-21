// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/library.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

#include "lqd/json.hpp"

namespace lqd {

// Carreau-Yasuda with a = 2: mu = muInf + (mu0 - muInf) (1 + (lambda gamma)^2)^((n-1)/2).
double LiquidDef::viscosityAt(double gamma) const {
    if (rheology == Rheology::Newtonian) return viscosity;
    return muInf + (mu0 - muInf) * std::pow(1.0 + sq(lambda * gamma), 0.5 * (powerIndex - 1.0));
}

// ---- Lookup ------------------------------------------------------------------

namespace {
template <class T>
const T* findId(const std::vector<T>& v, const std::string& id) {
    for (const auto& x : v)
        if (x.id == id) return &x;
    return nullptr;
}
}  // namespace

const LiquidDef* Library::liquid(const std::string& id) const { return findId(liquids, id); }
const SurfaceDef* Library::surface(const std::string& id) const { return findId(surfaces, id); }
const VesselDef* Library::vessel(const std::string& id) const { return findId(vessels, id); }
const ProjectileDef* Library::projectile(const std::string& id) const { return findId(projectiles, id); }
const DummyDef* Library::dummy(const std::string& id) const { return findId(dummies, id); }

// ---- JSON loading --------------------------------------------------------------

namespace {

// Reads the members of one definition object; every key must be known.
class Reader {
public:
    Reader(const json::Value& obj, std::string where, std::vector<std::string>& errors)
        : obj_(obj), where_(std::move(where)), errors_(errors) {}

    void num(const char* key, double scale, double& out, double lo, double hi) {
        const json::Value* v = take(key);
        if (!v) return;
        if (!v->isNumber()) return error(std::string(key) + " must be a number");
        const double x = v->asNumber();
        if (!(x >= lo && x <= hi))
            return error(std::string(key) + " = " + std::to_string(x) + " is out of range [" + std::to_string(lo) +
                         ", " + std::to_string(hi) + "]");
        out = x * scale;
    }
    void str(const char* key, std::string& out) {
        const json::Value* v = take(key);
        if (!v) return;
        if (!v->isString()) return error(std::string(key) + " must be a string");
        out = v->asString();
    }
    void color(const char* key, Vec3& out) {
        const json::Value* v = take(key);
        if (!v) return;
        if (!v->isArray() || v->items().size() != 3) return error(std::string(key) + " must be [r, g, b]");
        for (int i = 0; i < 3; ++i) {
            const auto& c = v->items()[i];
            if (!c.isNumber() || c.asNumber() < 0 || c.asNumber() > 1)
                return error(std::string(key) + " components must be numbers in [0, 1]");
            out[i] = c.asNumber();
        }
    }
    void provenance(Provenance& out) {
        const json::Value* v = take("provenance");
        if (!v) return;
        if (!v->isObject()) return error("provenance must be an object");
        for (const auto& [k, s] : v->members()) {
            if (!s.isString()) {
                error("provenance." + k + " must be a string");
                continue;
            }
            const std::string& t = s.asString();
            if (t.rfind("published", 0) != 0 && t.rfind("calibrated", 0) != 0 && t.rfind("estimated", 0) != 0 &&
                t.rfind("gameplay", 0) != 0)
                error("provenance." + k + " must start with published, calibrated, estimated or gameplay");
            bool known = false;
            for (const auto& [mk, mv] : obj_.members()) known = known || mk == k;
            if (!known && k.find('.') == std::string::npos) error("provenance." + k + " names no field");
            out[k] = t;
        }
    }
    const json::Value* take(const char* key) {
        used_.push_back(key);
        return obj_.find(key);
    }
    void finish() {
        for (const auto& [k, v] : obj_.members()) {
            bool known = false;
            for (const auto& u : used_) known = known || u == k;
            if (!known) error("unknown field \"" + k + "\"");
        }
    }
    void error(const std::string& m) {
        errors_.push_back(where_ + ": " + m);
        ++count_;
    }
    int errorCount() const { return count_; }

private:
    const json::Value& obj_;
    std::string where_;
    std::vector<std::string>& errors_;
    std::vector<std::string> used_;
    int count_ = 0;
};

template <class T>
void upsert(std::vector<T>& list, T item) {
    for (auto& x : list)
        if (x.id == item.id) {
            x = std::move(item);
            return;
        }
    list.push_back(std::move(item));
}

std::string idOf(const json::Value& o) {
    const json::Value* v = o.find("id");
    return v && v->isString() ? v->asString() : std::string("?");
}

void readLiquid(const json::Value& o, std::vector<LiquidDef>& out, std::vector<std::string>& errors) {
    LiquidDef d;
    Reader r(o, "liquid " + idOf(o), errors);
    r.str("id", d.id);
    r.str("name", d.name);
    r.str("category", d.category);
    r.str("description", d.description);
    r.num("density_kgm3", 1, d.density, 500, 3000);
    r.num("surface_tension_mNm", 1e-3, d.surfaceTension, 1, 1000);
    r.num("viscosity_mPas", 1e-3, d.viscosity, 0.1, 1e6);
    std::string rheo = "newtonian";
    r.str("rheology", rheo);
    if (rheo == "carreau") d.rheology = Rheology::Carreau;
    else if (rheo != "newtonian") r.error("rheology must be newtonian or carreau");
    r.num("mu0_mPas", 1e-3, d.mu0, 0.1, 1e6);
    r.num("mu_inf_mPas", 1e-3, d.muInf, 0.1, 1e6);
    r.num("lambda_s", 1, d.lambda, 0, 1000);
    r.num("power_index", 1, d.powerIndex, 0.05, 1);
    r.color("color", d.color);
    r.color("stain_color", d.stainColor);
    r.num("opacity", 1, d.opacity, 0, 1);
    r.provenance(d.provenance);
    r.finish();
    if (d.id.empty()) r.error("missing id");
    if (d.rheology == Rheology::Carreau) {
        if (d.mu0 <= 0 || d.muInf <= 0 || d.lambda <= 0) r.error("carreau needs mu0_mPas, mu_inf_mPas and lambda_s");
        d.viscosity = d.mu0;
    }
    if (r.errorCount() == 0) upsert(out, std::move(d));
}

void readSurface(const json::Value& o, std::vector<SurfaceDef>& out, std::vector<std::string>& errors) {
    SurfaceDef d;
    Reader r(o, "surface " + idOf(o), errors);
    r.str("id", d.id);
    r.str("name", d.name);
    r.str("description", d.description);
    r.num("youngs_GPa", 1e9, d.youngs, 1e-6, 1000);
    r.num("poisson", 1, d.poisson, 0, 0.499);
    r.num("layer_mm", 1e-3, d.layer, 0, 100);
    r.num("restitution", 1, d.restitution, 0, 1);
    r.num("friction", 1, d.friction, 0, 2);
    r.num("roughness_um", 1e-6, d.roughness, 0, 5000);
    r.num("absorbency", 1, d.absorbency, 0, 1);
    r.num("hysteresis", 1, d.hysteresis, 0, 2);
    r.num("splash_factor", 1, d.splashFactor, 0.1, 3);
    r.color("color", d.color);
    r.provenance(d.provenance);
    r.finish();
    if (d.id.empty()) r.error("missing id");
    if (r.errorCount() == 0) upsert(out, std::move(d));
}

void readVessel(const json::Value& o, std::vector<VesselDef>& out, std::vector<std::string>& errors) {
    VesselDef d;
    Reader r(o, "vessel " + idOf(o), errors);
    r.str("id", d.id);
    r.str("name", d.name);
    r.str("category", d.category);
    r.str("description", d.description);
    if (const json::Value* p = r.take("profile_mm")) {
        if (!p->isArray() || p->items().size() < 3) r.error("profile_mm must be an array of [r, y, t] points");
        else
            for (const auto& q : p->items()) {
                if (!q.isArray() || q.items().size() != 3 || !q.items()[0].isNumber() || !q.items()[1].isNumber() ||
                    !q.items()[2].isNumber()) {
                    r.error("profile_mm points must be [r, y, t] numbers");
                    break;
                }
                ProfilePoint pp{q.items()[0].asNumber() * 1e-3, q.items()[1].asNumber() * 1e-3,
                                q.items()[2].asNumber() * 1e-3};
                if (pp.r < 0 || pp.t <= 0) r.error("profile_mm: r must be >= 0 and t > 0");
                d.profile.push_back(pp);
            }
    }
    double ce = 0;
    r.num("cavity_end", 1, ce, 1, 1000);
    d.cavityEnd = static_cast<int>(ce);
    r.num("glass_density_kgm3", 1, d.glassDensity, 1000, 6000);
    r.num("youngs_GPa", 1e9, d.youngs, 1, 500);
    r.num("poisson", 1, d.poisson, 0, 0.499);
    r.num("strength_MPa", 1e6, d.strength, 1, 5000);
    r.num("weibull_modulus", 1, d.weibullModulus, 1, 60);
    r.num("capacity_ml", 1e-6, d.capacity, 1, 1e5);
    r.num("fill_ml", 1e-6, d.fill, 0, 1e5);
    r.str("liquid", d.liquid);
    r.provenance(d.provenance);
    r.finish();
    if (d.id.empty()) r.error("missing id");
    if (d.cavityEnd <= 0 || d.cavityEnd >= static_cast<int>(d.profile.size()))
        r.error("cavity_end must index a profile point after the rim");
    if (r.errorCount() == 0) upsert(out, std::move(d));
}

void readProjectile(const json::Value& o, std::vector<ProjectileDef>& out, std::vector<std::string>& errors) {
    ProjectileDef d;
    Reader r(o, "projectile " + idOf(o), errors);
    r.str("id", d.id);
    r.str("name", d.name);
    r.str("description", d.description);
    r.num("mass_g", 1e-3, d.mass, 0.01, 10000);
    r.num("diameter_mm", 1e-3, d.diameter, 0.1, 200);
    r.num("speed_ms", 1, d.speed, 1, 2000);
    r.num("cd_tissue", 1, d.cdTissue, 0.01, 5);
    r.str("glass_impactor", d.glassImpactor);
    r.provenance(d.provenance);
    r.finish();
    if (d.id.empty()) r.error("missing id");
    if (r.errorCount() == 0) upsert(out, std::move(d));
}

void readSpatter(Reader& parent, const json::Value& o, SpatterModel& s, const std::string& where,
                 std::vector<std::string>& errors) {
    if (!o.isObject()) return parent.error(where + " must be an object");
    Reader r(o, where, errors);
    r.num("volume_per_joule_ul", 1e-9, s.volumePerJoule, 0, 1e4);
    r.num("median_diameter_mm", 1e-3, s.medianDiameter, 0.01, 10);
    r.num("sigma_ln", 1, s.sigmaLn, 0, 3);
    r.num("min_diameter_mm", 1e-3, s.minDiameter, 0.005, 10);
    r.num("max_diameter_mm", 1e-3, s.maxDiameter, 0.01, 20);
    r.num("speed_factor", 1, s.speedFactor, 0, 1);
    r.num("max_speed_ms", 1, s.maxSpeed, 0, 500);
    double cone = 30;
    r.num("cone_half_angle_deg", 1, cone, 1, 89);
    s.coneHalfAngle = cone * kPi / 180;
    r.finish();
}

void readDummy(const json::Value& o, std::vector<DummyDef>& out, std::vector<std::string>& errors) {
    DummyDef d;
    Reader r(o, "dummy " + idOf(o), errors);
    r.str("id", d.id);
    r.str("name", d.name);
    r.str("description", d.description);
    if (const json::Value* segs = r.take("segments_m")) {
        if (!segs->isArray()) r.error("segments_m must be an array");
        else
            for (const auto& s : segs->items()) {
                const json::Value *n = s.find("name"), *a = s.find("a"), *b = s.find("b"), *rad = s.find("radius");
                auto vec = [](const json::Value* v, Vec3& o2) {
                    if (!v || !v->isArray() || v->items().size() != 3) return false;
                    for (int i = 0; i < 3; ++i) {
                        if (!v->items()[i].isNumber()) return false;
                        o2[i] = v->items()[i].asNumber();
                    }
                    return true;
                };
                Capsule c;
                if (!n || !n->isString() || !vec(a, c.a) || !vec(b, c.b) || !rad || !rad->isNumber()) {
                    r.error("segments_m entries need name, a [x,y,z], b [x,y,z] and radius");
                    break;
                }
                c.name = n->asString();
                c.radius = rad->asNumber();
                d.segments.push_back(c);
            }
    }
    r.num("tissue_density_kgm3", 1, d.tissueDensity, 500, 3000);
    r.num("tissue_strength_MPa", 1e6, d.tissueStrength, 0, 100);
    r.str("liquid", d.liquid);
    if (const json::Value* b = r.take("back_spatter")) readSpatter(r, *b, d.back, "dummy " + d.id + " back_spatter", errors);
    if (const json::Value* f = r.take("forward_spatter"))
        readSpatter(r, *f, d.forward, "dummy " + d.id + " forward_spatter", errors);
    r.num("bleed_ml_per_s", 1e-6, d.bleedRate, 0, 1000);
    r.provenance(d.provenance);
    r.finish();
    if (d.id.empty()) r.error("missing id");
    if (d.segments.empty()) r.error("a dummy needs segments_m");
    if (r.errorCount() == 0) upsert(out, std::move(d));
}

}  // namespace

std::vector<std::string> Library::load(const std::string& text) {
    std::vector<std::string> errors;
    const json::ParseResult pr = json::parse(text);
    if (!pr.ok()) {
        errors.push_back("JSON " + pr.error);
        return errors;
    }
    if (!pr.value.isObject()) {
        errors.push_back("the document must be a JSON object");
        return errors;
    }
    for (const auto& [k, v] : pr.value.members()) {
        if (k == "comment") continue;
        if (!v.isArray()) {
            errors.push_back("\"" + k + "\" must be an array");
            continue;
        }
        for (const auto& item : v.items()) {
            if (!item.isObject()) {
                errors.push_back("\"" + k + "\" entries must be objects");
                continue;
            }
            if (k == "liquids") readLiquid(item, liquids, errors);
            else if (k == "surfaces") readSurface(item, surfaces, errors);
            else if (k == "vessels") readVessel(item, vessels, errors);
            else if (k == "projectiles") readProjectile(item, projectiles, errors);
            else if (k == "dummies") readDummy(item, dummies, errors);
            else {
                errors.push_back("unknown section \"" + k + "\"");
                break;
            }
        }
    }
    return errors;
}

std::vector<std::string> Library::resolve() const {
    std::vector<std::string> errors;
    for (const auto& v : vessels)
        if (!v.liquid.empty() && !liquid(v.liquid)) errors.push_back("vessel " + v.id + ": unknown liquid " + v.liquid);
    for (const auto& d : dummies)
        if (!d.liquid.empty() && !liquid(d.liquid)) errors.push_back("dummy " + d.id + ": unknown liquid " + d.liquid);
    return errors;
}

bool loadDefaultLibrary(Library& lib, std::vector<std::string>* errors, const std::string& dirIn) {
    std::string dir = dirIn;
#ifdef LQD_DATA_DIR
    if (dir.empty()) dir = LQD_DATA_DIR;
#endif
    std::vector<std::string> errs;
    for (const char* f : {"/liquids.json", "/surfaces.json", "/vessels.json", "/projectiles.json", "/dummies.json"}) {
        std::ifstream in(dir + f, std::ios::binary);
        if (!in) {
            errs.push_back("cannot read " + dir + f);
            continue;
        }
        std::stringstream ss;
        ss << in.rdbuf();
        for (auto& e : lib.load(ss.str())) errs.push_back(std::string(f + 1) + ": " + e);
    }
    for (auto& e : lib.resolve()) errs.push_back(e);
    if (errors) *errors = errs;
    return errs.empty();
}

}  // namespace lqd
