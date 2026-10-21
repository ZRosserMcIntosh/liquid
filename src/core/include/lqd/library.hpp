// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// The definitions a scene is built from: liquids, surfaces, vessels,
// projectiles and dummies. All are JSON (data/*.json); every value carries a
// provenance note ("published: ...", "calibrated: ...", "estimated: ...").
#pragma once

#include <map>
#include <string>
#include <vector>

#include "lqd/json.hpp"
#include "lqd/math.hpp"

namespace lqd {

using Provenance = std::map<std::string, std::string>;

enum class Rheology { Newtonian, Carreau };

struct LiquidDef {
    std::string id, name, category, description;
    double density = 1000;          // kg/m^3
    double surfaceTension = 0.072;  // N/m
    double viscosity = 1e-3;        // Pa s (Newtonian; Carreau: the zero-shear value is mu0)
    Rheology rheology = Rheology::Newtonian;
    double mu0 = 0, muInf = 0, lambda = 0, powerIndex = 1;  // Carreau-Yasuda with a = 2
    Vec3 color{0.5, 0.05, 0.1};     // linear RGB of the bulk liquid
    Vec3 stainColor{0.3, 0.02, 0.06};
    double opacity = 0.9;
    Provenance provenance;

    // Dynamic viscosity at shear rate `gamma` (1/s).
    double viscosityAt(double gamma) const;
    // The viscosity a drop impact sees (shear rates of 10^3-10^5 1/s): the
    // infinite-shear value of a shear-thinning liquid.
    double impactViscosity() const { return rheology == Rheology::Carreau ? muInf : viscosity; }
};

struct SurfaceDef {
    std::string id, name, description;
    double youngs = 50e9;     // Pa, of the surface layer
    double poisson = 0.25;
    double layer = 0;         // m, compliant layer thickness (carpet pile); 0 = half-space
    double restitution = 0.3; // for a glass object bouncing on it
    double friction = 0.4;    // Coulomb, glass on it
    double roughness = 1e-6;  // m, Ra
    double absorbency = 0;    // 0 sealed .. 1 wicks the stain (carpet, fabric)
    double hysteresis = 0.3;  // contact-angle hysteresis cos(receding) - cos(advancing): pins spreading liquid
    double splashFactor = 1;  // multiplies the splash threshold K (rough < 1)
    Vec3 color{0.7, 0.7, 0.7};
    Provenance provenance;
};

// A revolved vessel: the wall's mid-surface profile from the rim or lip down
// the inside, through the bottom (and stem) to the centre of the base. Each
// point is (r, y, thickness) in metres with y up from the base.
struct ProfilePoint {
    double r = 0, y = 0, t = 0.001;
};

struct VesselDef {
    std::string id, name, category, description;
    std::vector<ProfilePoint> profile;
    int cavityEnd = 0;          // profile index of the cavity's bottom centre
    double glassDensity = 2500; // kg/m^3
    double youngs = 70e9, poisson = 0.22;
    double strength = 120e6;    // Pa, Weibull scale of the wall's bending strength
    double weibullModulus = 6;
    double capacity = 0;        // m^3 brimful, nominal (checked against the profile)
    double fill = 0;            // m^3 default fill
    std::string liquid;         // default liquid id
    Provenance provenance;
};

struct ProjectileDef {
    std::string id, name, description;
    double mass = 0.008, diameter = 0.009, speed = 360;  // kg, m, m/s
    double cdTissue = 0.3;      // drag coefficient in soft-tissue simulant
    std::string glassImpactor;  // id of the same round (one pellet of a shot load) in the Glass Dynamics library
    int pellets = 1;            // a shot load fires this many identical pellets
    double spread = 0;          // rad, rms angle of the pellets about the aim
    double shotString = 0;      // m, how far the pellets string out along the path
    Provenance provenance;
};

struct Capsule {
    std::string name;
    Vec3 a, b;
    double radius = 0.05;
};

// How a wound throws liquid: volume per joule deposited, drop sizes,
// speeds and the cone it leaves in.
struct SpatterModel {
    double volumePerJoule = 0;  // m^3/J
    double medianDiameter = 3e-4, sigmaLn = 0.6, minDiameter = 5e-5, maxDiameter = 3e-3;  // m
    double speedFactor = 0.05;  // drop speed per m/s of the round (fastest drops)
    double maxSpeed = 30;       // m/s
    double coneHalfAngle = 0.6; // rad
};

struct DummyDef {
    std::string id, name, description;
    std::vector<Capsule> segments;  // standing, feet at y = 0, facing +z
    double tissueDensity = 1060;    // kg/m^3 of the simulant the rounds cross
    double mass = 0;                // kg, the whole dummy (0: from the capsules' volume)
    double tissueStrength = 4.9e6;  // Pa, its Poncelet resistance (the constant part of the retarding stress)
    std::string liquid;
    SpatterModel back, forward;
    double bleedRate = 0;           // m^3/s from each wound afterwards
    Provenance provenance;
};

class Library {
public:
    // Adds the definitions of one JSON document (any of the top-level arrays
    // "liquids", "surfaces", "vessels", "projectiles", "dummies"). Returns errors.
    std::vector<std::string> load(const std::string& json);
    // Checks the cross references; returns errors.
    std::vector<std::string> resolve() const;

    const LiquidDef* liquid(const std::string& id) const;
    const SurfaceDef* surface(const std::string& id) const;
    const VesselDef* vessel(const std::string& id) const;
    const ProjectileDef* projectile(const std::string& id) const;
    const DummyDef* dummy(const std::string& id) const;

    std::vector<LiquidDef> liquids;
    std::vector<SurfaceDef> surfaces;
    std::vector<VesselDef> vessels;
    std::vector<ProjectileDef> projectiles;
    std::vector<DummyDef> dummies;
};

// The library shipped in data/ (read from LQD_DATA_DIR when the build has it).
// Returns false and the errors when a file is missing or malformed.
bool loadDefaultLibrary(Library& lib, std::vector<std::string>* errors = nullptr, const std::string& dir = "");

}  // namespace lqd
