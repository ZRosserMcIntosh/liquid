// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/spray.hpp"

#include <algorithm>
#include <cmath>

#include "lqd/fluid.hpp"

namespace lqd {

double dragCoefficient(double Re) {
    if (Re <= 1e-9) return 1e9;
    if (Re < 1000) return 24.0 / Re * (1.0 + 0.15 * std::pow(Re, 0.687));
    return 0.44;
}

double terminalSpeed(double d, double rho) {
    double v = 1;
    for (int i = 0; i < 100; ++i) {
        const double Re = kAirDensity * v * d / kAirViscosity;
        const double vn = std::sqrt(4.0 * rho * kGravity * d / (3.0 * kAirDensity * dragCoefficient(Re)));
        if (std::abs(vn - v) < 1e-9) break;
        v = 0.5 * (v + vn);
    }
    return v;
}

StainShape stainShape(const LiquidDef& L, double d, double vn, double vt, const SurfaceDef* surface) {
    StainShape s;
    const double mu = L.impactViscosity();
    vn = std::max(vn, 1e-4);
    s.alpha = std::atan2(vn, std::max(vt, 0.0));
    const double Re = L.density * vn * d / mu;
    const double We = L.density * vn * vn * d / L.surfaceTension;
    const double P = We * std::pow(Re, -0.4);
    const double beta = std::pow(Re, 0.2) * std::sqrt(P) / (1.24 + std::sqrt(P));
    // A slow drop still wets a cap about its own size.
    const double absorb = surface ? surface->absorbency : 0.0;
    s.width = d * std::max(beta, 1.0) * (1.0 + 0.5 * absorb);
    const double sa = std::sin(std::max(s.alpha, 4.0 * kPi / 180));
    s.length = s.width / sa;
    s.K = std::sqrt(We) * std::pow(Re, 0.25);
    s.Kc = 57.7 * (surface ? surface->splashFactor : 1.0);
    s.splash = s.K > s.Kc;
    if (s.splash) {
        const double over = s.K / s.Kc - 1.0;
        s.spines = static_cast<int>(std::round(clamp(6 + 14 * over, 6, 48)));
        s.satelliteFraction = (1.0 - absorb) * clamp(0.04 + 0.08 * over, 0.0, 0.25);
        s.satellites = s.satelliteFraction > 0 ? static_cast<int>(std::round(clamp(2 + 6 * over, 2, 14))) : 0;
    }
    return s;
}

int Spray::addSource(const std::string& kind, const Vec3& origin, double time) {
    sources.push_back({kind, origin, time, 0, 0});
    return static_cast<int>(sources.size()) - 1;
}

void Spray::emit(const Droplet& d) {
    if (static_cast<int>(drops.size()) >= maxDrops) {
        lostVolume += kPi / 6 * d.d * d.d * d.d;
        return;
    }
    drops.push_back(d);
    if (d.source >= 0 && d.source < static_cast<int>(sources.size())) {
        sources[d.source].volume += kPi / 6 * d.d * d.d * d.d;
        ++sources[d.source].drops;
    }
}

double Spray::airVolume() const {
    double v = 0;
    for (const auto& d : drops) v += kPi / 6 * d.d * d.d * d.d;
    return v;
}

void Spray::land(const Droplet& d, const Collider& c, const Vec3& at, const Vec3& n, double t,
                 std::vector<Droplet>& born, Rng& rng) {
    const double vol = kPi / 6 * d.d * d.d * d.d;
    const Vec3 vrel = d.v;
    const double vn = -dot(vrel, n);
    const Vec3 vtan = vrel + n * vn;
    const double vt = length(vtan);
    const StainShape s = stainShape(*liquid_, d.d, vn, vt, c.surface);
    double kept = vol;
    if (s.satellites > 0 && static_cast<int>(drops.size() + born.size()) + s.satellites <= maxDrops) {
        // Satellite drops thrown from the crown, mostly onward along the tail.
        const double each = vol * s.satelliteFraction / s.satellites;
        const double ds = std::cbrt(6 * each / kPi);
        for (int k = 0; k < s.satellites; ++k) {
            Droplet sd;
            sd.d = ds;
            sd.t = t;
            sd.source = d.source;
            sd.satellite = true;
            const Vec3 u = vt > 1e-6 ? vtan / vt : perpendicular(n);
            const Vec3 w = cross(n, u);
            const double phi = rng.uniform(-1.0, 1.0) * (vt > 0.5 * vn ? 0.9 : kPi);
            const Vec3 dir = normalize(u * std::cos(phi) + w * std::sin(phi));
            const double sp = rng.uniform(0.3, 0.8) * std::sqrt(vn * vn + vt * vt);
            sd.v = dir * sp + n * (rng.uniform(0.15, 0.45) * sp);
            sd.x = at + n * (0.5 * ds + 2e-4) + dir * (0.5 * s.width);
            born.push_back(sd);
            kept -= each;
        }
        ++splashes;
    }
    if (!c.stains) {
        // Wetting a glass body: it runs off and joins the liquid.
        if (c.body >= 0 && particleVolume > 0) {
            toBulk += kept;
            toBulkAt = at + n * 0.003;
            toBulkV = {};
        } else {
            absorbedVolume += kept;
        }
        return;
    }
    Stain st;
    st.pos = at;
    st.normal = n;
    st.dir = vt > 1e-6 ? vtan / vt : perpendicular(n);
    st.width = s.width;
    st.length = s.length;
    st.alpha = s.alpha;
    st.diameter = d.d;
    st.speed = length(d.v);
    st.volume = kept;
    st.time = t;
    st.collider = c.id;
    st.source = d.source;
    st.spines = s.spines;
    st.splashed = s.splash;
    st.satellite = d.satellite;
    stains.push_back(st);
    stainedVolume += kept;
    if (c.body < 0 && c.role != Role::Dummy) pool(c, at, n, kept);
}

void Spray::pool(const Collider& c, const Vec3& at, const Vec3& n, double volume) {
    if (particleVolume <= 0) return;
    const double cell = 0.01;
    auto q = [&](double v) { return static_cast<std::uint64_t>(static_cast<std::int64_t>(std::floor(v / cell)) & 0xFFFFF); };
    const std::uint64_t key = (static_cast<std::uint64_t>(c.id) << 60) ^ (q(at.x) << 40) ^ (q(at.y) << 20) ^ q(at.z);
    Cell& k = cells_[key];
    k.volume += volume;
    const double film = cell * cell * filmThickness;
    // At half a particle beyond the film it pools; the cell lends the rest of
    // its film (the totals stay exact).
    while (k.volume - film >= 0.5 * particleVolume) {
        k.volume -= particleVolume;
        stainedVolume -= particleVolume;
        pooled.push_back({at, n});
    }
}

void Spray::step(double tEnd, const std::vector<Collider>& colliders, const Fluid* fluid, std::vector<Droplet>* rejoin,
                 Rng& rng) {
    if (!liquid_) return;
    const double rhoL = liquid_->density, sigma = liquid_->surfaceTension;
    const Vec3 g{0, -kGravity, 0};
    std::vector<Droplet> born, keep;
    keep.reserve(drops.size());
    const bool haveFluid = fluid && fluid->size() > 0;
    for (std::size_t idx = 0; idx < drops.size(); ++idx) {
        Droplet d = drops[idx];
        bool alive = true;
        while (alive && d.t < tEnd - 1e-12) {
            const double sp = length(d.v);
            double h = tEnd - d.t;
            h = std::min(h, 0.004 / std::max(sp, 1e-3));  // at most 4 mm per substep
            h = std::min(h, 0.002);
            const double Re = kAirDensity * sp * d.d / kAirViscosity;
            const double k = 0.75 * kAirDensity * dragCoefficient(Re) * sp / (rhoL * d.d);
            const Vec3 vNew = (d.v + g * h) / (1.0 + k * h);
            const Vec3 xNew = d.x + vNew * h;
            // First surface crossed on the way.
            double bestF = 2;
            int bestC = -1;
            Vec3 bestN;
            const double rad = 0.5 * d.d;
            for (const auto& c : colliders) {
                if (!c.enabled) continue;
                if (c.kind != ShapeKind::Plane) {
                    // Segment vs bounding sphere.
                    const Vec3 mid = (d.x + xNew) * 0.5;
                    if (length(mid - c.pos) > c.boundRadius + 0.5 * length(xNew - d.x) + rad) continue;
                }
                Vec3 n1;
                const double d1 = c.sdf(xNew, &n1);
                if (d1 > rad) continue;
                const double d0 = c.sdf(d.x, nullptr);
                if (d0 <= rad) {
                    // Started touching (born on the surface): land at once.
                    if (0 < bestF) {
                        bestF = 0;
                        bestC = c.id;
                        bestN = n1;
                    }
                    continue;
                }
                double lo = 0, hi = 1;
                for (int it = 0; it < 12; ++it) {
                    const double mid = 0.5 * (lo + hi);
                    (c.sdf(d.x + (xNew - d.x) * mid, nullptr) > rad ? lo : hi) = mid;
                }
                if (hi < bestF) {
                    bestF = hi;
                    bestC = c.id;
                    c.sdf(d.x + (xNew - d.x) * hi, &bestN);
                }
            }
            if (bestC >= 0) {
                const double th = h * bestF;
                const Vec3 at = d.x + (xNew - d.x) * bestF;
                d.v = d.v + (vNew - d.v) * bestF;
                d.t += th;
                const Collider& c = colliders[bestC];
                Vec3 surfN = bestN;
                const double dd = c.sdf(at, &surfN);
                const Vec3 onSurface = at - surfN * dd;
                // A blob of the bulk liquid joins the puddle; when it splashes,
                // it is thrown off as spatter instead with the probability of
                // the splash's satellite fraction, so on average that share of
                // the liquid leaves as satellites and stains (volume exact).
                bool splash = false;
                // Away from the bulk liquid (a wall, the far floor, the table)
                // a blob stays where it lands, as a stain.
                if (d.fromBulk && rejoin && !(haveFluid && fluid->densityAt(onSurface + surfN * (0.5 * d.d)) > 0.08 * fluid->restDensity()))
                    splash = true;
                if (d.fromBulk && rejoin && !splash) {
                    const double vn = -dot(d.v, surfN);
                    const StainShape sh = stainShape(*liquid_, d.d, std::max(vn, 0.0), length(d.v + surfN * vn), c.surface);
                    splash = sh.satellites > 0 && rng.uniform() < sh.satelliteFraction;
                }
                if (d.fromBulk && rejoin && !splash) {
                    d.x = onSurface + surfN * (0.5 * d.d);
                    rejoin->push_back(d);
                    bulkReturned += kPi / 6 * d.d * d.d * d.d;
                } else {
                    land(d, c, onSurface, surfN, d.t, born, rng);
                }
                alive = false;
                break;
            }
            d.x = xNew;
            d.v = vNew;
            d.t += h;
            // Into the bulk liquid?
            if (haveFluid && fluid->densityAt(d.x) > 0.5 * fluid->restDensity()) {
                if (d.fromBulk && rejoin) {
                    rejoin->push_back(d);
                    bulkReturned += kPi / 6 * d.d * d.d * d.d;
                } else {
                    toBulk += kPi / 6 * d.d * d.d * d.d;
                    toBulkAt = d.x;
                    toBulkV = d.v;
                }
                alive = false;
                break;
            }
            // Aerodynamic breakup.
            const double We = kAirDensity * sp * sp * d.d / sigma;
            if (We > 12) {
                d.breakup += h;
                const double tstar = d.d / std::max(sp, 1e-6) * std::sqrt(rhoL / kAirDensity);
                const double tb = 6.0 * tstar * std::pow(We - 12.0, -0.25);
                if (d.breakup >= tb && static_cast<int>(drops.size() + born.size()) + 32 <= maxDrops) {
                    const double dc = std::max(6.0 * sigma / (kAirDensity * sp * sp), 0.25 * d.d);
                    const int N = std::clamp(static_cast<int>(std::round(std::pow(d.d / dc, 3))), 2, 32);
                    const double dch = d.d / std::cbrt(static_cast<double>(N));
                    const Vec3 u = normalize(d.v);
                    for (int k2 = 0; k2 < N; ++k2) {
                        Droplet c = d;
                        c.d = dch;
                        c.breakup = 0;
                        c.fromBulk = false;
                        c.v = d.v + rng.cone(u, 1.2) * (0.12 * sp) - u * (0.12 * sp * 0.5);
                        c.x = d.x + rng.cone(u, kPi) * (0.5 * d.d);
                        born.push_back(c);
                    }
                    ++breakups;
                    alive = false;
                    break;
                }
            } else {
                d.breakup = 0;
            }
            if (d.x.y < -20) {
                lostVolume += kPi / 6 * d.d * d.d * d.d;
                alive = false;
            }
        }
        if (alive) keep.push_back(d);
    }
    drops.swap(keep);
    // Drops born this step (satellites, breakup) fly the rest of it too.
    if (!born.empty()) {
        std::vector<Droplet> mine;
        mine.swap(drops);
        const std::size_t room = static_cast<std::size_t>(std::max(0, maxDrops - static_cast<int>(mine.size())));
        for (std::size_t k = room; k < born.size(); ++k) lostVolume += kPi / 6 * born[k].d * born[k].d * born[k].d;
        if (born.size() > room) born.resize(room);
        drops = born;
        for (auto& b : drops)
            if (b.source >= 0 && b.source < static_cast<int>(sources.size())) ++sources[b.source].drops;
        step(tEnd, colliders, fluid, rejoin, rng);
        drops.insert(drops.begin(), mine.begin(), mine.end());
    }
}

}  // namespace lqd
