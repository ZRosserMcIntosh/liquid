// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/fluid.hpp"

#include "lqd/rng.hpp"

#include <algorithm>
#include <cmath>

namespace lqd {

namespace {
constexpr int kMaxTouch = 3;
constexpr double kRelaxedPacking = 1.05;

// prod over m != k of (1 - phi_m).
double others(const double* phi, int n, int k) {
    double p = 1;
    for (int m = 0; m < n; ++m)
        if (m != k) p *= 1 - phi[m];
    return p;
}

}  // namespace

void Fluid::configure(const LiquidDef& liquid, const FluidSettings& s) {
    liquid_ = &liquid;
    set_ = s;
    rho0_ = liquid.density;
    H_ = 2.0 * s.spacing;
    kernel_.setSupport(H_);
    // Particle mass so that liquid at the rest spacing sums to rho0. A liquid
    // that has relaxed under the density constraint does not keep the lattice
    // it was laid out in: it packs about 5% more particles into the same
    // volume at the same kernel sum. The mass carries that factor
    // (calibrated on a resting column, tests/test_scenes.cpp), so liquid at
    // rest fills the volume it was poured as.
    double sum = 0;
    const int n = 3;
    for (int i = -n; i <= n; ++i)
        for (int j = -n; j <= n; ++j)
            for (int k = -n; k <= n; ++k) sum += kernel_.W(s.spacing * std::sqrt(double(i * i + j * j + k * k)));
    mass_ = kRelaxedPacking * rho0_ / sum;
    selfDensity_ = kernel_.W(0) * mass_ / rho0_;
    // The wall's kernel volume is shifted so that a lattice's first layer rests
    // half a spacing from it, as it would against a mirrored lattice: the
    // particles of a half-space lattice are short of rest density by exactly
    // what the wall must supply there.
    double half = 0;
    for (int i = -n; i <= n; ++i)
        for (int j = 0; j <= 2 * n; ++j)
            for (int k = -n; k <= n; ++k)
                half += kernel_.W(s.spacing * std::sqrt(double(i * i + j * j + k * k)));
    const double need = std::max(0.0, 1.0 - half * mass_ / rho0_);
    double lo = -H_, hi = H_;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        (kernel_.wallVolume(mid) > need ? lo : hi) = mid;
    }
    wallShift_ = 0.5 * s.spacing - 0.5 * (lo + hi);
}

int Fluid::add(const Vec3& p, const Vec3& vel) {
    x.push_back(p);
    v.push_back(vel);
    rho.push_back(rho0_);
    return static_cast<int>(x.size()) - 1;
}

void Fluid::remove(const std::vector<int>& idx) {
    if (idx.empty()) return;
    std::size_t w = 0, k = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (k < idx.size() && static_cast<int>(i) == idx[k]) {
            ++k;
            continue;
        }
        x[w] = x[i];
        v[w] = v[i];
        rho[w] = rho[i];
        ++w;
    }
    x.resize(w);
    v.resize(w);
    rho.resize(w);
}

double Fluid::maxSpeed() const {
    double m = 0;
    for (const auto& u : v) m = std::max(m, length2(u));
    return std::sqrt(m);
}

double Fluid::maxStep() const {
    // CFL on the particle spacing, and the explicit viscous limit. The speed
    // that matters is relative: liquid riding in a falling bottle moves with
    // its walls and its neighbours, so it is measured against the liquid's
    // mean, except where a particle touches something that is not moving.
    Vec3 mean;
    for (const auto& u : v) mean += u;
    if (!v.empty()) mean = mean / static_cast<double>(v.size());
    double v2 = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const bool fixedNear = i < staticTouch_.size() && staticTouch_[i];
        v2 = std::max(v2, length2(fixedNear ? v[i] : v[i] - mean));
    }
    const double vmax = std::sqrt(v2);
    double dt = 0.4 * set_.spacing / std::max(vmax, 0.3);
    const double nu = liquid_ ? liquid_->viscosity / rho0_ : 1e-6;
    dt = std::min(dt, 0.1 * set_.spacing * set_.spacing / std::max(nu, 1e-12));
    return dt;
}

double Fluid::kineticEnergy() const {
    double e = 0;
    for (const auto& u : v) e += 0.5 * mass_ * length2(u);
    return e;
}

double Fluid::densityAt(const Vec3& p) const {
    double r = 0;
    grid_.query(x, p, H_, [&](int j) { r += mass_ * kernel_.W(length(x[j] - p)); });
    return r;
}

void Fluid::findWalls(const std::vector<Collider>& colliders) {
    const std::size_t n = x.size();
    touch_.assign(n * kMaxTouch, WallTouch{});
    nTouch_.assign(n, 0);
    staticTouch_.assign(n, 0);
    for (const auto& c : colliders) {
        if (!c.enabled) continue;
        for (std::size_t i = 0; i < n; ++i) {
            const Vec3& p = xp_[i];
            if (c.kind != ShapeKind::Plane && length2(p - c.pos) > sq(c.boundRadius + 2 * H_)) continue;
            Vec3 nn;
            const double d = c.sdf(p, &nn);
            if (c.body < 0 && d < 2 * H_) staticTouch_[i] = 1;
            if (d >= H_) continue;
            WallTouch* t = &touch_[i * kMaxTouch];
            int& k = nTouch_[i];
            if (k < kMaxTouch) t[k++] = {c.id, d, nn};
            else {
                int worst = 0;
                for (int q = 1; q < kMaxTouch; ++q)
                    if (t[q].d > t[worst].d) worst = q;
                if (d < t[worst].d) t[worst] = {c.id, d, nn};
            }
        }
    }
}

double Fluid::coupling(const RigidBody& b) const {
    // A body only a few particles heavy cannot be resolved by the liquid round
    // it: it feels the liquid in proportion to how many particles it outweighs.
    const double n = b.mass / mass_;
    return clamp((n - 20.0) / 180.0, 0.0, 1.0);
}

void Fluid::step(double dt, const std::vector<Collider>& colliders, std::vector<RigidBody>& bodies,
                 std::vector<int>* spray) {
    const std::size_t n = x.size();
    if (spray) spray->clear();
    if (n == 0) return;
    const Vec3 g{0, -kGravity, 0};
    xp_.resize(n);
    dx_.resize(n);
    lambda_.resize(n);
    normal_.resize(n);
    mu_.assign(n, liquid_->viscosity);
    rho.resize(n);
    fluidRho_.assign(n, rho0_);
    for (std::size_t i = 0; i < n; ++i) {
        v[i] += g * dt;
        xp_[i] = x[i] + v[i] * dt;
    }
    grid_.build(xp_, H_);
    grid_.neighbourLists(xp_, H_, off_, nb_);
    findWalls(colliders);

    const double mr = mass_ / rho0_;
    const double W0 = kernel_.W(0);
    const double eps = 1.0 / (set_.spacing * set_.spacing) * 1e-3;
    std::vector<Vec3> wallPush(n * kMaxTouch), wallGrad(n * kMaxTouch);

    for (int it = 0; it < set_.iterations; ++it) {
        double errSum = 0, errMax = 0;
        for (std::size_t i = 0; i < n; ++i) {
            double r = mass_ * W0;
            Vec3 gi;
            double sum2 = 0;
            for (int k = off_[i]; k < off_[i + 1]; ++k) {
                const int j = nb_[k];
                const Vec3 rij = xp_[i] - xp_[j];
                const double l = length(rij);
                r += mass_ * kernel_.W(l);
                const Vec3 gw = kernel_.gradW(rij, l) * mr;
                gi += gw;
                sum2 += length2(gw);
            }
            // Walls combine by inclusion-exclusion, Phi = 1 - prod(1 - Phi_k), so a
            // corner does not count the solid where two walls overlap twice.
            double phi[kMaxTouch], keep = 1;
            for (int t = 0; t < nTouch_[i]; ++t) {
                phi[t] = kernel_.wallVolume(touch_[i * kMaxTouch + t].d - wallShift_);
                keep *= 1 - phi[t];
            }
            r += rho0_ * (1 - keep);
            for (int t = 0; t < nTouch_[i]; ++t) {
                const WallTouch& wt = touch_[i * kMaxTouch + t];
                wallGrad[i * kMaxTouch + t] = wt.n * (kernel_.wallSection(wt.d - wallShift_) * others(phi, nTouch_[i], t));
                gi -= wallGrad[i * kMaxTouch + t];
            }
            fluidRho_[i] = r - rho0_ * (1 - keep);
            rho[i] = r;
            const double C = std::max(r / rho0_ - 1.0, 0.0);
            errSum += C;
            errMax = std::max(errMax, C);
            lambda_[i] = -C / (length2(gi) + sum2 + eps);
        }
        for (std::size_t i = 0; i < n; ++i) {
            Vec3 d;
            for (int k = off_[i]; k < off_[i + 1]; ++k) {
                const int j = nb_[k];
                Vec3 rij = xp_[i] - xp_[j];
                double l = length(rij);
                if (l < 1e-6 * H_) {
                    // Exactly on top of each other: part them along a fixed direction.
                    const double a = 2.399963 * static_cast<double>(std::min<std::size_t>(i, j) % 97);
                    rij = Vec3{std::cos(a), 0.3, std::sin(a)} * ((i < static_cast<std::size_t>(j) ? 1.0 : -1.0) * 0.05 * H_);
                    l = length(rij);
                }
                d += kernel_.gradW(rij, l) * (mr * (lambda_[i] + lambda_[j]));
            }
            for (int t = 0; t < nTouch_[i]; ++t) {
                const Vec3 push = wallGrad[i * kMaxTouch + t] * -lambda_[i];
                d += push;
                wallPush[i * kMaxTouch + t] += push;
            }
            // A correction is bounded per iteration (Macklin & Mueller's
            // practice), in length and as a speed (the constraint's own "sound
            // speed"): a particle found deep in a fast wall, or squeezed between
            // two shards, comes out over a few steps instead of taking a kick
            // that would grow as the step shrinks.
            const double cap = std::min(0.2 * set_.spacing, set_.pushSpeed * dt), dl = length(d);
            if (dl > cap) {
                const double f = cap / dl;
                d *= f;
                for (int t = 0; t < nTouch_[i]; ++t) {
                    wallPush[i * kMaxTouch + t] -= wallGrad[i * kMaxTouch + t] * (-lambda_[i] * (1 - f));
                }
            }
            dx_[i] = d;
        }
        for (std::size_t i = 0; i < n; ++i) {
            xp_[i] += dx_[i];
            for (int t = 0; t < nTouch_[i]; ++t) {
                WallTouch& wt = touch_[i * kMaxTouch + t];
                wt.d += dot(wt.n, dx_[i]);
            }
        }
        meanErr_ = errSum / static_cast<double>(n);
        maxErr_ = errMax;
    }

    // Hard contact: nothing ends closer than 0.4 spacing to a wall. The
    // velocity comes from the position before this projection, so a particle
    // pushed out of a wall does not fly off with the push.
    for (std::size_t i = 0; i < n; ++i) dx_[i] = xp_[i];
    const double rp = 0.4 * set_.spacing;
    std::vector<char> hit(n * kMaxTouch, 0);
    std::vector<Vec3> hitN(n * kMaxTouch), hitV(n * kMaxTouch);
    for (std::size_t i = 0; i < n; ++i) {
        for (int t = 0; t < nTouch_[i]; ++t) {
            const WallTouch& wt = touch_[i * kMaxTouch + t];
            const Collider& c = colliders[wt.collider];
            Vec3 nn;
            const double d = c.sdf(xp_[i], &nn);
            if (d < rp) {
                xp_[i] += nn * (rp - d);
                const std::size_t k = i * kMaxTouch + t;
                hit[k] = 1;
                hitN[k] = nn;
                hitV[k] = c.body >= 0 ? bodies[c.body].velocityAt(xp_[i]) : c.role == Role::Impactor ? c.vel : Vec3{};
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 vNew = (dx_[i] - x[i]) / dt;
        // The constraint's wall pushes are momentum taken from the wall's body.
        for (int t = 0; t < nTouch_[i]; ++t) {
            const WallTouch& wt = touch_[i * kMaxTouch + t];
            const Collider& c = colliders[wt.collider];
            if (c.body < 0) continue;
            const Vec3 J = wallPush[i * kMaxTouch + t] * (-mass_ / dt);
            RigidBody& b = bodies[c.body];
            const Vec3 at = xp_[i] - wt.n * wt.d;
            const double share = coupling(b);
            if (share <= 0) continue;
            b.fluidImpulse += J * share;
            b.fluidAngular += cross(at - b.x, J * share);
            if (wt.d < 0.5 * H_) b.fluidMass += mass_ * share;
        }
        v[i] = vNew;
        // Each wall the particle was pushed out of takes away the velocity
        // into it, and a little of the slip along it.
        for (int t = 0; t < nTouch_[i]; ++t) {
            const std::size_t k = i * kMaxTouch + t;
            if (!hit[k]) continue;
            const Vec3 rel = v[i] - hitV[k];
            const double vn = dot(rel, hitN[k]);
            Vec3 relN = rel;
            if (vn < 0) relN -= hitN[k] * vn;
            const Vec3 tang = relN - hitN[k] * dot(relN, hitN[k]);
            relN -= tang * set_.wallFriction;
            const Vec3 vOut = hitV[k] + relN;
            const Collider& c = colliders[touch_[k].collider];
            if (c.body >= 0) {
                const Vec3 J = (v[i] - vOut) * (mass_ * coupling(bodies[c.body]));
                bodies[c.body].fluidImpulse += J;
                bodies[c.body].fluidAngular += cross(xp_[i] - bodies[c.body].x, J);
            }
            v[i] = vOut;
        }
        // Contact-angle hysteresis pins a spreading edge: the force per length
        // sigma (cos theta_r - cos theta_a) on a particle's width of contact
        // line slows its slip along the surface, so a spill stops spreading
        // instead of thinning out for ever.
        if (nTouch_[i] > 0) {
            int nearest = 0;
            for (int t = 1; t < nTouch_[i]; ++t)
                if (touch_[i * kMaxTouch + t].d < touch_[i * kMaxTouch + nearest].d) nearest = t;
            const WallTouch& wt = touch_[i * kMaxTouch + nearest];
            const Collider& c = colliders[wt.collider];
            if (wt.d < 0.75 * set_.spacing && c.surface) {
                const Vec3 vc = c.body >= 0 ? bodies[c.body].velocityAt(xp_[i]) : c.role == Role::Impactor ? c.vel : Vec3{};
                const Vec3 rel = v[i] - vc;
                const Vec3 tang = rel - wt.n * dot(rel, wt.n);
                const double vt = length(tang);
                double pin = c.surface->hysteresis * liquid_->surfaceTension * set_.spacing / mass_ * dt;
                // Nearly alone on the surface, a particle is a landed drop: it
                // has spread to its stain and stays (spray.hpp).
                if (fluidRho_[i] < 0.55 * rho0_) pin = vt;
                if (vt > 1e-12) v[i] -= tang * (std::min(vt, pin) / vt);
            }
        }
        x[i] = xp_[i];
    }

    // Viscosity (with the local shear rate for shear-thinning liquids),
    // surface tension and XSPH, from the corrected positions.
    //
    // Surface tension is the continuum surface force (Brackbill, Kothe &
    // Zemach 1992, in SPH after Morris 2000): the colour gradient n = grad c
    // marks the surface, its curvature is kappa = -div(n / |n|), and the force
    // per unit volume is sigma kappa n. A flat surface feels nothing; a drop of
    // radius R is squeezed by 2 sigma / R.
    const bool thinning = liquid_->rheology == Rheology::Carreau;
    const double nMin = 0.1 / H_;
    std::vector<Vec3> unit(n);
    for (std::size_t i = 0; i < n; ++i) {
        Vec3 nrm;
        Mat3 L;
        for (int k = off_[i]; k < off_[i + 1]; ++k) {
            const int j = nb_[k];
            const Vec3 rij = x[i] - x[j];
            const Vec3 gw = kernel_.gradW(rij, length(rij));
            const double vj = mass_ / std::max(rho[j], 0.1 * rho0_);
            nrm += gw * vj;
            if (thinning) L = L + outer(v[j] - v[i], gw) * vj;
        }
        // grad c = sum V_j (c_j - c_i) grad W_ij with c = 1 in the liquid: the
        // missing neighbours outside count as c = 0.
        normal_[i] = -nrm;
        const double nl = length(normal_[i]);
        unit[i] = nl > nMin && nTouch_[i] == 0 ? normal_[i] / nl : Vec3{};
        if (thinning) {
            const Mat3 D = (L + L.transposed()) * 0.5;
            double dd = 0;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) dd += D.m[a][b] * D.m[a][b];
            mu_[i] = liquid_->viscosityAt(std::sqrt(2 * dd));
        }
    }
    const double sigma = liquid_->surfaceTension * set_.surfaceTension;
    for (std::size_t i = 0; i < n; ++i) {
        Vec3 acc, smooth;
        const double ri = std::max(rho[i], 0.1 * rho0_);
        double div = 0, wsum = 0;
        for (int k = off_[i]; k < off_[i + 1]; ++k) {
            const int j = nb_[k];
            const Vec3 rij = x[i] - x[j];
            const double l = length(rij);
            const double rj = std::max(rho[j], 0.1 * rho0_);
            const Vec3 gw = kernel_.gradW(rij, l);
            const Vec3 vij = v[i] - v[j];
            const double muij = 0.5 * (mu_[i] + mu_[j]);
            acc += gw * (10.0 * (mass_ / rj) * (muij / ri) * dot(vij, rij) / (l * l + 0.01 * H_ * H_));
            smooth -= vij * ((mass_ / rj) * kernel_.W(l));
            if (length2(unit[j]) > 0) {
                div += (mass_ / rj) * dot(unit[j] - unit[i], gw);
                wsum += (mass_ / rj) * kernel_.W(l);
            }
        }
        if (sigma > 0 && length2(unit[i]) > 0 && wsum > 0) {
            // Normalised by the neighbours that carry a normal (Morris 2000).
            const double kappa = -div / std::max(wsum, 0.3);
            acc += normal_[i] * (sigma * kappa / ri);
        }
        dx_[i] = acc * dt + smooth * set_.xsph;
    }
    for (std::size_t i = 0; i < n; ++i) v[i] += dx_[i];

    if (spray) {
        for (std::size_t i = 0; i < n; ++i) {
            if (nTouch_[i] > 0) {
                bool near = false;
                for (int t = 0; t < nTouch_[i]; ++t) near = near || touch_[i * kMaxTouch + t].d < 0.75 * H_;
                if (near) continue;
            }
            double r = mass_ * W0;
            for (int k = off_[i]; k < off_[i + 1]; ++k) r += mass_ * kernel_.W(length(x[i] - x[nb_[k]]));
            if (r < (selfDensity_ + set_.sprayDensity) * rho0_) spray->push_back(static_cast<int>(i));
        }
    }
}

}  // namespace lqd
