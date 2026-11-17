// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/rigid.hpp"

#include <algorithm>
#include <cmath>

namespace lqd {

namespace {

double combinedModulus(double E1, double nu1, double E2, double nu2) {
    return 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
}

}  // namespace

ImpactAssessment assessImpact(const VesselDef& v, const WallSample& s, const SurfaceDef& surface, double meff,
                              double speed, double alongNormal) {
    ImpactAssessment a;
    a.speed = speed;
    a.effectiveMass = meff;
    a.edge = s.edge;
    const double energy = 0.5 * meff * speed * speed;
    const double R = std::max(s.radius, 1e-4);
    const double t = std::max(s.t, 2e-4);
    const double cn = clamp(alongNormal, 0.0, 1.0);
    // The wall's point-load stiffness for the normal part of the strike; struck
    // end-on the wall is stiff in its own plane and only the contact yields.
    const double kShell = 4 * v.youngs * t * t / (std::max(s.shellRadius, t) * std::sqrt(3 * (1 - v.poisson * v.poisson)));
    const double kWall = kShell / std::max(cn * cn, 1e-4);
    const double Es = combinedModulus(v.youngs, v.poisson, surface.youngs, surface.poisson);
    // Energy stored at force F in the contact (Hertz or Winkler layer) and the wall in series.
    double layerSqueeze = 0;
    auto stored = [&](double F) {
        double e = F * F / (2 * kWall);
        if (surface.layer > 0) {
            // The pile spreads a point load over about its own thickness.
            const double Rw = std::max(R, surface.layer);
            const double k = surface.youngs / surface.layer;
            const double dl = std::sqrt(F / (k * kPi * Rw));  // F = k pi R d^2
            if (dl <= surface.layer) {
                layerSqueeze = dl;
                return e + k * kPi * Rw * dl * dl * dl / 3;
            }
            // Bottomed out: the layer's full energy, the rest in the substrate.
            layerSqueeze = surface.layer;
            const double Fl = k * kPi * Rw * surface.layer * surface.layer;
            const double Esub = combinedModulus(v.youngs, v.poisson, 30e9, 0.2);
            const double extra = std::max(0.0, F - Fl);
            const double dh = std::pow(extra / (4.0 / 3.0 * Esub * std::sqrt(R)), 2.0 / 3.0);
            return e + k * kPi * Rw * std::pow(surface.layer, 3) / 3 + 0.4 * extra * dh;
        }
        const double dh = std::pow(F / (4.0 / 3.0 * Es * std::sqrt(R)), 2.0 / 3.0);
        layerSqueeze = dh;
        return e + 0.4 * F * dh;  // Hertz: energy = (2/5) F delta
    };
    double lo = 0, hi = 1;
    while (stored(hi) < energy && hi < 1e9) hi *= 2;
    for (int it = 0; it < 80; ++it) {
        const double mid = 0.5 * (lo + hi);
        (stored(mid) < energy ? lo : hi) = mid;
    }
    const double force = 0.5 * (lo + hi);
    stored(force);
    a.force = force;
    a.indentation = layerSqueeze;
    a.contactRadius = surface.layer > 0 ? std::sqrt(2.0 * std::max(R, surface.layer) * layerSqueeze)
                                        : std::sqrt(R * std::min(layerSqueeze, R));
    // Bending under the normal part of the force spread over the contact
    // (Timoshenko & Woinowsky-Krieger, load over a small circle), the
    // logarithmic term cut off at the bending reach; membrane stress from the
    // part along the wall, spread over the reach.
    const double ac = std::max(a.contactRadius, 0.5 * t);
    const double reach = std::max(s.reach, ac);
    double bend = 3.0 * (1 + v.poisson) * force * cn / (2 * kPi * t * t) * (1.0 + std::log(reach / ac));
    if (s.edge) bend *= 1.3;  // a free edge has no wall on one side to share the load
    const double membrane = force * std::sqrt(std::max(0.0, 1 - cn * cn)) / (t * 2 * reach);
    a.stress = bend + membrane;
    a.probability = 1.0 - std::exp(-std::pow(a.stress / v.strength, v.weibullModulus));
    return a;
}

void stepBody(RigidBody& b, const std::vector<Collider>& colliders, double dt, double strikeSpeed,
              std::vector<Strike>* strikes) {
    const double fluidKick = length(b.fluidImpulse) / b.mass;
    if (b.sleeping && fluidKick < 0.05) {
        b.fluidImpulse = b.fluidAngular = {};
        return;
    }
    b.sleeping = false;
    // The liquid's push, relaxed by the liquid in contact (its added mass):
    // the coupling is explicit, and a light shard under a heavy puddle would
    // otherwise take the whole puddle's correction in one step.
    const double relax = b.mass / (b.mass + b.fluidMass);
    b.v += Vec3{0, -kGravity, 0} * dt + b.fluidImpulse * (relax / b.mass);
    b.w += b.invInertiaWorld() * b.fluidAngular * relax;
    b.fluidImpulse = b.fluidAngular = {};
    b.fluidMass = 0;

    double reachR = 0;
    for (const auto& s : b.samples) reachR = std::max(reachR, length(s.p - b.comLocal));
    const double fast = length(b.v) + length(b.w) * reachR;
    const int nsub = std::clamp(static_cast<int>(std::ceil(fast * dt / 0.0008)), 1, 24);
    const double h = dt / nsub;

    struct Contact {
        Vec3 r, n;
        double d, vn0, jn = 0, jt = 0, mu, e;
        int collider, sample;
    };
    std::vector<Contact> contacts;
    for (int sub = 0; sub < nsub; ++sub) {
        contacts.clear();
        const Mat3 Iw = b.invInertiaWorld();
        for (int k = 0; k < static_cast<int>(b.samples.size()); ++k) {
            const Vec3 p = b.fromVessel(b.samples[k].p);
            for (const auto& c : colliders) {
                if (c.body >= 0 || !c.enabled) continue;
                if (c.kind != ShapeKind::Plane && length2(p - c.pos) > sq(c.boundRadius + 0.05)) continue;
                Vec3 n;
                const double d = c.sdf(p, &n);
                const Vec3 vp = b.velocityAt(p);
                const double vn = dot(vp, n);
                const double margin = std::max(0.0, -vn) * h + 2e-4;
                if (d > margin) continue;
                Contact ct;
                ct.r = p - b.x;
                ct.n = n;
                ct.d = d;
                ct.vn0 = vn;
                ct.mu = c.surface ? c.surface->friction : 0.4;
                ct.e = (c.surface ? c.surface->restitution : 0.3) * b.restitution;
                ct.collider = c.id;
                ct.sample = k;
                contacts.push_back(ct);
            }
        }
        if (strikes && !contacts.empty()) {
            // The fastest arriving point per collider is the strike.
            for (const auto& ct : contacts) {
                // A strike: the point arrives within this substep, at the speed
                // it has before the solver touches it.
                if (ct.d > -ct.vn0 * h + 1e-4 || -ct.vn0 < strikeSpeed) continue;
                bool dup = false;
                for (auto& s : *strikes)
                    if (s.body == b.id && s.collider == ct.collider) {
                        dup = true;
                        if (-ct.vn0 > s.speed) {
                            s.speed = -ct.vn0;
                            s.sample = ct.sample;
                            s.point = b.x + ct.r;
                            s.normal = ct.n;
                        }
                    }
                if (!dup) strikes->push_back({b.id, ct.collider, ct.sample, b.x + ct.r, ct.n, -ct.vn0});
            }
        }
        // The deepest-approach contact gets the bounce; the others just stop.
        double fastest = 0;
        for (const auto& ct : contacts) fastest = std::min(fastest, ct.vn0);
        for (int it = 0; it < 10 && !contacts.empty(); ++it) {
            for (auto& ct : contacts) {
                const Vec3 vrel = b.v + cross(b.w, ct.r);
                const double vn = dot(vrel, ct.n);
                double target;
                if (ct.d > 0) target = -ct.d / h;  // speculative: may close the gap, not more
                else {
                    target = 0;  // penetration is removed by moving the body, below
                    if (ct.vn0 < -0.3 && ct.vn0 <= fastest * 0.7) target = -ct.e * ct.vn0;
                }
                const Vec3 rn = cross(ct.r, ct.n);
                const double kn = 1.0 / b.mass + dot(rn, Iw * rn);
                const double jnew = std::max(ct.jn + (target - vn) / kn, 0.0);
                const double dj = jnew - ct.jn;
                ct.jn = jnew;
                b.v += ct.n * (dj / b.mass);
                b.w += Iw * (rn * dj);
                // Friction.
                const Vec3 v2 = b.v + cross(b.w, ct.r);
                Vec3 vt = v2 - ct.n * dot(v2, ct.n);
                const double vtl = length(vt);
                if (vtl > 1e-9) {
                    const Vec3 tdir = vt / vtl;
                    const Vec3 rt = cross(ct.r, tdir);
                    const double kt = 1.0 / b.mass + dot(rt, Iw * rt);
                    const double lim = ct.mu * ct.jn;
                    const double jtNew = clamp(ct.jt - vtl / kt, -lim, lim);
                    const double djt = jtNew - ct.jt;
                    ct.jt = jtNew;
                    b.v += tdir * (djt / b.mass);
                    b.w += Iw * (rt * djt);
                }
            }
        }
        b.x += b.v * h;
        b.q = b.q.integrated(b.w, h);
        // Positional correction: out of each collider by most of the deepest
        // penetration, without adding velocity.
        for (const auto& c : colliders) {
            if (c.body >= 0 || !c.enabled) continue;
            double deepest = 0;
            Vec3 dn;
            bool any = false;
            for (const auto& ct : contacts) any = any || ct.collider == c.id;
            if (!any) continue;
            for (const auto& smp : b.samples) {
                Vec3 n;
                const double d = c.sdf(b.fromVessel(smp.p), &n);
                if (-d > deepest) {
                    deepest = -d;
                    dn = n;
                }
            }
            if (deepest > 1e-4) b.x += dn * (0.8 * (deepest - 1e-4));
        }
    }
    // Sleep once it has been still for a while on something.
    if (!contacts.empty() && length(b.v) < 0.03 && length(b.w) < 0.6) b.still += dt;
    else b.still = 0;
    if (b.still > 0.3) {
        b.sleeping = true;
        b.v = b.w = {};
    }
}

}  // namespace lqd
