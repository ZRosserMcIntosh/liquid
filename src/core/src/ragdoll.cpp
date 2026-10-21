// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/ragdoll.hpp"

#include <algorithm>
#include <cmath>

namespace lqd {

void Ragdoll::build(const DummyDef& def, const std::vector<int>& ids, const std::vector<Collider>& all) {
    x_.clear();
    w_.clear();
    seg_.clear();
    links_.clear();
    mass_ = 0;
    for (std::size_t k = 0; k < ids.size(); ++k) {
        const Collider& c = all[ids[k]];
        const double L = length(c.b - c.a), r = c.radius;
        const double m = def.tissueDensity * (kPi * r * r * L + 4.0 / 3.0 * kPi * r * r * r);
        mass_ += m;
        Seg s;
        s.i = static_cast<int>(x_.size());
        s.j = s.i + 1;
        s.r = r;
        s.collider = ids[k];
        s.name = c.name;
        x_.push_back(c.a);
        x_.push_back(c.b);
        w_.push_back(2.0 / m);
        w_.push_back(2.0 / m);
        seg_.push_back(s);
        links_.push_back({s.i, s.j, std::max(L, 1e-4), false});
        if (s.name == "left chest") chestL_ = static_cast<int>(k);
        if (s.name == "right chest") chestR_ = static_cast<int>(k);
        if (s.name == "pelvis") pelvis_ = static_cast<int>(k);
        if (s.name == "neck") neck_ = static_cast<int>(k);
    }
    if (def.mass > 0 && mass_ > 0) {
        // The capsules overlap at the joints: scale to the dummy's real mass.
        const double f = mass_ / def.mass;
        for (auto& w : w_) w *= f;
        mass_ = def.mass;
    }
    prev_ = x_;
    vel_.assign(x_.size(), {});
    auto trunk = [](const std::string& n) {
        return n == "left chest" || n == "right chest" || n == "pelvis" || n == "neck" || n == "head";
    };
    // The trunk is rigid: every pair of its points held.
    std::vector<int> tp;
    for (const auto& s : seg_)
        if (trunk(s.name)) {
            tp.push_back(s.i);
            tp.push_back(s.j);
        }
    for (std::size_t p = 0; p < tp.size(); ++p)
        for (std::size_t q = p + 1; q < tp.size(); ++q)
            if (!(tp[p] / 2 == tp[q] / 2)) links_.push_back({tp[p], tp[q], length(x_[tp[p]] - x_[tp[q]]), false});
    // Joints where capsules meet: the nearest ends of two capsules within 0.2 m.
    for (std::size_t A = 0; A < seg_.size(); ++A)
        for (std::size_t B = A + 1; B < seg_.size(); ++B) {
            if (trunk(seg_[A].name) && trunk(seg_[B].name)) continue;
            int bi = -1, bj = -1;
            double best = 0.2;
            for (int pa : {seg_[A].i, seg_[A].j})
                for (int pb : {seg_[B].i, seg_[B].j}) {
                    const double d = length(x_[pa] - x_[pb]);
                    if (d < best) {
                        best = d;
                        bi = pa;
                        bj = pb;
                    }
                }
            if (bi < 0) continue;
            links_.push_back({bi, bj, best, false});
            // A gap (shoulder, hip): pin the limb's end to the other capsule's line too.
            if (best > 0.01) {
                const bool aTrunk = trunk(seg_[A].name);
                const int limbEnd = aTrunk ? bj : bi;
                const Seg& host = aTrunk ? seg_[A] : seg_[B];
                for (int h : {host.i, host.j}) links_.push_back({limbEnd, h, length(x_[limbEnd] - x_[h]), false});
            } else {
                // A knee or an elbow: it may bend, not fold flat.
                const int farA = bi == seg_[A].i ? seg_[A].j : seg_[A].i, farB = bj == seg_[B].i ? seg_[B].j : seg_[B].i;
                links_.push_back({farA, farB, 0.55 * length(x_[farA] - x_[farB]), true});
            }
        }
    released_ = false;
    releaseTime_ = -1;
}

void Ragdoll::release(double time) {
    if (released_ || seg_.empty()) return;
    released_ = true;
    releaseTime_ = time;
}

int Ragdoll::segmentOf(int collider) const {
    for (int s = 0; s < static_cast<int>(seg_.size()); ++s)
        if (seg_[s].collider == collider) return s;
    return -1;
}

void Ragdoll::impulse(int s, const Vec3& at, const Vec3& J) {
    if (s < 0 || s >= static_cast<int>(seg_.size())) return;
    const Seg& g = seg_[s];
    const Vec3 A = x_[g.i], B = x_[g.j];
    const double t = clamp(dot(at - A, B - A) / std::max(length2(B - A), 1e-12), 0.0, 1.0);
    still_ = 0;
    // Two point masses on a lever: momentum J shared as (1 - t) and t.
    vel_[g.i] += J * ((1 - t) * w_[g.i]);
    vel_[g.j] += J * (t * w_[g.j]);
}

void Ragdoll::solveLinks() {
    for (const auto& l : links_) {
        const Vec3 d = x_[l.j] - x_[l.i];
        const double len = length(d);
        if (len < 1e-12) continue;
        if (l.minOnly && len >= l.rest) continue;
        const double ws = w_[l.i] + w_[l.j];
        if (ws <= 0) continue;
        const Vec3 corr = d * ((len - l.rest) / (len * ws));
        x_[l.i] += corr * w_[l.i];
        x_[l.j] -= corr * w_[l.j];
    }
}

void Ragdoll::collide(const std::vector<Collider>& colliders) {
    for (const auto& c : colliders) {
        if (!c.enabled || c.role == Role::Dummy || c.role == Role::Shard || c.body >= 0) continue;
        const double mu = c.surface ? c.surface->friction : 0.5;
        for (const auto& s : seg_) {
            // The two ends and the middle of the capsule.
            for (int k = 0; k < 3; ++k) {
                const Vec3 p = k == 0 ? x_[s.i] : k == 1 ? x_[s.j] : (x_[s.i] + x_[s.j]) * 0.5;
                if (c.kind != ShapeKind::Plane && length(p - c.pos) > c.boundRadius + s.r) continue;
                Vec3 n;
                const double d = c.sdf(p, &n) - s.r;
                if (d >= 0) continue;
                auto push = [&](int idx, double share) {
                    x_[idx] += n * (-d * share);
                    // Coulomb friction on the slip this step.
                    const Vec3 slip = x_[idx] - prev_[idx];
                    const Vec3 tang = slip - n * dot(slip, n);
                    const double tl = length(tang), cap = mu * (-d * share);
                    if (tl > 1e-12) prev_[idx] += tang * (std::min(tl, cap) / tl);
                };
                if (k == 0) push(s.i, 1);
                else if (k == 1) push(s.j, 1);
                else {
                    push(s.i, 1);
                    push(s.j, 1);
                }
            }
        }
    }
}

void Ragdoll::step(double dt, const std::vector<Collider>& colliders) {
    if (seg_.empty()) return;
    if (!released_) {
        // On the stand: still, but a round's push is held by the stand.
        std::fill(vel_.begin(), vel_.end(), Vec3{});
        prev_ = x_;
        return;
    }
    const int n = std::clamp(static_cast<int>(std::ceil(dt * 480)), 1, 8);
    const double h = dt / n;
    h_ = h;
    const Vec3 g{0, -kGravity, 0};
    for (int sub = 0; sub < n; ++sub) {
        for (std::size_t i = 0; i < x_.size(); ++i) {
            const Vec3 v = (x_[i] - prev_[i]) / h * 0.9995 + vel_[i];  // a little air drag
            vel_[i] = {};
            prev_[i] = x_[i];
            x_[i] += v * h + g * (h * h);
        }
        for (int it = 0; it < 10; ++it) {
            solveLinks();
            collide(colliders);
        }
    }
    // Lying still on the floor, static friction holds it: below 3 cm/s for
    // half a second, it settles.
    double vmax = 0;
    for (std::size_t i = 0; i < x_.size(); ++i) vmax = std::max(vmax, length(x_[i] - prev_[i]) / h);
    still_ = vmax < 0.03 ? still_ + dt : 0;
    if (still_ > 0.5) prev_ = x_;
}

bool Ragdoll::moving() const {
    if (!released_ || still_ > 0.5) return false;
    for (std::size_t i = 0; i < x_.size(); ++i)
        if (length2(x_[i] - prev_[i]) > sq(0.02 * h_)) return true;
    return false;
}

void Ragdoll::pose(std::vector<Collider>& colliders) const {
    for (const auto& s : seg_) {
        Collider& c = colliders[s.collider];
        c.a = x_[s.i];
        c.b = x_[s.j];
        c.pos = (c.a + c.b) * 0.5;
        c.boundRadius = 0.5 * length(c.b - c.a) + c.radius;
    }
}

Vec3 Ragdoll::facing() const {
    if (chestL_ < 0 || chestR_ < 0) return {0, 0, 1};
    const Vec3 lateral = x_[seg_[chestR_].i] - x_[seg_[chestL_].i];
    const Vec3 spine = x_[seg_[chestL_].j] - x_[seg_[chestL_].i];
    return normalize(cross(lateral, spine), {0, 0, 1});
}

void Ragdoll::frame(int s, Vec3& e, Vec3& ref, Vec3& side) const {
    e = normalize(x_[seg_[s].j] - x_[seg_[s].i], {0, 1, 0});
    Vec3 f = facing();
    f = f - e * dot(f, e);
    if (length2(f) < 1e-6) f = perpendicular(e);
    ref = normalize(f);
    side = cross(e, ref);
}

Attachment Ragdoll::attach(int s, const Vec3& p, const Vec3& nrm) const {
    Attachment a;
    a.segment = s;
    if (s < 0) return a;
    Vec3 e, ref, side;
    frame(s, e, ref, side);
    const Vec3 A = x_[seg_[s].i];
    const double L = std::max(length(x_[seg_[s].j] - A), 1e-9);
    const double along = dot(p - A, e);
    a.t = along / L;
    const Vec3 off = p - A - e * along;
    a.radial = length(off);
    a.angle = std::atan2(dot(off, side), dot(off, ref));
    a.normalAxial = dot(nrm, e);
    const Vec3 nn = nrm - e * a.normalAxial;
    a.normalAngle = std::atan2(dot(nn, side), dot(nn, ref));
    return a;
}

void Ragdoll::place(const Attachment& a, Vec3& p, Vec3& nrm) const {
    if (a.segment < 0) return;
    Vec3 e, ref, side;
    frame(a.segment, e, ref, side);
    const Vec3 A = x_[seg_[a.segment].i];
    const double L = length(x_[seg_[a.segment].j] - A);
    p = A + e * (a.t * L) + (ref * std::cos(a.angle) + side * std::sin(a.angle)) * a.radial;
    const double r = std::sqrt(std::max(0.0, 1 - a.normalAxial * a.normalAxial));
    nrm = normalize(e * a.normalAxial + (ref * std::cos(a.normalAngle) + side * std::sin(a.normalAngle)) * r);
}

Vec3 Ragdoll::velocityOf(const Attachment& a) const {
    if (a.segment < 0 || !released_) return {};
    const Seg& s = seg_[a.segment];
    const Vec3 vi = (x_[s.i] - prev_[s.i]) / h_, vj = (x_[s.j] - prev_[s.j]) / h_;
    return vi * (1 - a.t) + vj * a.t;
}

}  // namespace lqd
