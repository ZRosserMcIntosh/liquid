// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/ragdoll.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace lqd {

namespace {
// Closest points of segments ab and cd: parameters s, t and the distance.
// An angle into [lo, hi] by the nearer end round the circle (never the far
// one: an angle just past +180 is just past +180, not near -180).
double clampAngle(double a, double lo, double hi) {
    if (a >= lo && a <= hi) return a;
    auto gap = [](double x, double y) {
        double d = std::fmod(std::abs(x - y), 2 * kPi);
        return d > kPi ? 2 * kPi - d : d;
    };
    return gap(a, lo) < gap(a, hi) ? lo : hi;
}

double segmentDistance(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, double& s, double& t) {
    const Vec3 u = b - a, v = d - c, w = a - c;
    const double A = dot(u, u), B = dot(u, v), C = dot(v, v), D = dot(u, w), E = dot(v, w);
    const double den = A * C - B * B;
    s = den > 1e-12 ? clamp((B * E - C * D) / den, 0.0, 1.0) : 0.0;
    t = C > 1e-12 ? clamp((B * s + E) / C, 0.0, 1.0) : 0.0;
    s = A > 1e-12 ? clamp((B * t - D) / A, 0.0, 1.0) : 0.0;
    return length(a + u * s - (c + v * t));
}
}  // namespace

void Ragdoll::build(const DummyDef& def, const std::vector<int>& ids, const std::vector<Collider>& all) {
    x_.clear();
    w_.clear();
    seg_.clear();
    links_.clear();
    limits_.clear();
    pairs_.clear();
    mass_ = 0;
    auto rigid = [&](int i, int j) {
        const double d = length(x_[i] - x_[j]);
        links_.push_back({i, j, d, d});
    };
    auto range = [&](int i, int j, double f) {
        const double d = length(x_[i] - x_[j]);
        links_.push_back({i, j, d * (1 - f), d * (1 + f)});
    };
    std::map<std::string, int> by;
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
        links_.push_back({s.i, s.j, std::max(L, 1e-4), std::max(L, 1e-4)});
        by[s.name] = static_cast<int>(k);
    }
    // A hinge folds no further than its two far ends half its length apart
    // (about 120 degrees): a link that only keeps them apart.
    auto fold = [&](int a, int f) {
        const double d = length(x_[f] - x_[a]);
        links_.push_back({a, f, 0.5 * d, 1e9, 0.3});
    };
    auto has = [&](const char* n) { return by.count(n) > 0; };
    auto S = [&](const char* n) -> const Seg& { return seg_[by.at(n)]; };
    if (has("left chest")) chestL_ = by["left chest"];
    if (has("right chest")) chestR_ = by["right chest"];
    if (has("pelvis")) pelvis_ = by["pelvis"];
    if (has("neck")) neck_ = by["neck"];
    if (def.mass > 0 && mass_ > 0) {
        // The capsules overlap at the joints: scale to the dummy's real mass.
        const double f = mass_ / def.mass;
        for (auto& w : w_) w *= f;
        mass_ = def.mass;
    }
    prev_ = x_;
    vel_.assign(x_.size(), {});
    const bool human = chestL_ >= 0 && chestR_ >= 0 && pelvis_ >= 0 && neck_ >= 0 && has("head");
    std::vector<std::pair<int, int>> joined;  // segments that meet at a joint
    if (human) {
        // The thorax (both chests) is rigid; the neck's base is fixed to it, its
        // top and the head nod and turn a little (about +-25 degrees); the
        // pelvis flexes against it at the lumbar spine (about +-15 degrees).
        const Seg &cl = S("left chest"), &cr = S("right chest"), &pv = S("pelvis"), &nk = S("neck"), &hd = S("head");
        const int th[4] = {cl.i, cl.j, cr.i, cr.j};
        for (int p = 0; p < 4; ++p)
            for (int q = p + 1; q < 4; ++q)
                if (th[p] / 2 != th[q] / 2) rigid(th[p], th[q]);
        for (int t : th) {
            rigid(nk.i, t);
            range(nk.j, t, 0.06);
            range(hd.i, t, 0.07);
            range(hd.j, t, 0.07);
            range(pv.i, t, 0.08);
            range(pv.j, t, 0.08);
        }
        rigid(nk.j, hd.i);
        rigid(nk.i, hd.j);
        rigid(nk.j, hd.j);
        rigid(nk.i, hd.i);
        for (const char* a : {"left chest", "right chest", "pelvis", "neck", "head"})
            for (const char* b : {"left chest", "right chest", "pelvis", "neck", "head"})
                if (std::string(a) < b) joined.push_back({by[a], by[b]});
        // Limbs: shoulders and hips pinned to the trunk, knees, elbows, ankles.
        auto pin = [&](int end, const Seg& host) {
            rigid(end, host.i);
            rigid(end, host.j);
        };
        const double deg = kPi / 180;
        for (int side = 0; side < 2; ++side) {
            const std::string L = side ? "right " : "left ";
            const double out = side ? 1 : -1;
            const Seg* chest = side ? &cr : &cl;
            if (by.count(L + "upper arm")) {
                const Seg& ua = seg_[by[L + "upper arm"]];
                pin(ua.i, *chest);
                joined.push_back({by[L + "upper arm"], side ? chestR_ : chestL_});
                joined.push_back({by[L + "upper arm"], side ? chestL_ : chestR_});
                // Flexion -50..170, abduction -30..150 (limited to 85 in this parametrisation).
                Limit c;
                c.near = ua.i;
                c.far = ua.j;
                c.outward = out;
                c.lo1 = -50 * deg;
                c.hi1 = 170 * deg;
                c.lo2 = -30 * deg;
                c.hi2 = 85 * deg;
                limits_.push_back(c);
                if (by.count(L + "forearm")) {
                    const Seg& fa = seg_[by[L + "forearm"]];
                    rigid(fa.i, ua.j);
                    joined.push_back({by[L + "upper arm"], by[L + "forearm"]});
                    // The elbow: 0..145 degrees forward; the forearm turns with the humerus (+-35 out of plane).
                    Limit h;
                    h.hinge = true;
                    h.near = fa.i;
                    h.far = fa.j;
                    h.parentA = ua.i;
                    h.parentB = ua.j;
                    h.bendSign = 1;
                    h.lo1 = -3 * deg;
                    h.hi1 = 145 * deg;
                    h.lo2 = -40 * deg;
                    h.hi2 = 40 * deg;
                    limits_.push_back(h);
                    fold(ua.i, fa.j);
                }
            }
            if (by.count(L + "thigh")) {
                const Seg& th2 = seg_[by[L + "thigh"]];
                pin(th2.i, pv);
                joined.push_back({by[L + "thigh"], pelvis_});
                // Flexion -20..120, abduction -20..45.
                Limit c;
                c.near = th2.i;
                c.far = th2.j;
                c.pelvisFrame = true;
                c.outward = out;
                c.lo1 = -20 * deg;
                c.hi1 = 120 * deg;
                c.lo2 = -20 * deg;
                c.hi2 = 45 * deg;
                limits_.push_back(c);
                if (by.count(L + "shin")) {
                    const Seg& sh = seg_[by[L + "shin"]];
                    rigid(sh.i, th2.j);
                    joined.push_back({by[L + "thigh"], by[L + "shin"]});
                    // The knee: 0..140 degrees backward, a few degrees of play out of plane.
                    Limit h;
                    h.hinge = true;
                    h.near = sh.i;
                    h.far = sh.j;
                    h.parentA = th2.i;
                    h.parentB = th2.j;
                    h.pelvisFrame = true;
                    h.bendSign = -1;
                    h.lo1 = -3 * deg;
                    h.hi1 = 140 * deg;
                    h.lo2 = -25 * deg;  // the thigh cannot twist here: the hip's rotation shows as play
                    h.hi2 = 25 * deg;
                    limits_.push_back(h);
                    fold(th2.i, sh.j);
                    if (by.count(L + "foot")) {
                        const Seg& ft = seg_[by[L + "foot"]];
                        // The ankle: the heel hangs off the shin's end, the toe rocks
                        // through about +-25 degrees.
                        rigid(sh.j, ft.i);
                        range(sh.i, ft.j, 0.08);
                        range(sh.i, ft.i, 0.02);
                        joined.push_back({by[L + "shin"], by[L + "foot"]});
                    }
                }
            }
        }
    } else {
        // Any other figure: the nearest ends of two capsules within 0.2 m are joined.
        for (std::size_t A = 0; A < seg_.size(); ++A)
            for (std::size_t B = A + 1; B < seg_.size(); ++B) {
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
                rigid(bi, bj);
                joined.push_back({static_cast<int>(A), static_cast<int>(B)});
            }
    }
    // Limbs collide with the rest of the body (not with what they are joined to).
    for (int A = 0; A < static_cast<int>(seg_.size()); ++A)
        for (int B = A + 1; B < static_cast<int>(seg_.size()); ++B) {
            if (std::find(joined.begin(), joined.end(), std::make_pair(A, B)) != joined.end() ||
                std::find(joined.begin(), joined.end(), std::make_pair(B, A)) != joined.end())
                continue;
            double s, t;
            const double d = segmentDistance(x_[seg_[A].i], x_[seg_[A].j], x_[seg_[B].i], x_[seg_[B].j], s, t);
            pairs_.push_back({A, B, std::min(seg_[A].r + seg_[B].r, d)});
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
        const double target = clamp(len, l.lo, l.hi);
        if (target == len) continue;
        const double ws = w_[l.i] + w_[l.j];
        if (ws <= 0) continue;
        const Vec3 corr = d * (l.stiffness * (len - target) / (len * ws));
        x_[l.i] += corr * w_[l.i];
        x_[l.j] -= corr * w_[l.j];
    }
}

void Ragdoll::thorax(Vec3& up, Vec3& lat, Vec3& fwd) const {
    const Seg &l = seg_[chestL_], &r = seg_[chestR_];
    lat = normalize(x_[r.i] + x_[r.j] - x_[l.i] - x_[l.j], {1, 0, 0});
    up = x_[l.j] + x_[r.j] - x_[l.i] - x_[r.i];
    up = normalize(up - lat * dot(up, lat), {0, 1, 0});
    fwd = cross(lat, up);
}

void Ragdoll::pelvisFrame(Vec3& up, Vec3& lat, Vec3& fwd) const {
    Vec3 tu, tl, tf;
    thorax(tu, tl, tf);
    const Seg& p = seg_[pelvis_];
    lat = normalize(x_[p.j] - x_[p.i], tl);
    up = normalize(tu - lat * dot(tu, lat), tu);
    fwd = cross(lat, up);
}

void Ragdoll::solveLimits() {
    if (chestL_ < 0 || chestR_ < 0) return;
    for (const auto& L : limits_) {
        const Vec3 N = x_[L.near];
        Vec3 d = x_[L.far] - N;
        const double len = length(d);
        if (len < 1e-9) continue;
        d = d / len;
        Vec3 nd;
        if (!L.hinge) {
            Vec3 up, lat, fwd;
            if (L.pelvisFrame) pelvisFrame(up, lat, fwd);
            else thorax(up, lat, fwd);
            const Vec3 down = -up, out = lat * L.outward;
            const double c = clamp(dot(d, out), -1.0, 1.0);
            const double phi = std::atan2(dot(d, fwd), dot(d, down)), psi = std::asin(c);
            const double p2 = clampAngle(phi, L.lo1, L.hi1), s2 = clamp(psi, L.lo2, L.hi2);
            if (p2 == phi && s2 == psi) continue;
            nd = (down * std::cos(p2) + fwd * std::sin(p2)) * std::cos(s2) + out * std::sin(s2);
        } else {
            const Vec3 e = normalize(x_[L.parentB] - x_[L.parentA], d);
            // The hinge turns about the left-right axis of what the limb hangs
            // from (the pelvis for a knee, the thorax for an elbow), so it bends
            // the right way however the body lies.
            Vec3 up, lat, fwd;
            if (L.pelvisFrame) pelvisFrame(up, lat, fwd);
            else thorax(up, lat, fwd);
            Vec3 bend = cross(lat, e) * -L.bendSign;
            if (length2(bend) < 0.04) {
                // The limb lies along that axis (an arm held out sideways): fall back on the facing.
                bend = fwd * L.bendSign;
                bend = bend - e * dot(bend, e);
            }
            if (length2(bend) < 1e-8) continue;
            bend = normalize(bend);
            const Vec3 axis = cross(e, bend);
            const double c = clamp(dot(d, axis), -1.0, 1.0);
            // Softly (hard hinges fight the floor and the other joints and pump
            // energy in): not backwards past straight, nor far out of its
            // plane. How far it may fold is a link (fold, in build()).
            const double beta = std::atan2(dot(d, bend), dot(d, e)), psi = std::asin(c);
            const double b2 = beta < L.lo1 && beta > -kPi / 2 ? L.lo1 : beta, s2 = clamp(psi, L.lo2, L.hi2);
            if (b2 == beta && s2 == psi) continue;
            nd = (e * std::cos(b2) + bend * std::sin(b2)) * std::cos(s2) + axis * std::sin(s2);
            const Vec3 soft = (N + nd * len - x_[L.far]) * 0.15;
            const double wf = w_[L.far], wa = w_[L.parentA], ws2 = wf + wa;
            if (ws2 <= 0) continue;
            // Shared with the parent's far end, so the chain keeps its momentum.
            x_[L.far] += soft * (wf / ws2);
            x_[L.parentA] -= soft * (wa / ws2);
            continue;
        }
        // A third of the way back into range each pass, mostly by the child's far end.
        const Vec3 corr = (N + nd * len - x_[L.far]) * 0.35;
        const double wf = w_[L.far], wn = w_[L.near], ws = wf + wn;
        if (ws <= 0) continue;
        x_[L.far] += corr * (wf / ws + 0.5 * wn / ws);
        x_[L.near] -= corr * (0.5 * wn / ws);
    }
}

void Ragdoll::selfCollide() {
    for (const auto& P : pairs_) {
        const Seg &A = seg_[P.a], &B = seg_[P.b];
        double s, t;
        const double d = segmentDistance(x_[A.i], x_[A.j], x_[B.i], x_[B.j], s, t);
        if (d >= P.gap || d < 1e-9) continue;
        const Vec3 pa = x_[A.i] + (x_[A.j] - x_[A.i]) * s, pb = x_[B.i] + (x_[B.j] - x_[B.i]) * t;
        const Vec3 n = (pa - pb) / d;
        const double wa = (1 - s) * (1 - s) * w_[A.i] + s * s * w_[A.j], wb = (1 - t) * (1 - t) * w_[B.i] + t * t * w_[B.j];
        if (wa + wb <= 0) continue;
        const double lam = 0.5 * (P.gap - d) / (wa + wb);  // softly: a limit may be pressing them together
        x_[A.i] += n * (lam * (1 - s) * w_[A.i]);
        x_[A.j] += n * (lam * s * w_[A.j]);
        x_[B.i] -= n * (lam * (1 - t) * w_[B.i]);
        x_[B.j] -= n * (lam * t * w_[B.j]);
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
                    // An inelastic contact: no speed along the normal survives
                    // being pushed out (else the push itself is a bounce, and
                    // with the joints pressing back, a growing one).
                    const double vn = dot(x_[idx] - prev_[idx], n);
                    prev_[idx] += n * vn;
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
    // The body's mechanical energy (kinetic and potential, particle masses 1/w).
    auto energy = [&](const std::vector<Vec3>* v) {
        double e = 0;
        for (std::size_t i = 0; i < x_.size(); ++i) {
            if (w_[i] <= 0) continue;
            const double m = 1 / w_[i];
            const Vec3 u = v ? (*v)[i] : (x_[i] - prev_[i]) / h;
            e += 0.5 * m * dot(u, u) + m * kGravity * x_[i].y;
        }
        return e;
    };
    std::vector<Vec3> v0(x_.size());
    for (int sub = 0; sub < n; ++sub) {
        // Energy in: what it had, and what a hit has just given it.
        for (std::size_t i = 0; i < x_.size(); ++i) v0[i] = (x_[i] - prev_[i]) / h * 0.9995 + vel_[i];
        const double e0 = energy(&v0);
        for (std::size_t i = 0; i < x_.size(); ++i) {
            Vec3 v = (x_[i] - prev_[i]) / h * 0.9995 + vel_[i];  // a little air drag
            // Constraints that fight (a joint limit against a limb in the way)
            // must not pump energy in: no part of a falling body moves at 15 m/s.
            const double sp = length(v);
            if (sp > 15) v = v * (15 / sp);
            vel_[i] = {};
            prev_[i] = x_[i];
            x_[i] += v * h + g * (h * h);
        }
        for (int it = 0; it < 10; ++it) {
            solveLinks();
            solveLimits();
            selfCollide();
            collide(colliders);
        }
        // A limp body with friction only loses energy: if the constraints,
        // pressing against each other and the floor, made some, take it back
        // out of the motion (position-based constraints can pump energy in
        // where they conflict; this keeps the sum honest).
        const double e1 = energy(nullptr);
        if (e1 > e0 + 1e-9) {
            double ke = 0, pe = 0;
            for (std::size_t i = 0; i < x_.size(); ++i) {
                if (w_[i] <= 0) continue;
                const double m = 1 / w_[i];
                const Vec3 u = (x_[i] - prev_[i]) / h;
                ke += 0.5 * m * dot(u, u);
                pe += m * kGravity * x_[i].y;
            }
            const double f = ke > 1e-12 ? std::sqrt(std::max(0.0, e0 - pe) / ke) : 0.0;
            for (std::size_t i = 0; i < x_.size(); ++i) prev_[i] = x_[i] - (x_[i] - prev_[i]) * f;
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
