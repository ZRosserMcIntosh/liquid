// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/collider.hpp"

#include <cmath>

namespace lqd {

const char* roleName(Role r) {
    switch (r) {
        case Role::Floor: return "floor";
        case Role::Table: return "table";
        case Role::Wall: return "wall";
        case Role::Pane: return "pane";
        case Role::Dummy: return "dummy";
        case Role::Vessel: return "vessel";
        case Role::Shard: return "shard";
        default: return "other";
    }
}

int PaneMap::pieceAt(double x, double y, bool frontFace) const {
    if (nx <= 0 || ny <= 0 || cell <= 0) return -1;
    const int i = static_cast<int>(std::floor((x - x0) / cell)), j = static_cast<int>(std::floor((y - y0) / cell));
    if (i < 0 || j < 0 || i >= nx || j >= ny) return -1;
    const std::vector<int>& m = frontFace || back.empty() ? front : back;
    return m.empty() ? -1 : m[static_cast<std::size_t>(j) * nx + i];
}

bool PaneMap::open(double x, double y) const {
    if (laminated || nx <= 0) return false;
    for (int f = 0; f < 2; ++f) {
        const int k = pieceAt(x, y, f == 0);
        if (k >= 0 && k < static_cast<int>(held.size()) && held[k]) return false;
    }
    return true;
}

namespace {
// Signed distance in the plane to a polygon (x, y pairs; either winding) and
// the outward unit direction from its nearest edge.
double polygonDistance(const std::vector<double>& poly, double px, double py, double& nx, double& ny) {
    const std::size_t n = poly.size() / 2;
    double best = 1e30;
    bool inside = false;
    nx = 1;
    ny = 0;
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const double ax = poly[2 * j], ay = poly[2 * j + 1], bx = poly[2 * i], by = poly[2 * i + 1];
        const double ex = bx - ax, ey = by - ay, wx = px - ax, wy = py - ay;
        const double t = clamp((wx * ex + wy * ey) / std::max(ex * ex + ey * ey, 1e-300), 0.0, 1.0);
        const double dx = wx - ex * t, dy = wy - ey * t, d2 = dx * dx + dy * dy;
        if (d2 < best) {
            best = d2;
            const double l = std::sqrt(d2);
            if (l > 1e-12) {
                nx = dx / l;
                ny = dy / l;
            } else {
                const double el = std::sqrt(ex * ex + ey * ey) + 1e-300;
                nx = ey / el;
                ny = -ex / el;
            }
        }
        if ((ay > py) != (by > py) && px < ax + (py - ay) * ex / (by - ay)) inside = !inside;
    }
    const double d = std::sqrt(best);
    if (inside) {
        nx = -nx;
        ny = -ny;
        return -d;
    }
    return d;
}
}  // namespace

// Ray o + t d (|d| = 1) against a sphere: entry and exit.
bool raySphere(const Vec3& o, const Vec3& d, const Vec3& c, double r, double& t0, double& t1) {
    const Vec3 oc = o - c;
    const double b = dot(oc, d), cc = length2(oc) - r * r, disc = b * b - cc;
    if (disc < 0) return false;
    const double s = std::sqrt(disc);
    t0 = -b - s;
    t1 = -b + s;
    return true;
}

// Ray against a capsule (convex: the hull of the cylinder and end spheres).
bool rayCapsule(const Vec3& o, const Vec3& d, const Vec3& a, const Vec3& b, double r, double& tin, double& tout) {
    double lo = 1e30, hi = -1e30, s0, s1;
    if (raySphere(o, d, a, r, s0, s1)) {
        lo = std::min(lo, s0);
        hi = std::max(hi, s1);
    }
    if (raySphere(o, d, b, r, s0, s1)) {
        lo = std::min(lo, s0);
        hi = std::max(hi, s1);
    }
    const Vec3 ax = b - a;
    const double L = length(ax);
    if (L > 1e-9) {
        const Vec3 u = ax / L;
        const Vec3 oa = o - a;
        const Vec3 dp = d - u * dot(d, u), op = oa - u * dot(oa, u);
        const double A = length2(dp), B = 2 * dot(dp, op), C = length2(op) - r * r;
        double c0 = -1e30, c1 = 1e30;
        bool ok = true;
        if (A > 1e-12) {
            const double disc = B * B - 4 * A * C;
            if (disc < 0) ok = false;
            else {
                const double s = std::sqrt(disc);
                c0 = (-B - s) / (2 * A);
                c1 = (-B + s) / (2 * A);
            }
        } else if (C > 0) ok = false;
        if (ok) {
            // Slab between the end planes.
            const double du = dot(d, u), ou = dot(oa, u);
            double p0 = -1e30, p1 = 1e30;
            if (std::abs(du) > 1e-12) {
                p0 = (0 - ou) / du;
                p1 = (L - ou) / du;
                if (p0 > p1) std::swap(p0, p1);
            } else if (ou < 0 || ou > L) ok = false;
            const double e0 = std::max(c0, p0), e1 = std::min(c1, p1);
            if (ok && e0 <= e1) {
                lo = std::min(lo, e0);
                hi = std::max(hi, e1);
            }
        }
    }
    if (lo > hi) return false;
    tin = lo;
    tout = hi;
    return true;
}

// Ray against an oriented box (slabs in the box frame).
bool rayBox(const Collider& c, const Vec3& o, const Vec3& d, double& tin, double& tout) {
    const Vec3 lo = c.toLocal(o), ld = c.rot.unrotate(d);
    double t0 = -1e30, t1 = 1e30;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(ld[i]) < 1e-12) {
            if (std::abs(lo[i]) > c.half[i]) return false;
            continue;
        }
        double a = (-c.half[i] - lo[i]) / ld[i], b = (c.half[i] - lo[i]) / ld[i];
        if (a > b) std::swap(a, b);
        t0 = std::max(t0, a);
        t1 = std::min(t1, b);
    }
    if (t0 > t1) return false;
    tin = t0;
    tout = t1;
    return true;
}

bool Collider::inHole(const Vec3& l) const {
    if (map && map->open(l.x, l.y)) return true;
    for (const auto& h : holes)
        if (sq(l.x - h.x) + sq(l.y - h.y) < h.radius * h.radius) return true;
    return false;
}

double Collider::sdf(const Vec3& p, Vec3* normal) const {
    switch (kind) {
        case ShapeKind::Plane: {
            const Vec3 n = rot.rotate({0, 1, 0});
            if (normal) *normal = n;
            return dot(p - pos, n);
        }
        case ShapeKind::Box: {
            const Vec3 l = toLocal(p);
            const Vec3 q{std::abs(l.x) - half.x, std::abs(l.y) - half.y, std::abs(l.z) - half.z};
            double d;
            Vec3 nl;
            if (q.x > 0 || q.y > 0 || q.z > 0) {
                const Vec3 o{std::max(q.x, 0.0), std::max(q.y, 0.0), std::max(q.z, 0.0)};
                d = length(o);
                nl = Vec3{l.x < 0 ? -o.x : o.x, l.y < 0 ? -o.y : o.y, l.z < 0 ? -o.z : o.z} / std::max(d, 1e-300);
            } else {
                // Inside: the nearest face.
                if (q.x >= q.y && q.x >= q.z) {
                    d = q.x;
                    nl = {l.x < 0 ? -1.0 : 1.0, 0, 0};
                } else if (q.y >= q.z) {
                    d = q.y;
                    nl = {0, l.y < 0 ? -1.0 : 1.0, 0};
                } else {
                    d = q.z;
                    nl = {0, 0, l.z < 0 ? -1.0 : 1.0};
                }
            }
            if (map && map->open(l.x, l.y)) {
                // Where the pieces have gone there is no glass at all.
                if (normal) *normal = rot.rotate({0, 0, l.z < 0 ? -1.0 : 1.0});
                return std::max(std::abs(l.z) - half.z, 0.0) + 1e3;
            }
            if (!holes.empty() && inHole(l)) {
                // Inside a hole's bore: the distance to the bore's wall.
                double best = 1e30;
                Vec3 bn;
                for (const auto& h : holes) {
                    const double dx = l.x - h.x, dy = l.y - h.y, rr = std::sqrt(dx * dx + dy * dy);
                    if (rr >= h.radius) continue;
                    const double dd = h.radius - rr;
                    if (dd < best) {
                        best = dd;
                        bn = rr > 1e-12 ? Vec3{-dx / rr, -dy / rr, 0} : Vec3{1, 0, 0};
                    }
                }
                // Out beyond the faces the plain box distance is larger and wins.
                const double faceGap = std::abs(l.z) - half.z;
                if (faceGap > best) {
                    if (normal) *normal = rot.rotate({0, 0, l.z < 0 ? -1.0 : 1.0});
                    return faceGap;
                }
                if (normal) *normal = rot.rotate(bn);
                return best;
            }
            if (normal) *normal = rot.rotate(nl);
            return d;
        }
        case ShapeKind::Capsule: {
            const Vec3 ab = b - a;
            const double f = clamp(dot(p - a, ab) / std::max(length2(ab), 1e-300), 0.0, 1.0);
            const Vec3 c = a + ab * f, dpc = p - c;
            const double l = length(dpc);
            if (normal) *normal = l > 1e-12 ? dpc / l : perpendicular(normalize(ab));
            return l - radius;
        }
        case ShapeKind::Slab: {
            const Vec3 l = toLocal(p);
            double ex, ey;
            double dp = polygonDistance(outline, l.x, l.y, ex, ey);
            for (const auto& h : outlineHoles) {
                // Inside a hole is outside the glass.
                double hx, hy;
                const double dh = -polygonDistance(h, l.x, l.y, hx, hy);
                if (dh > dp) {
                    dp = dh;
                    ex = -hx;
                    ey = -hy;
                }
            }
            const double dz = std::abs(l.z) - half.z;
            const double sz = l.z < 0 ? -1.0 : 1.0;
            Vec3 nl;
            double d;
            if (dp > 0 && dz > 0) {
                d = std::sqrt(dp * dp + dz * dz);
                nl = Vec3{ex * dp, ey * dp, sz * dz} / d;
            } else if (dp > dz) {
                d = dp;
                nl = {ex, ey, 0};
            } else {
                d = dz;
                nl = {0, 0, sz};
            }
            if (normal) *normal = rot.rotate(nl);
            return d;
        }
        case ShapeKind::Shell: {
            const Vec3 q = toVessel(p);
            WallPoint w;
            Vec3 n;
            const double d = shell->sdf(q, &n, &w);
            if (shards && shard >= 0 && shards->shardAt(*shell, w.s, w.theta) != shard) {
                if (normal) *normal = rot.rotate(n);
                return 1e3;  // that part of the wall belongs to another shard
            }
            if (normal) *normal = rot.rotate(n);
            return d;
        }
    }
    return 1e30;
}

}  // namespace lqd
