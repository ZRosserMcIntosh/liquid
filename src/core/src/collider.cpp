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

bool Collider::inHole(const Vec3& l) const {
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
