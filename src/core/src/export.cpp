// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/export.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

namespace lqd {

// ---- Recording -------------------------------------------------------------------

void Recorder::captureNow(const World& w) {
    if (frames.size() >= max_) return;
    RecordedFrame f;
    f.time = w.time();
    f.particles.reserve(w.fluid.size() * 3);
    for (const auto& p : w.fluid.x) f.particles.insert(f.particles.end(), {float(p.x), float(p.y), float(p.z)});
    for (const auto& d : w.spray.drops) f.drops.insert(f.drops.end(), {float(d.x.x), float(d.x.y), float(d.x.z), float(d.d)});
    for (const auto& b : w.bodies)
        f.poses.insert(f.poses.end(), {float(b.x.x), float(b.x.y), float(b.x.z), float(b.q.x), float(b.q.y), float(b.q.z), float(b.q.w)});
    if (w.ragdoll.active())
        for (int k = 0; k < w.ragdoll.segments(); ++k) {
            const Vec3 a = w.ragdoll.a(k), b = w.ragdoll.b(k);
            f.dummy.insert(f.dummy.end(), {float(a.x), float(a.y), float(a.z), float(b.x), float(b.y), float(b.z)});
        }
    for (const auto& bp : w.brokenPanes()) {
        std::vector<float> v;
        v.reserve(bp.pieces.size() * 8);
        for (const auto& p : bp.pieces)
            v.insert(v.end(), {float(p.pos.x), float(p.pos.y), float(p.pos.z), float(p.rot.x), float(p.rot.y),
                               float(p.rot.z), float(p.rot.w), p.posed && p.state != -1 ? 1.0f : 0.0f});
        f.pieces.push_back(std::move(v));
    }
    for (const auto& im : w.impactors) {
        const Vec3 p = w.colliders[im.collider].pos;
        f.impactors.insert(f.impactors.end(), {float(p.x), float(p.y), float(p.z)});
    }
    f.stains = w.spray.stains.size();
    spacing = w.fluid.settings().spacing;
    frames.push_back(std::move(f));
}

bool Recorder::capture(const World& w) {
    if (next_ < 0) next_ = w.time();
    if (w.time() + 1e-12 < next_) return false;
    captureNow(w);
    next_ += 1.0 / fps_;
    if (next_ < w.time()) next_ = w.time() + 1.0 / fps_;
    return true;
}

VesselMesh bodyMesh(const World& w, int body) {
    VesselMesh m;
    if (body < 0 || body >= static_cast<int>(w.bodies.size())) return m;
    const RigidBody& b = w.bodies[body];
    if (!b.geom) return m;
    if (b.shard >= 0) {
        const ShardPattern* P = w.shardsOf(body);
        if (!P) return m;
        const int k = b.shard;
        m = b.geom->mesh([&](int c) { return P->cellShard[c] == k; });
    } else {
        m = b.geom->mesh();
    }
    for (auto& p : m.pos) p -= b.comLocal;
    return m;
}

// ---- glTF writer --------------------------------------------------------------------

namespace {

std::string num(double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.9g", std::isfinite(v) ? v : 0.0);  // round-trips a float32
    return b;
}

struct Prim {
    int pos = -1, nrm = -1, idx = -1, material = -1;
};

class Glb {
public:
    int view(const void* data, std::size_t bytes, int target) {
        while (bin.size() % 4) bin.push_back(0);
        const std::size_t off = bin.size();
        const auto* p = static_cast<const std::uint8_t*>(data);
        bin.insert(bin.end(), p, p + bytes);
        std::string v = "{\"buffer\":0,\"byteOffset\":" + std::to_string(off) + ",\"byteLength\":" + std::to_string(bytes);
        if (target) v += ",\"target\":" + std::to_string(target);
        v += "}";
        views.push_back(v);
        return static_cast<int>(views.size()) - 1;
    }
    int floats(const std::vector<float>& d, int comps, int target, bool bounds) {
        const int bv = view(d.data(), d.size() * 4, target);
        static const char* type[] = {"", "SCALAR", "VEC2", "VEC3", "VEC4"};
        const std::size_t count = d.size() / comps;
        std::string a = "{\"bufferView\":" + std::to_string(bv) + ",\"componentType\":5126,\"count\":" +
                        std::to_string(count) + ",\"type\":\"" + type[comps] + "\"";
        if (bounds && count > 0) {
            std::vector<float> lo(comps, 1e30f), hi(comps, -1e30f);
            for (std::size_t i = 0; i < count; ++i)
                for (int c = 0; c < comps; ++c) {
                    lo[c] = std::min(lo[c], d[i * comps + c]);
                    hi[c] = std::max(hi[c], d[i * comps + c]);
                }
            a += ",\"min\":[";
            for (int c = 0; c < comps; ++c) a += (c ? "," : "") + num(lo[c]);
            a += "],\"max\":[";
            for (int c = 0; c < comps; ++c) a += (c ? "," : "") + num(hi[c]);
            a += "]";
        }
        a += "}";
        accessors.push_back(a);
        return static_cast<int>(accessors.size()) - 1;
    }
    int indices(const std::vector<std::uint32_t>& d) {
        const int bv = view(d.data(), d.size() * 4, 34963);
        accessors.push_back("{\"bufferView\":" + std::to_string(bv) + ",\"componentType\":5125,\"count\":" +
                            std::to_string(d.size()) + ",\"type\":\"SCALAR\"}");
        return static_cast<int>(accessors.size()) - 1;
    }
    int material(const std::string& name, const Vec3& c, double alpha, double rough, double metal) {
        std::string m = "{\"name\":\"" + name + "\",\"pbrMetallicRoughness\":{\"baseColorFactor\":[" + num(c.x) + "," +
                        num(c.y) + "," + num(c.z) + "," + num(alpha) + "],\"metallicFactor\":" + num(metal) +
                        ",\"roughnessFactor\":" + num(rough) + "},\"doubleSided\":true";
        if (alpha < 0.999) m += ",\"alphaMode\":\"BLEND\"";
        m += "}";
        materials.push_back(m);
        return static_cast<int>(materials.size()) - 1;
    }
    int mesh(const std::string& name, const std::vector<Prim>& prims) {
        std::string m = "{\"name\":\"" + name + "\",\"primitives\":[";
        for (std::size_t i = 0; i < prims.size(); ++i) {
            const Prim& p = prims[i];
            m += (i ? "," : "") + std::string("{\"attributes\":{\"POSITION\":") + std::to_string(p.pos);
            if (p.nrm >= 0) m += ",\"NORMAL\":" + std::to_string(p.nrm);
            m += "},\"indices\":" + std::to_string(p.idx) + ",\"material\":" + std::to_string(p.material) + "}";
        }
        m += "]}";
        meshes.push_back(m);
        return static_cast<int>(meshes.size()) - 1;
    }
    Prim prim(const std::vector<float>& pos, const std::vector<float>& nrm, const std::vector<std::uint32_t>& idx,
              int material) {
        Prim p;
        p.pos = floats(pos, 3, 34962, true);
        p.nrm = nrm.empty() ? -1 : floats(nrm, 3, 34962, false);
        p.idx = indices(idx);
        p.material = material;
        return p;
    }
    int node(const std::string& name, int mesh, const Vec3& t = {}, const Quat& r = {}, double scale = 1,
             const std::vector<int>& children = {}) {
        std::string n = "{\"name\":\"" + name + "\"";
        if (mesh >= 0) n += ",\"mesh\":" + std::to_string(mesh);
        n += ",\"translation\":[" + num(t.x) + "," + num(t.y) + "," + num(t.z) + "]";
        n += ",\"rotation\":[" + num(r.x) + "," + num(r.y) + "," + num(r.z) + "," + num(r.w) + "]";
        n += ",\"scale\":[" + num(scale) + "," + num(scale) + "," + num(scale) + "]";
        if (!children.empty()) {
            n += ",\"children\":[";
            for (std::size_t i = 0; i < children.size(); ++i) n += (i ? "," : "") + std::to_string(children[i]);
            n += "]";
        }
        n += "}";
        nodes.push_back(n);
        return static_cast<int>(nodes.size()) - 1;
    }
    void channel(int node, const char* path, int input, int output, const char* interp) {
        samplers.push_back("{\"input\":" + std::to_string(input) + ",\"output\":" + std::to_string(output) +
                           ",\"interpolation\":\"" + interp + "\"}");
        channels.push_back("{\"sampler\":" + std::to_string(samplers.size() - 1) + ",\"target\":{\"node\":" +
                           std::to_string(node) + ",\"path\":\"" + path + "\"}}");
    }
    std::vector<std::uint8_t> finish(const std::vector<int>& roots) {
        std::string j = "{\"asset\":{\"version\":\"2.0\",\"generator\":\"Liquid Dynamics\"},\"scene\":0,\"scenes\":[{\"nodes\":[";
        for (std::size_t i = 0; i < roots.size(); ++i) j += (i ? "," : "") + std::to_string(roots[i]);
        j += "]}]";
        auto list = [&](const char* key, const std::vector<std::string>& v) {
            if (v.empty()) return;
            j += ",\"" + std::string(key) + "\":[";
            for (std::size_t i = 0; i < v.size(); ++i) j += (i ? "," : "") + v[i];
            j += "]";
        };
        list("nodes", nodes);
        list("meshes", meshes);
        list("materials", materials);
        list("accessors", accessors);
        list("bufferViews", views);
        if (!channels.empty()) {
            j += ",\"animations\":[{\"name\":\"simulation\",\"channels\":[";
            for (std::size_t i = 0; i < channels.size(); ++i) j += (i ? "," : "") + channels[i];
            j += "],\"samplers\":[";
            for (std::size_t i = 0; i < samplers.size(); ++i) j += (i ? "," : "") + samplers[i];
            j += "]}]";
        }
        while (bin.size() % 4) bin.push_back(0);
        j += ",\"buffers\":[{\"byteLength\":" + std::to_string(bin.size()) + "}]}";
        while (j.size() % 4) j += ' ';
        std::vector<std::uint8_t> out;
        auto u32 = [&](std::uint32_t v) {
            for (int i = 0; i < 4; ++i) out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
        };
        u32(0x46546C67);  // "glTF"
        u32(2);
        u32(static_cast<std::uint32_t>(12 + 8 + j.size() + 8 + bin.size()));
        u32(static_cast<std::uint32_t>(j.size()));
        u32(0x4E4F534A);  // "JSON"
        out.insert(out.end(), j.begin(), j.end());
        u32(static_cast<std::uint32_t>(bin.size()));
        u32(0x004E4942);  // "BIN\0"
        out.insert(out.end(), bin.begin(), bin.end());
        return out;
    }

    std::vector<std::uint8_t> bin;
    std::vector<std::string> views, accessors, materials, meshes, nodes, samplers, channels;
};

struct Geo {
    std::vector<float> pos, nrm;
    std::vector<std::uint32_t> idx;
    void vert(const Vec3& p, const Vec3& n) {
        pos.insert(pos.end(), {float(p.x), float(p.y), float(p.z)});
        nrm.insert(nrm.end(), {float(n.x), float(n.y), float(n.z)});
    }
    std::uint32_t base() const { return static_cast<std::uint32_t>(pos.size() / 3); }
};

void addQuad(Geo& g, const Vec3& c, const Vec3& u, const Vec3& v, const Vec3& n) {
    const std::uint32_t b = g.base();
    g.vert(c - u - v, n);
    g.vert(c + u - v, n);
    g.vert(c + u + v, n);
    g.vert(c - u + v, n);
    g.idx.insert(g.idx.end(), {b, b + 1, b + 2, b, b + 2, b + 3});
}

void addBox(Geo& g, const Collider& c) {
    for (int ax = 0; ax < 3; ++ax)
        for (int s = -1; s <= 1; s += 2) {
            Vec3 n, u, v;
            n[ax] = s;
            u[(ax + 1) % 3] = c.half[(ax + 1) % 3];
            v[(ax + 2) % 3] = c.half[(ax + 2) % 3];
            Vec3 ctr;
            ctr[ax] = s * c.half[ax];
            addQuad(g, c.pos + c.rot.rotate(ctr), c.rot.rotate(u), c.rot.rotate(v), c.rot.rotate(n));
        }
}

void addCapsule(Geo& g, const Vec3& a, const Vec3& b, double r) {
    const Vec3 ax = normalize(b - a), u = perpendicular(ax), w = cross(ax, u);
    const int seg = 14, ring = 5;
    // Rings from the sphere at a, along the cylinder, to the sphere at b.
    std::vector<std::pair<double, double>> prof;  // (axial offset, radius factor) with the centre they belong to
    std::vector<int> end;
    for (int i = 0; i <= ring; ++i) {
        const double t = -kPi / 2 + kPi / 2 * i / ring;
        prof.push_back({std::sin(t), std::cos(t)});
        end.push_back(0);
    }
    for (int i = 0; i <= ring; ++i) {
        const double t = kPi / 2 * i / ring;
        prof.push_back({std::sin(t), std::cos(t)});
        end.push_back(1);
    }
    const std::uint32_t b0 = g.base();
    for (std::size_t i = 0; i < prof.size(); ++i) {
        const Vec3 cen = end[i] ? b : a;
        for (int j = 0; j <= seg; ++j) {
            const double ph = 2 * kPi * j / seg;
            const Vec3 n = ax * prof[i].first + (u * std::cos(ph) + w * std::sin(ph)) * prof[i].second;
            g.vert(cen + n * r, n);
        }
    }
    for (std::size_t i = 0; i + 1 < prof.size(); ++i)
        for (int j = 0; j < seg; ++j) {
            const std::uint32_t p = b0 + static_cast<std::uint32_t>(i * (seg + 1) + j), q = p + seg + 1;
            g.idx.insert(g.idx.end(), {p, q, p + 1, p + 1, q, q + 1});
        }
}

void addDrop(Geo& g, const Vec3& c, double r) {
    static const Vec3 dirs[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    const std::uint32_t b = g.base();
    for (const auto& d : dirs) g.vert(c + d * r, d);
    static const std::uint32_t f[8][3] = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4}, {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
    for (const auto& t : f) g.idx.insert(g.idx.end(), {b + t[0], b + t[1], b + t[2]});
}

// Triangles of a simple polygon with holes (x, y pairs): the holes are bridged
// into the outline, then ears are clipped. Indices into the merged loop `pts`.
void triangulate(const std::vector<double>& outer, const std::vector<std::vector<double>>& holes,
                 std::vector<std::pair<double, double>>& pts, std::vector<std::uint32_t>& tri) {
    auto loop = [](const std::vector<double>& v, bool ccw) {
        std::vector<std::pair<double, double>> l;
        for (std::size_t i = 0; i + 1 < v.size(); i += 2) l.push_back({v[i], v[i + 1]});
        double a = 0;
        for (std::size_t i = 0, j = l.size() - 1; i < l.size(); j = i++) a += l[j].first * l[i].second - l[i].first * l[j].second;
        if ((a > 0) != ccw) std::reverse(l.begin(), l.end());
        return l;
    };
    auto cross2 = [](std::pair<double, double> o, std::pair<double, double> a, std::pair<double, double> b) {
        return (a.first - o.first) * (b.second - o.second) - (a.second - o.second) * (b.first - o.first);
    };
    auto crosses = [&](std::pair<double, double> a, std::pair<double, double> b, std::pair<double, double> c,
                       std::pair<double, double> d) {
        const double d1 = cross2(a, b, c), d2 = cross2(a, b, d), d3 = cross2(c, d, a), d4 = cross2(c, d, b);
        return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) && std::abs(d1) > 1e-14 && std::abs(d2) > 1e-14;
    };
    pts = loop(outer, true);
    if (pts.size() < 3) return;
    std::vector<std::vector<std::pair<double, double>>> hs;
    for (const auto& h : holes)
        if (h.size() >= 6) hs.push_back(loop(h, false));
    std::sort(hs.begin(), hs.end(), [](const auto& a, const auto& b) {
        double ma = -1e30, mb = -1e30;
        for (const auto& p : a) ma = std::max(ma, p.first);
        for (const auto& p : b) mb = std::max(mb, p.first);
        return ma > mb;
    });
    for (const auto& h : hs) {
        std::size_t hi = 0;
        for (std::size_t i = 1; i < h.size(); ++i)
            if (h[i].first > h[hi].first) hi = i;
        // The nearest outline vertex the bridge to it crosses nothing on the way.
        std::size_t best = 0;
        double bd = 1e30;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            const double d = std::hypot(pts[i].first - h[hi].first, pts[i].second - h[hi].second);
            if (d >= bd) continue;
            bool clear = true;
            for (std::size_t j = 0, k = pts.size() - 1; j < pts.size() && clear; k = j++)
                if (j != i && k != i && crosses(h[hi], pts[i], pts[k], pts[j])) clear = false;
            for (std::size_t j = 0, k = h.size() - 1; j < h.size() && clear; k = j++)
                if (j != hi && k != hi && crosses(h[hi], pts[i], h[k], h[j])) clear = false;
            if (clear) {
                bd = d;
                best = i;
            }
        }
        std::vector<std::pair<double, double>> merged(pts.begin(), pts.begin() + best + 1);
        for (std::size_t k = 0; k <= h.size(); ++k) merged.push_back(h[(hi + k) % h.size()]);
        merged.insert(merged.end(), pts.begin() + best, pts.end());
        pts.swap(merged);
    }
    std::vector<std::uint32_t> idx(pts.size());
    for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = static_cast<std::uint32_t>(i);
    std::size_t guard = 0;
    while (idx.size() > 3 && guard++ < 4 * pts.size() * pts.size()) {
        bool cut = false;
        for (std::size_t i = 0; i < idx.size() && !cut; ++i) {
            const std::uint32_t a = idx[(i + idx.size() - 1) % idx.size()], b = idx[i], c = idx[(i + 1) % idx.size()];
            if (cross2(pts[a], pts[b], pts[c]) <= 1e-16) continue;  // reflex or flat
            bool empty = true;
            for (std::uint32_t q : idx) {
                if (q == a || q == b || q == c || pts[q] == pts[a] || pts[q] == pts[b] || pts[q] == pts[c]) continue;
                if (cross2(pts[a], pts[b], pts[q]) >= 0 && cross2(pts[b], pts[c], pts[q]) >= 0 && cross2(pts[c], pts[a], pts[q]) >= 0) {
                    empty = false;
                    break;
                }
            }
            if (!empty) continue;
            tri.insert(tri.end(), {a, b, c});
            idx.erase(idx.begin() + static_cast<std::ptrdiff_t>(i));
            cut = true;
        }
        if (!cut) {
            // Degenerate leftovers: fan them rather than loop for ever.
            for (std::size_t i = 1; i + 1 < idx.size(); ++i) tri.insert(tri.end(), {idx[0], idx[i], idx[i + 1]});
            idx.clear();
        }
    }
    if (idx.size() == 3) tri.insert(tri.end(), {idx[0], idx[1], idx[2]});
}

// A flat piece of glass: its outline extruded +-h along z, in its own frame.
void addSlab(Geo& g, const std::vector<double>& outline, const std::vector<std::vector<double>>& holes, double h) {
    std::vector<std::pair<double, double>> pts;
    std::vector<std::uint32_t> tri;
    triangulate(outline, holes, pts, tri);
    if (tri.empty()) return;
    for (int side = -1; side <= 1; side += 2) {
        const std::uint32_t b = g.base();
        for (const auto& p : pts) g.vert({p.first, p.second, side * h}, {0, 0, double(side)});
        for (std::size_t i = 0; i + 2 < tri.size(); i += 3)
            if (side > 0) g.idx.insert(g.idx.end(), {b + tri[i], b + tri[i + 1], b + tri[i + 2]});
            else g.idx.insert(g.idx.end(), {b + tri[i], b + tri[i + 2], b + tri[i + 1]});
    }
    auto rim = [&](const std::vector<double>& l) {
        const std::size_t n = l.size() / 2;
        for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
            const Vec3 a{l[2 * j], l[2 * j + 1], 0}, c{l[2 * i], l[2 * i + 1], 0};
            const Vec3 e = c - a, nn = normalize(Vec3{e.y, -e.x, 0}, {1, 0, 0});
            const std::uint32_t b = g.base();
            g.vert(a + Vec3{0, 0, -h}, nn);
            g.vert(c + Vec3{0, 0, -h}, nn);
            g.vert(c + Vec3{0, 0, h}, nn);
            g.vert(a + Vec3{0, 0, h}, nn);
            g.idx.insert(g.idx.end(), {b, b + 1, b + 2, b, b + 2, b + 3});
        }
    };
    rim(outline);
    for (const auto& hl : holes) rim(hl);
}

void addStain(Geo& g, const Stain& s) {
    const Vec3 n = s.normal, u = normalize(s.dir - n * dot(s.dir, n), perpendicular(n)), v = cross(n, u);
    const Vec3 c = s.pos + n * 4e-4;
    const int seg = 16;
    const std::uint32_t b = g.base();
    g.vert(c, n);
    for (int k = 0; k < seg; ++k) {
        const double a = 2 * kPi * k / seg;
        g.vert(c + u * (0.5 * s.length * std::cos(a)) + v * (0.5 * s.width * std::sin(a)), n);
    }
    for (int k = 0; k < seg; ++k)
        g.idx.insert(g.idx.end(), {b, b + 1 + static_cast<std::uint32_t>(k), b + 1 + static_cast<std::uint32_t>((k + 1) % seg)});
}

}  // namespace

std::vector<std::uint8_t> exportGlb(const World& w, const Recorder& r, const ExportOptions& o) {
    Glb G;
    std::vector<int> roots;
    const LiquidDef* L = w.fluid.liquid() ? w.fluid.liquid() : w.spray.liquid();
    const Vec3 liqColor = L ? L->color : Vec3{0.4, 0.02, 0.08};
    const int mLiquid = G.material("liquid", liqColor, L ? std::max(0.5, L->opacity) : 0.9, 0.08, 0);
    const int mStain = G.material("stain", L ? L->stainColor : Vec3{0.3, 0.02, 0.05}, 1, 0.5, 0);
    const int mGlass = G.material("glass", {0.85, 0.92, 0.95}, 0.25, 0.04, 0);
    std::vector<const RecordedFrame*> frames;
    for (std::size_t i = 0; i < r.frames.size(); i += std::max(1, o.stride)) frames.push_back(&r.frames[i]);
    const std::size_t F = frames.size();
    const double t0 = F ? frames[0]->time : 0;
    std::vector<float> times;
    for (const auto* f : frames) times.push_back(static_cast<float>(f->time - t0));
    const int timeAcc = F ? G.floats(times, 1, 0, true) : -1;

    // Liquid flipbook.
    std::vector<int> liquidNodes;
    for (std::size_t k = 0; k < F; ++k) {
        const RecordedFrame& f = *frames[k];
        std::vector<Vec3> pts;
        for (std::size_t i = 0; i + 2 < f.particles.size(); i += 3) pts.push_back({f.particles[i], f.particles[i + 1], f.particles[i + 2]});
        std::vector<Prim> prims;
        MesherSettings ms;
        ms.cell = 0.75 * r.spacing;  // a little coarser than on screen: one mesh per frame adds up
        const SurfaceMesh sm = meshLiquid(pts, r.spacing, ms);
        if (sm.triangles() > 0) prims.push_back(G.prim(sm.pos, sm.nrm, sm.tri, mLiquid));
        if (o.drops && !f.drops.empty()) {
            Geo g;
            for (std::size_t i = 0; i + 3 < f.drops.size(); i += 4) addDrop(g, {f.drops[i], f.drops[i + 1], f.drops[i + 2]}, 0.5 * f.drops[i + 3]);
            prims.push_back(G.prim(g.pos, g.nrm, g.idx, mLiquid));
        }
        const int mesh = prims.empty() ? -1 : G.mesh("liquid_" + std::to_string(k), prims);
        const int nd = G.node("liquid_" + std::to_string(k), mesh, {}, {}, k == 0 ? 1 : 0);
        liquidNodes.push_back(nd);
        std::vector<float> sc;
        for (std::size_t q = 0; q < F; ++q) {
            const float s = q == k ? 1.0f : 0.0f;
            sc.insert(sc.end(), {s, s, s});
        }
        if (F > 1) G.channel(nd, "scale", timeAcc, G.floats(sc, 3, 0, false), "STEP");
    }
    roots.push_back(G.node("liquid", -1, {}, {}, 1, liquidNodes));

    // Glass bodies and shards.
    std::vector<int> glassNodes;
    for (int b = 0; b < static_cast<int>(w.bodies.size()); ++b) {
        const RigidBody& B = w.bodies[b];
        const VesselMesh vm = bodyMesh(w, b);
        if (vm.tri.empty()) continue;
        Geo g;
        for (std::size_t i = 0; i < vm.pos.size(); ++i) g.vert(vm.pos[i], vm.nrm[i]);
        g.idx = vm.tri;
        const int mesh = G.mesh(B.name, {G.prim(g.pos, g.nrm, g.idx, mGlass)});
        const int nd = G.node(B.name, mesh, B.x, B.q, 1);
        glassNodes.push_back(nd);
        if (F == 0) continue;
        // The first recorded pose stands in for the frames before it existed.
        std::size_t first = 0;
        while (first < F && frames[first]->poses.size() < static_cast<std::size_t>(7 * (b + 1))) ++first;
        std::vector<float> tr, ro, sc;
        for (std::size_t k = 0; k < F; ++k) {
            const RecordedFrame& f = *frames[std::max(k, std::min(first, F - 1))];
            const float* p = f.poses.size() >= static_cast<std::size_t>(7 * (b + 1)) ? &f.poses[7 * b] : nullptr;
            if (p) {
                tr.insert(tr.end(), {p[0], p[1], p[2]});
                ro.insert(ro.end(), {p[3], p[4], p[5], p[6]});
            } else {
                tr.insert(tr.end(), {float(B.x.x), float(B.x.y), float(B.x.z)});
                ro.insert(ro.end(), {float(B.q.x), float(B.q.y), float(B.q.z), float(B.q.w)});
            }
            const double t = frames[k]->time;
            const bool alive = t + 1e-9 >= B.born && (B.brokenAt < 0 || t < B.brokenAt);
            const float s = alive ? 1.0f : 0.0f;
            sc.insert(sc.end(), {s, s, s});
        }
        G.channel(nd, "translation", timeAcc, G.floats(tr, 3, 0, false), "LINEAR");
        G.channel(nd, "rotation", timeAcc, G.floats(ro, 4, 0, false), "LINEAR");
        G.channel(nd, "scale", timeAcc, G.floats(sc, 3, 0, false), "STEP");
    }
    if (!glassNodes.empty()) roots.push_back(G.node("glass", -1, {}, {}, 1, glassNodes));

    // The frame a moment falls in (the stains and pieces that appear then).
    auto frameOf = [&](double t) {
        std::size_t k = 0;
        while (k + 1 < F && frames[k + 1]->time <= t) ++k;
        return k;
    };
    auto appear = [&](int nd, std::size_t k) {
        if (k == 0 || F < 2) return;
        std::vector<float> sc;
        for (std::size_t q = 0; q < F; ++q) {
            const float v = q >= k ? 1.0f : 0.0f;
            sc.insert(sc.end(), {v, v, v});
        }
        G.channel(nd, "scale", timeAcc, G.floats(sc, 3, 0, false), "STEP");
    };
    // The pieces of broken panes, each keyed every frame, with the stains on it
    // as its children (they fall with it).
    std::vector<unsigned char> onPiece(w.spray.stains.size(), 0);
    for (const auto& ps : w.pieceStains())
        if (ps.piece >= 0) onPiece[ps.stain] = 1;
    for (std::size_t bi = 0; bi < w.brokenPanes().size(); ++bi) {
        const PanePieces& bp = w.brokenPanes()[bi];
        std::vector<int> pieceNodes;
        for (int k = 0; k < static_cast<int>(bp.pieces.size()); ++k) {
            const PanePiece& pc = bp.pieces[k];
            if (!pc.posed || pc.outline.size() < 6) continue;
            Geo g;
            addSlab(g, pc.outline, pc.holes, 0.5 * pc.thickness);
            if (g.idx.empty()) continue;
            std::vector<int> kids;
            if (o.stains) {
                std::vector<std::pair<std::size_t, Geo>> byFrame;
                for (const auto& ps : w.pieceStains()) {
                    if (ps.pane != bp.pane || ps.piece != k) continue;
                    Stain st = w.spray.stains[ps.stain];
                    if (st.gone) continue;
                    st.pos = ps.p;
                    st.normal = ps.n;
                    st.dir = ps.d;
                    const std::size_t fk = frameOf(st.time);
                    auto it = std::find_if(byFrame.begin(), byFrame.end(), [&](const auto& e) { return e.first == fk; });
                    if (it == byFrame.end()) {
                        byFrame.push_back({fk, Geo{}});
                        it = byFrame.end() - 1;
                    }
                    addStain(it->second, st);
                }
                for (auto& [fk, sg] : byFrame) {
                    const std::string nm = "piece_" + std::to_string(k) + "_stains_" + std::to_string(fk);
                    const int nd = G.node(nm, G.mesh(nm, {G.prim(sg.pos, sg.nrm, sg.idx, mStain)}), {}, {}, fk == 0 ? 1 : 0);
                    appear(nd, fk);
                    kids.push_back(nd);
                }
            }
            const std::string nm = w.colliders[bp.pane].name + " piece " + std::to_string(k);
            const int nd = G.node(nm, G.mesh(nm, {G.prim(g.pos, g.nrm, g.idx, mGlass)}), pc.pos, pc.rot, 1, kids);
            pieceNodes.push_back(nd);
            if (F == 0) continue;
            std::vector<float> tr, ro, sc;
            // Before it was recorded it stands where it was first seen.
            const float* firstSeen = nullptr;
            for (std::size_t q = 0; q < F && !firstSeen; ++q) {
                const auto& fp = frames[q]->pieces;
                if (bi < fp.size() && fp[bi].size() >= static_cast<std::size_t>(8 * (k + 1)) && fp[bi][8 * k + 7] > 0) firstSeen = &fp[bi][8 * k];
            }
            const float* last = firstSeen;
            for (std::size_t q = 0; q < F; ++q) {
                const auto& fp = frames[q]->pieces;
                const float* p = bi < fp.size() && fp[bi].size() >= static_cast<std::size_t>(8 * (k + 1)) ? &fp[bi][8 * k] : nullptr;
                const bool shown = p && p[7] > 0;
                if (shown) last = p;
                const float* use = last;
                if (use) {
                    tr.insert(tr.end(), {use[0], use[1], use[2]});
                    ro.insert(ro.end(), {use[3], use[4], use[5], use[6]});
                } else {
                    tr.insert(tr.end(), {float(pc.pos.x), float(pc.pos.y), float(pc.pos.z)});
                    ro.insert(ro.end(), {float(pc.rot.x), float(pc.rot.y), float(pc.rot.z), float(pc.rot.w)});
                }
                const float s = shown ? 1.0f : 0.0f;
                sc.insert(sc.end(), {s, s, s});
            }
            G.channel(nd, "translation", timeAcc, G.floats(tr, 3, 0, false), "LINEAR");
            G.channel(nd, "rotation", timeAcc, G.floats(ro, 4, 0, false), "LINEAR");
            G.channel(nd, "scale", timeAcc, G.floats(sc, 3, 0, false), "STEP");
        }
        if (!pieceNodes.empty()) roots.push_back(G.node(w.colliders[bp.pane].name + " pieces", -1, {}, {}, 1, pieceNodes));
    }

    // Impactors, keyed every frame (they travel without turning), with the
    // stains on them as their children.
    std::vector<unsigned char> onImpactor(w.spray.stains.size(), 0);
    for (std::size_t si = 0; si < w.spray.stains.size(); ++si) {
        const int c = w.spray.stains[si].collider;
        if (c >= 0 && w.colliders[c].role == Role::Impactor) onImpactor[si] = 1;
    }
    for (std::size_t ii = 0; ii < w.impactors.size(); ++ii) {
        const Collider& c = w.colliders[w.impactors[ii].collider];
        Collider local = c;
        local.pos = {};
        Geo g;
        if (c.kind == ShapeKind::Box) addBox(g, local);
        else addCapsule(g, c.a - c.pos, c.b - c.pos, c.radius);
        std::vector<int> kids;
        if (o.stains) {
            std::vector<std::pair<std::size_t, Geo>> byFrame;
            for (std::size_t si = 0; si < w.spray.stains.size(); ++si) {
                if (!onImpactor[si] || w.spray.stains[si].collider != c.id) continue;
                Stain st = w.spray.stains[si];
                st.pos -= c.pos;
                const std::size_t fk = frameOf(st.time);
                auto it = std::find_if(byFrame.begin(), byFrame.end(), [&](const auto& e) { return e.first == fk; });
                if (it == byFrame.end()) {
                    byFrame.push_back({fk, Geo{}});
                    it = byFrame.end() - 1;
                }
                addStain(it->second, st);
            }
            for (auto& [fk, sg] : byFrame) {
                const std::string nm = c.name + " stains " + std::to_string(fk);
                const int nd = G.node(nm, G.mesh(nm, {G.prim(sg.pos, sg.nrm, sg.idx, mStain)}), {}, {}, fk == 0 ? 1 : 0);
                appear(nd, fk);
                kids.push_back(nd);
            }
        }
        const int mat = G.material(c.name, c.color, 1, 0.4, 0.6);
        const int nd = G.node(c.name, G.mesh(c.name, {G.prim(g.pos, g.nrm, g.idx, mat)}), c.pos, {}, 1, kids);
        if (F > 1) {
            std::vector<float> tr;
            for (std::size_t k = 0; k < F; ++k) {
                const std::vector<float>& p = frames[k]->impactors;
                if (p.size() >= 3 * (ii + 1)) tr.insert(tr.end(), {p[3 * ii], p[3 * ii + 1], p[3 * ii + 2]});
                else tr.insert(tr.end(), {float(c.pos.x), float(c.pos.y), float(c.pos.z)});
            }
            G.channel(nd, "translation", timeAcc, G.floats(tr, 3, 0, false), "LINEAR");
        }
        roots.push_back(nd);
    }

    // Stains, grouped by the frame they appear in.
    if (o.stains && !w.spray.stains.empty()) {
        std::vector<int> stainNodes;
        const std::size_t groups = std::max<std::size_t>(F, 1);
        std::size_t from = 0;
        for (std::size_t k = 0; k < groups; ++k) {
            const double until = k + 1 < F ? frames[k + 1]->time : 1e30;
            Geo g;
            std::size_t to = from;
            while (to < w.spray.stains.size() && (w.spray.stains[to].time < until || k + 1 == groups)) {
                const std::size_t si = to++;
                const Stain& st = w.spray.stains[si];
                if (st.collider >= 0 && w.colliders[st.collider].role == Role::Dummy && w.ragdoll.released()) continue;
                if (onPiece[si] || onImpactor[si] || st.gone) continue;  // with its piece of glass, or its impactor
                addStain(g, st);
            }
            if (to == from) continue;
            from = to;
            if (g.idx.empty()) continue;
            const int mesh = G.mesh("stains_" + std::to_string(k), {G.prim(g.pos, g.nrm, g.idx, mStain)});
            const int nd = G.node("stains_" + std::to_string(k), mesh, {}, {}, k == 0 ? 1 : 0);
            stainNodes.push_back(nd);
            if (k > 0) {
                std::vector<float> sc;
                for (std::size_t q = 0; q < F; ++q) {
                    const float v = q >= k ? 1.0f : 0.0f;
                    sc.insert(sc.end(), {v, v, v});
                }
                G.channel(nd, "scale", timeAcc, G.floats(sc, 3, 0, false), "STEP");
            }
        }
        roots.push_back(G.node("stains", -1, {}, {}, 1, stainNodes));
    }

    // The static scene.
    if (o.scene) {
        std::vector<int> sceneNodes;
        for (const auto& c : w.colliders) {
            if (c.body >= 0 || c.role == Role::Impactor) continue;  // impactors are drawn above
            Geo g;
            if (c.kind == ShapeKind::Plane) {
                const Vec3 n = c.rot.rotate({0, 1, 0}), u = perpendicular(n), v = cross(n, u);
                const double half = c.role == Role::Floor ? 2.5 : 2.0;
                Vec3 ctr = c.pos;
                if (c.role == Role::Wall) ctr.y = std::max(ctr.y, 1.5);
                addQuad(g, ctr, u * half, v * half, n);
            } else if (c.kind == ShapeKind::Box) {
                const PanePieces* bp = w.panePieces(c.id);
                if (bp && !bp->pieces.empty()) continue;  // drawn as its pieces
                addBox(g, c);
            }
            else if (c.kind == ShapeKind::Capsule) {
                if (c.role == Role::Dummy && w.ragdoll.active() && F > 1) {
                    // A ragdoll capsule: modelled along +y about its centre, its pose keyed every frame.
                    const double len = length(c.b - c.a);
                    addCapsule(g, {0, -0.5 * len, 0}, {0, 0.5 * len, 0}, c.radius);
                    const int seg = w.ragdoll.segmentOf(c.id);
                    const int mat = G.material(c.name, c.color, 1, 0.6, 0);
                    const int mesh = G.mesh(c.name, {G.prim(g.pos, g.nrm, g.idx, mat)});
                    std::vector<float> tr, ro;
                    for (std::size_t k = 0; k < F; ++k) {
                        const std::vector<float>& d = frames[k]->dummy;
                        Vec3 a = c.a, b = c.b;
                        if (seg >= 0 && d.size() >= static_cast<std::size_t>(6 * (seg + 1))) {
                            a = {d[6 * seg], d[6 * seg + 1], d[6 * seg + 2]};
                            b = {d[6 * seg + 3], d[6 * seg + 4], d[6 * seg + 5]};
                        }
                        const Vec3 m = (a + b) * 0.5, e = normalize(b - a);
                        const Vec3 ax = cross(Vec3{0, 1, 0}, e);
                        const double ang = std::acos(clamp(e.y, -1.0, 1.0));
                        const Quat q = length(ax) > 1e-9 ? Quat::axisAngle(ax, ang) : (e.y > 0 ? Quat{} : Quat::axisAngle({1, 0, 0}, kPi));
                        tr.insert(tr.end(), {float(m.x), float(m.y), float(m.z)});
                        ro.insert(ro.end(), {float(q.x), float(q.y), float(q.z), float(q.w)});
                    }
                    const int nd = G.node(c.name, mesh, (c.a + c.b) * 0.5);
                    G.channel(nd, "translation", timeAcc, G.floats(tr, 3, 0, false), "LINEAR");
                    G.channel(nd, "rotation", timeAcc, G.floats(ro, 4, 0, false), "LINEAR");
                    sceneNodes.push_back(nd);
                    continue;
                }
                addCapsule(g, c.a, c.b, c.radius);
            } else continue;
            const bool glass = c.role == Role::Pane;
            const int mat = glass ? mGlass : G.material(c.name, c.color, 1, 0.7, 0);
            const int mesh = G.mesh(c.name, {G.prim(g.pos, g.nrm, g.idx, mat)});
            sceneNodes.push_back(G.node(c.name, mesh));
        }
        if (!sceneNodes.empty()) roots.push_back(G.node("scene", -1, {}, {}, 1, sceneNodes));
    }
    return G.finish(roots);
}

std::vector<std::uint8_t> exportPly(const RecordedFrame& f, double spacing) {
    const std::size_t nb = f.particles.size() / 3, nd = f.drops.size() / 4;
    std::string h = "ply\nformat binary_little_endian 1.0\ncomment Liquid Dynamics point cache, metres, +y up, t = " +
                    num(f.time) + " s\nelement vertex " + std::to_string(nb + nd) +
                    "\nproperty float x\nproperty float y\nproperty float z\nproperty float pscale\nproperty uchar kind\nend_header\n";
    std::vector<std::uint8_t> out(h.begin(), h.end());
    auto put = [&](float v) {
        std::uint8_t b[4];
        std::memcpy(b, &v, 4);
        out.insert(out.end(), b, b + 4);
    };
    const float rb = static_cast<float>(0.5 * spacing);
    for (std::size_t i = 0; i < nb; ++i) {
        put(f.particles[3 * i]);
        put(f.particles[3 * i + 1]);
        put(f.particles[3 * i + 2]);
        put(rb);
        out.push_back(0);
    }
    for (std::size_t i = 0; i < nd; ++i) {
        put(f.drops[4 * i]);
        put(f.drops[4 * i + 1]);
        put(f.drops[4 * i + 2]);
        put(0.5f * f.drops[4 * i + 3]);
        out.push_back(1);
    }
    return out;
}

std::string exportObj(const SurfaceMesh& m) {
    std::string o = "# Liquid Dynamics liquid surface, metres, +y up\n";
    char b[96];
    for (std::size_t i = 0; i < m.vertices(); ++i) {
        std::snprintf(b, sizeof b, "v %.6f %.6f %.6f\n", m.pos[3 * i], m.pos[3 * i + 1], m.pos[3 * i + 2]);
        o += b;
    }
    for (std::size_t i = 0; i < m.vertices(); ++i) {
        std::snprintf(b, sizeof b, "vn %.4f %.4f %.4f\n", m.nrm[3 * i], m.nrm[3 * i + 1], m.nrm[3 * i + 2]);
        o += b;
    }
    for (std::size_t t = 0; t < m.triangles(); ++t) {
        const unsigned a = m.tri[3 * t] + 1, c = m.tri[3 * t + 1] + 1, d = m.tri[3 * t + 2] + 1;
        std::snprintf(b, sizeof b, "f %u//%u %u//%u %u//%u\n", a, a, c, c, d, d);
        o += b;
    }
    return o;
}

std::string exportStainsJson(const World& w) {
    const LiquidDef* L = w.spray.liquid();
    std::string o = "{\"units\":\"m\",\"frame\":\"glTF: right-handed, +y up\",\"liquid\":\"" + std::string(L ? L->id : "") + "\"";
    if (L)
        o += ",\"stain_color\":[" + num(L->stainColor.x) + "," + num(L->stainColor.y) + "," + num(L->stainColor.z) + "]";
    o += ",\"stains\":[";
    for (std::size_t i = 0; i < w.spray.stains.size(); ++i) {
        const Stain& s = w.spray.stains[i];
        const std::string surf = s.collider >= 0 ? w.colliders[s.collider].name : "";
        o += (i ? "," : "") + std::string("{\"p\":[") + num(s.pos.x) + "," + num(s.pos.y) + "," + num(s.pos.z) +
             "],\"n\":[" + num(s.normal.x) + "," + num(s.normal.y) + "," + num(s.normal.z) + "],\"tail\":[" +
             num(s.dir.x) + "," + num(s.dir.y) + "," + num(s.dir.z) + "],\"width\":" + num(s.width) +
             ",\"length\":" + num(s.length) + ",\"alpha_deg\":" + num(s.alpha * 180 / kPi) + ",\"drop_mm\":" +
             num(s.diameter * 1e3) + ",\"speed\":" + num(s.speed) + ",\"t\":" + num(s.time) + ",\"surface\":\"" +
             surf + "\",\"satellite\":" + (s.satellite ? "true" : "false") + ",\"spines\":" + std::to_string(s.spines) +
             ",\"source\":" + std::to_string(s.source) + "}";
    }
    o += "]}";
    return o;
}

}  // namespace lqd
