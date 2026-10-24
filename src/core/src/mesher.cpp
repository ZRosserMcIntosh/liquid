// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/mesher.hpp"

#include <algorithm>
#include <cmath>

namespace lqd {

namespace {

// Open-addressing hash map from 64-bit keys to 32-bit slots.
class FlatMap {
public:
    explicit FlatMap(std::size_t expected) {
        std::size_t n = 1024;
        while (n < expected * 2) n <<= 1;
        keys_.assign(n, kEmpty);
        vals_.assign(n, 0);
        mask_ = n - 1;
    }
    // Returns the slot for key, inserting it with `fresh` if absent.
    std::uint32_t& at(std::uint64_t key, std::uint32_t fresh, bool* inserted = nullptr) {
        if (count_ * 2 >= keys_.size()) grow();
        std::size_t i = hash(key) & mask_;
        while (keys_[i] != kEmpty && keys_[i] != key) i = (i + 1) & mask_;
        if (keys_[i] == kEmpty) {
            keys_[i] = key;
            vals_[i] = fresh;
            ++count_;
            if (inserted) *inserted = true;
        } else if (inserted) *inserted = false;
        return vals_[i];
    }
    const std::uint32_t* find(std::uint64_t key) const {
        std::size_t i = hash(key) & mask_;
        while (keys_[i] != kEmpty) {
            if (keys_[i] == key) return &vals_[i];
            i = (i + 1) & mask_;
        }
        return nullptr;
    }

private:
    static constexpr std::uint64_t kEmpty = ~0ull;
    static std::size_t hash(std::uint64_t k) {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdull;
        k ^= k >> 33;
        return static_cast<std::size_t>(k);
    }
    void grow() {
        std::vector<std::uint64_t> ok = std::move(keys_);
        std::vector<std::uint32_t> ov = std::move(vals_);
        keys_.assign(ok.size() * 2, kEmpty);
        vals_.assign(ok.size() * 2, 0);
        mask_ = keys_.size() - 1;
        count_ = 0;
        for (std::size_t i = 0; i < ok.size(); ++i)
            if (ok[i] != kEmpty) at(ok[i], ov[i]);
    }
    std::vector<std::uint64_t> keys_;
    std::vector<std::uint32_t> vals_;
    std::size_t mask_ = 0, count_ = 0;
};

constexpr int kBits = 19;
constexpr int kOff = 1 << (kBits - 1);
std::uint64_t pack(int i, int j, int k) {
    return (static_cast<std::uint64_t>(i + kOff) << (2 * kBits)) | (static_cast<std::uint64_t>(j + kOff) << kBits) |
           static_cast<std::uint64_t>(k + kOff);
}

constexpr int kCorner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
constexpr int kTet[6][4] = {{0, 6, 1, 2}, {0, 6, 2, 3}, {0, 6, 3, 7}, {0, 6, 7, 4}, {0, 6, 4, 5}, {0, 6, 5, 1}};

}  // namespace

SurfaceMesh meshLiquid(const std::vector<Vec3>& particles, double spacing, const MesherSettings& set) {
    SurfaceMesh out;
    if (particles.empty()) return out;
    const double c = set.cell > 0 ? set.cell : 0.6 * spacing;
    const double R = set.radius > 0 ? set.radius : 1.8 * spacing;
    const double inv = 1.0 / c;
    auto splat = [&](double q2) {
        const double a = 1 - q2;
        return a * a * a;
    };
    // The field of a resting cubic lattice, averaged over the node's offset.
    double rest = 0;
    {
        int samples = 0;
        for (double ox : {0.0, 0.5})
            for (double oy : {0.0, 0.5}) {
                double f = 0;
                const int n = static_cast<int>(std::ceil(R / spacing)) + 1;
                for (int i = -n; i <= n; ++i)
                    for (int j = -n; j <= n; ++j)
                        for (int k = -n; k <= n; ++k) {
                            const double q2 = (sq(i + ox) + sq(j + oy) + sq(k + 0.25)) * spacing * spacing / (R * R);
                            if (q2 < 1) f += splat(q2);
                        }
                rest += f;
                ++samples;
            }
        rest /= samples;
    }
    const float iso = static_cast<float>(set.iso * rest);

    // Splat onto the sparse grid.
    FlatMap nodes(particles.size() * 64);
    std::vector<float> value;
    std::vector<int> ni, nj, nk;
    const int span = static_cast<int>(std::ceil(R * inv));
    for (const auto& p : particles) {
        const int ci = static_cast<int>(std::floor(p.x * inv)), cj = static_cast<int>(std::floor(p.y * inv)),
                  ck = static_cast<int>(std::floor(p.z * inv));
        for (int i = ci - span + 1; i <= ci + span; ++i)
            for (int j = cj - span + 1; j <= cj + span; ++j)
                for (int k = ck - span + 1; k <= ck + span; ++k) {
                    const double q2 = (sq(i * c - p.x) + sq(j * c - p.y) + sq(k * c - p.z)) / (R * R);
                    if (q2 >= 1) continue;
                    bool fresh = false;
                    std::uint32_t& slot = nodes.at(pack(i, j, k), static_cast<std::uint32_t>(value.size()), &fresh);
                    if (fresh) {
                        value.push_back(0);
                        ni.push_back(i);
                        nj.push_back(j);
                        nk.push_back(k);
                    }
                    value[slot] += static_cast<float>(splat(q2));
                }
    }
    auto field = [&](int i, int j, int k) {
        const std::uint32_t* s = nodes.find(pack(i, j, k));
        return s ? value[*s] : 0.0f;
    };

    FlatMap cubes(value.size() * 2);
    FlatMap edges(value.size());
    std::vector<Vec3> acc;  // area-weighted normals
    auto vertexOn = [&](int ai, int aj, int ak, float fa, int bi, int bj, int bk, float fb) {
        // The edge key: its lower node and the direction to the other one.
        if (pack(bi, bj, bk) < pack(ai, aj, ak)) {
            std::swap(ai, bi);
            std::swap(aj, bj);
            std::swap(ak, bk);
            std::swap(fa, fb);
        }
        const std::uint64_t key = (pack(ai, aj, ak) << 5) |
                                  static_cast<std::uint64_t>((bi - ai + 1) * 9 + (bj - aj + 1) * 3 + (bk - ak + 1));
        bool fresh = false;
        std::uint32_t& slot = edges.at(key, static_cast<std::uint32_t>(out.pos.size() / 3), &fresh);
        if (fresh) {
            const float t = std::abs(fb - fa) > 1e-12f ? (iso - fa) / (fb - fa) : 0.5f;
            const double td = t;
            out.pos.push_back(static_cast<float>((ai + td * (bi - ai)) * c));
            out.pos.push_back(static_cast<float>((aj + td * (bj - aj)) * c));
            out.pos.push_back(static_cast<float>((ak + td * (bk - ak)) * c));
            acc.push_back({});
        }
        return slot;
    };
    auto emit = [&](std::uint32_t a, std::uint32_t b, std::uint32_t d, const Vec3& outward) {
        auto P = [&](std::uint32_t v) {
            return Vec3{out.pos[3 * v], out.pos[3 * v + 1], out.pos[3 * v + 2]};
        };
        Vec3 n = cross(P(b) - P(a), P(d) - P(a));
        if (dot(n, outward) < 0) {
            std::swap(b, d);
            n = -n;
        }
        if (length2(n) < 1e-30) return;
        out.tri.insert(out.tri.end(), {a, b, d});
        acc[a] += n;
        acc[b] += n;
        acc[d] += n;
    };

    for (std::size_t s = 0; s < value.size(); ++s)
        for (int dz = -1; dz <= 0; ++dz)
            for (int dy = -1; dy <= 0; ++dy)
                for (int dx = -1; dx <= 0; ++dx) {
                    const int i0 = ni[s] + dx, j0 = nj[s] + dy, k0 = nk[s] + dz;
                    bool fresh = false;
                    cubes.at(pack(i0, j0, k0), 0, &fresh);
                    if (!fresh) continue;
                    float f[8];
                    int inside = 0;
                    for (int q = 0; q < 8; ++q) {
                        f[q] = field(i0 + kCorner[q][0], j0 + kCorner[q][1], k0 + kCorner[q][2]);
                        inside += f[q] > iso ? 1 : 0;
                    }
                    if (inside == 0 || inside == 8) continue;
                    // Outward direction: down the field's gradient across the cube.
                    Vec3 grad{(f[1] + f[2] + f[5] + f[6]) - (f[0] + f[3] + f[4] + f[7]),
                              (f[2] + f[3] + f[6] + f[7]) - (f[0] + f[1] + f[4] + f[5]),
                              (f[4] + f[5] + f[6] + f[7]) - (f[0] + f[1] + f[2] + f[3])};
                    const Vec3 outward = -grad;
                    for (const auto& T : kTet) {
                        int in[4], ins = 0, outs = 0, o[4];
                        for (int q = 0; q < 4; ++q) {
                            if (f[T[q]] > iso) in[ins++] = T[q];
                            else o[outs++] = T[q];
                        }
                        if (ins == 0 || ins == 4) continue;
                        auto V = [&](int a, int b) {
                            return vertexOn(i0 + kCorner[a][0], j0 + kCorner[a][1], k0 + kCorner[a][2], f[a],
                                            i0 + kCorner[b][0], j0 + kCorner[b][1], k0 + kCorner[b][2], f[b]);
                        };
                        if (ins == 1) emit(V(in[0], o[0]), V(in[0], o[1]), V(in[0], o[2]), outward);
                        else if (ins == 3) emit(V(o[0], in[0]), V(o[0], in[1]), V(o[0], in[2]), outward);
                        else {
                            const std::uint32_t a = V(in[0], o[0]), b = V(in[0], o[1]), d = V(in[1], o[1]),
                                                e = V(in[1], o[0]);
                            emit(a, b, d, outward);
                            emit(a, d, e, outward);
                        }
                    }
                }
    out.nrm.resize(out.pos.size());
    for (std::size_t v = 0; v < acc.size(); ++v) {
        const Vec3 n = normalize(acc[v]);
        out.nrm[3 * v] = static_cast<float>(n.x);
        out.nrm[3 * v + 1] = static_cast<float>(n.y);
        out.nrm[3 * v + 2] = static_cast<float>(n.z);
    }
    return out;
}

}  // namespace lqd
