// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Uniform-grid neighbour search by spatial hashing (Teschner et al. 2003),
// rebuilt every step with a counting sort, so it costs O(N) and has no upper
// bound on the domain.
#pragma once

#include <cstdint>
#include <vector>

#include "lqd/math.hpp"

namespace lqd {

class NeighbourGrid {
public:
    // Bins the points with cells of size `cell` (normally the kernel support).
    void build(const std::vector<Vec3>& x, double cell);
    double cell() const { return cell_; }
    bool empty() const { return order_.empty(); }

    // Calls f(j) for every binned point j with |x_j - p| < radius (radius <= cell).
    template <class F>
    void query(const std::vector<Vec3>& x, const Vec3& p, double radius, F&& f) const {
        if (order_.empty()) return;
        const double r2 = radius * radius;
        const int cx = coord(p.x), cy = coord(p.y), cz = coord(p.z);
        std::uint32_t seen[27];
        int nseen = 0;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const std::uint32_t b = bucket(cx + dx, cy + dy, cz + dz);
                    bool dup = false;
                    for (int s = 0; s < nseen; ++s) dup = dup || seen[s] == b;
                    if (dup) continue;
                    seen[nseen++] = b;
                    for (std::uint32_t k = start_[b], e = start_[b + 1]; k < e; ++k) {
                        const int j = order_[k];
                        if (length2(x[j] - p) < r2) f(j);
                    }
                }
    }

    // Compressed neighbour lists of every binned point (itself excluded).
    void neighbourLists(const std::vector<Vec3>& x, double radius, std::vector<int>& offsets,
                        std::vector<int>& list) const;

private:
    int coord(double v) const { return static_cast<int>(std::floor(v * inv_)); }
    std::uint32_t bucket(int x, int y, int z) const {
        const std::uint32_t h = (static_cast<std::uint32_t>(x) * 73856093u) ^ (static_cast<std::uint32_t>(y) * 19349663u) ^
                                (static_cast<std::uint32_t>(z) * 83492791u);
        // Fibonacci hashing: the high bits of the product are well mixed.
        return static_cast<std::uint32_t>((h * 2654435769u) >> shift_);
    }

    double cell_ = 1, inv_ = 1;
    std::uint32_t mask_ = 0;
    int shift_ = 26;
    std::vector<std::uint32_t> start_;
    std::vector<int> order_;
    std::vector<std::uint32_t> keys_;
};

}  // namespace lqd
