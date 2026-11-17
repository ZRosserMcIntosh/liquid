// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/grid.hpp"

namespace lqd {

void NeighbourGrid::build(const std::vector<Vec3>& x, double cell) {
    cell_ = cell;
    inv_ = 1.0 / cell;
    const std::size_t n = x.size();
    std::uint32_t size = 64;
    int bits = 6;
    while (size < 2 * n) {
        size <<= 1;
        ++bits;
    }
    mask_ = size - 1;
    shift_ = 32 - bits;
    start_.assign(size + 1, 0);
    keys_.resize(n);
    order_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        keys_[i] = bucket(coord(x[i].x), coord(x[i].y), coord(x[i].z));
        ++start_[keys_[i] + 1];
    }
    for (std::uint32_t b = 0; b < size; ++b) start_[b + 1] += start_[b];
    std::vector<std::uint32_t> fill(start_.begin(), start_.end() - 1);
    for (std::size_t i = 0; i < n; ++i) order_[fill[keys_[i]]++] = static_cast<int>(i);
}

void NeighbourGrid::neighbourLists(const std::vector<Vec3>& x, double radius, std::vector<int>& offsets,
                                   std::vector<int>& list) const {
    const int n = static_cast<int>(x.size());
    offsets.assign(n + 1, 0);
    list.clear();
    for (int i = 0; i < n; ++i) {
        query(x, x[i], radius, [&](int j) {
            if (j != i) list.push_back(j);
        });
        offsets[i + 1] = static_cast<int>(list.size());
    }
}

}  // namespace lqd
