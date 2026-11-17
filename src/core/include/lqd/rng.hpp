// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Deterministic random numbers (SplitMix64). Every scenario takes a seed, and
// the same seed gives the same splash, spray and shards on every platform.
#pragma once

#include <cmath>
#include <cstdint>

#include "lqd/math.hpp"

namespace lqd {

class Rng {
public:
    explicit Rng(std::uint64_t seed = 1) : s_(seed * 0x9E3779B97F4A7C15ull + 0x2545F4914F6CDD1Dull) {}
    std::uint64_t next() {
        std::uint64_t z = (s_ += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    // Uniform in [0, 1).
    double uniform() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }
    double uniform(double lo, double hi) { return lo + (hi - lo) * uniform(); }
    double normal() {
        const double u1 = std::max(uniform(), 1e-300), u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
    }
    // Log-normal with the given median and log standard deviation.
    double lognormal(double median, double sigmaLn) { return median * std::exp(sigmaLn * normal()); }
    // A unit vector inside the cone of half-angle `halfAngle` around unit `axis`,
    // uniform over the cone's solid angle.
    Vec3 cone(const Vec3& axis, double halfAngle) {
        const double c = 1.0 - uniform() * (1.0 - std::cos(halfAngle));
        const double s = std::sqrt(std::max(0.0, 1.0 - c * c)), phi = 2.0 * kPi * uniform();
        const Vec3 u = perpendicular(axis), v = cross(axis, u);
        return axis * c + u * (s * std::cos(phi)) + v * (s * std::sin(phi));
    }

private:
    std::uint64_t s_;
};

}  // namespace lqd
