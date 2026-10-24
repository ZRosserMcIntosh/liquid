// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// The SPH smoothing kernel and the volume a solid wall fills inside it.
//
// The kernel is the cubic B-spline (Monaghan 1992) written with its support
// radius H, normalised so that its integral over the ball is 1. A wall is
// handled analytically rather than with layers of wall particles: a fluid
// particle at distance d from a solid surface sees the fraction Phi(d) of its
// kernel inside the solid (the half-space integral of W, tabulated once), so
// the wall adds rho0 Phi(d) to its density and pushes along the surface normal
// with -A(d), the kernel's cross-section at depth d. This is the
// continuum limit of Akinci et al. (2012) boundary particles; it is exact for
// flat walls and good to O(H/R) on walls of curvature radius R.
#pragma once

#include <array>
#include <cmath>

#include "lqd/math.hpp"

namespace lqd {

class Kernel {
public:
    explicit Kernel(double supportRadius = 0.01) { setSupport(supportRadius); }

    void setSupport(double H) {
        H_ = H;
        k_ = 8.0 / (kPi * H * H * H);
        l_ = 48.0 / (kPi * H * H * H);
        buildWallTable();
    }
    double support() const { return H_; }

    double W(double r) const {
        const double q = r / H_;
        if (q >= 1) return 0;
        if (q <= 0.5) return k_ * (6 * q * q * q - 6 * q * q + 1);
        const double a = 1 - q;
        return k_ * 2 * a * a * a;
    }
    // Gradient with respect to the first point, for rij = xi - xj, r = |rij|.
    Vec3 gradW(const Vec3& rij, double r) const {
        if (r <= 1e-12 * H_ || r >= H_) return {};
        const double q = r / H_;
        const double f = q <= 0.5 ? l_ * q * (3 * q - 2) : -l_ * (1 - q) * (1 - q);
        return rij * (f / (r * H_));
    }

    // Fraction of the kernel inside a solid half-space whose surface is at
    // signed distance d (positive: the particle is outside the solid).
    double wallVolume(double d) const { return lookup(phi_, d, 1.0, 0.0); }
    // The kernel's area-weighted cross-section A(d) at depth d: d Phi / d d = -A(d).
    double wallSection(double d) const { return lookup(section_, d, 0.0, 0.0); }

private:
    static constexpr int kTable = 512;
    double lookup(const std::array<double, kTable + 1>& t, double d, double below, double above) const {
        if (d <= -H_) return below;
        if (d >= H_) return above;
        const double u = (d + H_) / (2 * H_) * kTable;
        const int i = static_cast<int>(u);
        if (i >= kTable) return t[kTable];
        const double f = u - i;
        return t[i] * (1 - f) + t[i + 1] * f;
    }
    // A(z) = integral of W(sqrt(rho^2 + z^2)) 2 pi rho d rho over the disc at height z.
    double crossSection(double z) const {
        const double R2 = H_ * H_ - z * z;
        if (R2 <= 0) return 0;
        const double R = std::sqrt(R2);
        const int n = 64;
        double s = 0;
        for (int i = 0; i <= n; ++i) {
            const double rho = R * i / n;
            const double w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
            s += w * W(std::sqrt(rho * rho + z * z)) * 2 * kPi * rho;
        }
        return s * (R / n) / 3;
    }
    void buildWallTable() {
        // section_[i] = A(d_i); phi_[i] = integral from z = -H to -d of A(z) dz.
        for (int i = 0; i <= kTable; ++i) section_[i] = crossSection(-H_ + 2 * H_ * i / kTable);
        // Phi(d) for d from -H (fully inside, 1) to H (0): cumulative trapezoid from the top.
        phi_[kTable] = 0;
        const double dz = 2 * H_ / kTable;
        for (int i = kTable - 1; i >= 0; --i) phi_[i] = phi_[i + 1] + 0.5 * (section_[i] + section_[i + 1]) * dz;
        const double total = phi_[0];
        if (total > 0)
            for (auto& p : phi_) p /= total;
        for (auto& a : section_) a /= total;
    }

    double H_ = 0, k_ = 0, l_ = 0;
    std::array<double, kTable + 1> phi_{}, section_{};
};

}  // namespace lqd
