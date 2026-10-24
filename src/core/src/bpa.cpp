// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/bpa.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace lqd {

namespace {
double quantile(std::vector<double> v, double q) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const double pos = q * static_cast<double>(v.size() - 1);
    const std::size_t i = static_cast<std::size_t>(pos);
    const double f = pos - static_cast<double>(i);
    return i + 1 < v.size() ? v[i] * (1 - f) + v[i + 1] * f : v[i];
}
}  // namespace

BpaResult analyseStains(const Spray& spray, const std::vector<Collider>& colliders, const BpaOptions& o) {
    BpaResult R;
    if (o.collider < 0 || o.collider >= static_cast<int>(colliders.size())) {
        R.summary = "no such surface";
        return R;
    }
    // The face with most stains (a pane has two).
    std::vector<int> pick;
    Vec3 nsum;
    for (int i = 0; i < static_cast<int>(spray.stains.size()); ++i) {
        const Stain& s = spray.stains[i];
        if (s.collider != o.collider) continue;
        if (o.source >= 0 && s.source != o.source) continue;
        if (s.satellite && !o.satellites) continue;
        // Only stains still where they landed: not those on glass that has fallen.
        if (s.carried || s.gone) continue;
        pick.push_back(i);
    }
    if (pick.empty()) {
        R.summary = "no stains on that surface";
        return R;
    }
    // Majority normal.
    Vec3 best = spray.stains[pick[0]].normal;
    int bestCount = 0;
    for (int a : pick) {
        int c = 0;
        for (int b : pick) c += dot(spray.stains[a].normal, spray.stains[b].normal) > 0.95 ? 1 : 0;
        if (c > bestCount) {
            bestCount = c;
            best = spray.stains[a].normal;
        }
        if (bestCount * 2 > static_cast<int>(pick.size())) break;
    }
    (void)nsum;
    const Vec3 n = normalize(best);
    const Vec3 P0 = spray.stains[pick[0]].pos;
    const Vec3 u = perpendicular(n), w = cross(n, u);
    R.planePoint = P0;
    R.planeNormal = n;
    // Lines.
    double A[2][2] = {{0, 0}, {0, 0}}, c[2] = {0, 0};
    std::vector<double> heights;
    for (int i : pick) {
        const Stain& s = spray.stains[i];
        if (dot(s.normal, n) < 0.95) continue;
        ++R.considered;
        double W = s.width, L = s.length;
        if (o.readError > 0) {
            W = std::round(W / o.readError) * o.readError;
            L = std::round(L / o.readError) * o.readError;
        }
        if (L < o.minLength || W <= 0) continue;
        const double ratio = std::min(W / L, 1.0);
        if (ratio > o.maxRatio) continue;
        const double alpha = std::asin(ratio);
        if (alpha < o.minAlpha) continue;
        BpaLine l;
        l.stain = i;
        l.pos = s.pos;
        l.back = normalize(-(s.dir - n * dot(s.dir, n)));
        l.ratio = ratio;
        l.alpha = alpha;
        l.alphaTrue = s.alpha;
        R.lines.push_back(l);
        const double px = dot(s.pos - P0, u), py = dot(s.pos - P0, w);
        const double bx = dot(l.back, u), by = dot(l.back, w);
        // (I - b b^T)
        const double m00 = 1 - bx * bx, m01 = -bx * by, m11 = 1 - by * by;
        A[0][0] += m00;
        A[0][1] += m01;
        A[1][1] += m11;
        c[0] += m00 * px + m01 * py;
        c[1] += m01 * px + m11 * py;
    }
    R.used = static_cast<int>(R.lines.size());
    if (R.used < 3) {
        R.summary = "too few measurable stains (" + std::to_string(R.used) + ")";
        return R;
    }
    A[1][0] = A[0][1];
    const double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (std::abs(det) < 1e-9) {
        R.summary = "the stain lines are parallel: no convergence";
        return R;
    }
    const double qx = (c[0] * A[1][1] - c[1] * A[0][1]) / det, qy = (A[0][0] * c[1] - A[1][0] * c[0]) / det;
    R.convergence = P0 + u * qx + w * qy;
    double rms = 0, angErr = 0;
    for (auto& l : R.lines) {
        const Vec3 d = R.convergence - l.pos;
        const Vec3 dIn = d - n * dot(d, n);
        // Perpendicular miss of this line.
        const Vec3 miss = dIn - l.back * dot(dIn, l.back);
        rms += length2(miss);
        l.distance = std::max(0.0, dot(dIn, l.back));
        l.height = l.distance * std::tan(l.alpha);
        heights.push_back(l.height);
        angErr += std::abs(l.alpha - l.alphaTrue);
    }
    R.convergenceRms = std::sqrt(rms / R.used);
    R.meanAngleError = angErr / R.used;
    R.height = quantile(heights, 0.5);
    R.heightLow = quantile(heights, 0.25);
    R.heightHigh = quantile(heights, 0.75);
    R.origin = R.convergence + n * R.height;
    R.ok = true;
    // Truth: the source of most of the used stains.
    std::vector<int> count(spray.sources.size(), 0);
    for (const auto& l : R.lines) {
        const int src = spray.stains[l.stain].source;
        if (src >= 0 && src < static_cast<int>(count.size())) ++count[src];
    }
    int src = -1, cnt = 0;
    for (int k = 0; k < static_cast<int>(count.size()); ++k)
        if (count[k] > cnt) {
            cnt = count[k];
            src = k;
        }
    if (src >= 0) {
        R.truth = true;
        R.trueOrigin = spray.sources[src].origin;
        const Vec3 d = R.trueOrigin - R.convergence;
        R.convergenceError = length(d - n * dot(d, n));
        R.heightError = R.height - dot(R.trueOrigin - R.convergence, n);
    }
    char buf[400];
    if (R.truth)
        std::snprintf(buf, sizeof buf,
                      "%d of %d stains measured; convergence %.0f mm from the true source's foot (lines miss it by "
                      "%.0f mm rms); tangent-method origin %.2f m from the surface (IQR %.2f-%.2f), true %.2f m: %s by "
                      "%.0f mm; impact angles read to %.1f deg on average.",
                      R.used, R.considered, R.convergenceError * 1e3, R.convergenceRms * 1e3, R.height, R.heightLow,
                      R.heightHigh, R.height - R.heightError, R.heightError >= 0 ? "too far" : "too near",
                      std::abs(R.heightError) * 1e3, R.meanAngleError * 180 / kPi);
    else
        std::snprintf(buf, sizeof buf, "%d of %d stains measured; origin %.2f m from the surface (IQR %.2f-%.2f).",
                      R.used, R.considered, R.height, R.heightLow, R.heightHigh);
    R.summary = buf;
    return R;
}

}  // namespace lqd
