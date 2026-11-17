// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
#include "lqd/vessel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace lqd {

namespace {

// Volume and radial centroid of a wall element of thickness t and length ds
// at mid radius r whose section is tilted so it spans `e` either side
// radially; the part that would cross the axis is clipped.
void wallElement(double r, double t, double ds, double e, double& volume, double& rc) {
    if (r >= e || e <= 0) {
        volume = 2 * kPi * r * t * ds;
        rc = r;
        return;
    }
    const double frac = (r + e) / (2 * e);
    rc = 0.5 * (r + e);
    volume = 2 * kPi * (t * ds * frac) * rc;
}

double wrapAngle(double a) {
    while (a > kPi) a -= 2 * kPi;
    while (a < -kPi) a += 2 * kPi;
    return a;
}

}  // namespace

VesselGeometry::VesselGeometry(const VesselDef& def, int columns) : def_(&def), p_(def.profile), columns_(columns) {
    s_.assign(p_.size(), 0);
    for (std::size_t k = 1; k < p_.size(); ++k)
        s_[k] = s_[k - 1] + std::hypot(p_[k].r - p_[k - 1].r, p_[k].y - p_[k - 1].y);
    for (const auto& q : p_) {
        height_ = std::max(height_, q.y + 0.5 * q.t);
        maxR_ = std::max(maxR_, q.r + 0.5 * q.t);
    }
    resample();
    cavityPolygon();
    massProperties();
    buildSamples();
}

// ---- Profile ------------------------------------------------------------------

void VesselGeometry::at(double s, double& r, double& y, double& nr, double& ny, double& t) const {
    const int n = static_cast<int>(p_.size());
    s = clamp(s, 0.0, s_.back());
    int k = static_cast<int>(std::upper_bound(s_.begin(), s_.end(), s) - s_.begin()) - 1;
    k = std::clamp(k, 0, n - 2);
    const double len = std::max(s_[k + 1] - s_[k], 1e-12);
    const double f = clamp((s - s_[k]) / len, 0.0, 1.0);
    r = p_[k].r + f * (p_[k + 1].r - p_[k].r);
    y = p_[k].y + f * (p_[k + 1].y - p_[k].y);
    t = p_[k].t + f * (p_[k + 1].t - p_[k].t);
    // Segment normals, averaged at the vertices and blended along the segment.
    auto segNormal = [&](int j, double& a, double& b) {
        j = std::clamp(j, 0, n - 2);
        const double tx = p_[j + 1].r - p_[j].r, ty = p_[j + 1].y - p_[j].y, l = std::max(std::hypot(tx, ty), 1e-12);
        a = -ty / l;
        b = tx / l;
    };
    double a0, b0, a1, b1, am, bm, ap, bp;
    segNormal(k, a0, b0);
    segNormal(k - 1, am, bm);
    segNormal(k + 1, ap, bp);
    a1 = a0;
    b1 = b0;
    double va = k > 0 ? a0 + am : a0, vb = k > 0 ? b0 + bm : b0;
    double wa = k + 1 < n - 1 ? a1 + ap : a1, wb = k + 1 < n - 1 ? b1 + bp : b1;
    double l = std::max(std::hypot(va, vb), 1e-12);
    va /= l;
    vb /= l;
    l = std::max(std::hypot(wa, wb), 1e-12);
    wa /= l;
    wb /= l;
    nr = va + f * (wa - va);
    ny = vb + f * (wb - vb);
    l = std::max(std::hypot(nr, ny), 1e-12);
    nr /= l;
    ny /= l;
}

void VesselGeometry::resample() {
    const double target = 0.0025;  // m between rows
    rowS_.clear();
    for (std::size_t k = 0; k + 1 < p_.size(); ++k) {
        const double len = s_[k + 1] - s_[k];
        const int n = std::max(1, static_cast<int>(std::ceil(len / target)));
        for (int i = 0; i < n; ++i) rowS_.push_back(s_[k] + len * i / n);
    }
    rowS_.push_back(s_.back());
    const std::size_t m = rowS_.size();
    rowR_.resize(m);
    rowY_.resize(m);
    rowNr_.resize(m);
    rowNy_.resize(m);
    rowT_.resize(m);
    for (std::size_t i = 0; i < m; ++i) at(rowS_[i], rowR_[i], rowY_[i], rowNr_[i], rowNy_[i], rowT_[i]);
}

WallPoint VesselGeometry::nearest(const Vec3& p) const {
    WallPoint w;
    const double rho = std::sqrt(p.x * p.x + p.z * p.z);
    const double theta = std::atan2(p.z, p.x);
    double best = std::numeric_limits<double>::max();
    int bk = 0;
    double bf = 0;
    bool bmirror = false;
    for (std::size_t k = 0; k + 1 < p_.size(); ++k) {
        const double ar = p_[k].r, ay = p_[k].y, br = p_[k + 1].r, by = p_[k + 1].y;
        const double dr = br - ar, dy = by - ay, l2 = dr * dr + dy * dy;
        for (int mirror = 0; mirror < 2; ++mirror) {
            const double qr = mirror ? -rho : rho;
            const double f = l2 > 0 ? clamp(((qr - ar) * dr + (p.y - ay) * dy) / l2, 0.0, 1.0) : 0.0;
            const double cr = ar + f * dr, cy = ay + f * dy;
            const double d2 = sq(qr - cr) + sq(p.y - cy);
            if (d2 < best) {
                best = d2;
                bk = static_cast<int>(k);
                bf = f;
                bmirror = mirror != 0;
            }
            if (rho > 4 * maxR_) break;  // far from the axis the mirror image never wins
        }
    }
    const ProfilePoint &A = p_[bk], &B = p_[bk + 1];
    const double cr = A.r + bf * (B.r - A.r), cy = A.y + bf * (B.y - A.y);
    w.dist = std::sqrt(best);
    w.s = s_[bk] + bf * (s_[bk + 1] - s_[bk]);
    w.t = A.t + bf * (B.t - A.t);
    w.theta = bmirror ? wrapAngle(theta + kPi) : theta;
    // 3D closest point and direction to p.
    const double ct = std::cos(w.theta), st = std::sin(w.theta);
    const Vec3 c{cr * ct, cy, cr * st};
    const double sdx = B.r - A.r, sdy = B.y - A.y, sl = std::max(std::hypot(sdx, sdy), 1e-12);
    const double snr = -sdy / sl, sny = sdx / sl;  // segment outward normal in the meridian plane
    const Vec3 segN{snr * ct, sny, snr * st};
    if (w.dist > 1e-9) {
        w.normal = (p - c) / w.dist;
        w.side = dot(w.normal, segN) >= 0 ? 1 : -1;
    } else {
        w.normal = segN;
        w.side = 1;
    }
    return w;
}

double VesselGeometry::sdf(const Vec3& p, Vec3* normal, WallPoint* wp) const {
    const WallPoint w = nearest(p);
    if (normal) *normal = w.normal;
    if (wp) *wp = w;
    return w.dist - 0.5 * w.t;
}

// ---- Cavity -------------------------------------------------------------------

void VesselGeometry::cavityPolygon() {
    polyR_.clear();
    polyY_.clear();
    const int ce = def_->cavityEnd;
    for (int k = 0; k <= ce; ++k) {
        double r, y, nr, ny, t;
        at(s_[k], r, y, nr, ny, t);
        polyR_.push_back(std::max(0.0, r - nr * 0.5 * t));
        polyY_.push_back(y - ny * 0.5 * t);
    }
    if (polyR_.back() > 0) {
        polyR_.push_back(0);
        polyY_.push_back(polyY_.back());
    }
    rimY_ = polyY_.front();
    polyR_.push_back(0);
    polyY_.push_back(rimY_);
}

bool VesselGeometry::inCavity(double r, double y, double margin) const {
    const std::size_t n = polyR_.size();
    bool inside = false;
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        if ((polyY_[i] > y) != (polyY_[j] > y)) {
            const double xr = polyR_[j] + (y - polyY_[j]) * (polyR_[i] - polyR_[j]) / (polyY_[i] - polyY_[j]);
            if (r < xr) inside = !inside;
        }
    }
    if (!inside) return false;
    // Clearance from the wall edges (the axis and the open top need none).
    const std::size_t wallEdges = static_cast<std::size_t>(def_->cavityEnd);
    for (std::size_t i = 0; i < wallEdges && i + 1 < n; ++i) {
        const double ar = polyR_[i], ay = polyY_[i], dr = polyR_[i + 1] - ar, dy = polyY_[i + 1] - ay;
        const double l2 = dr * dr + dy * dy;
        const double f = l2 > 0 ? clamp(((r - ar) * dr + (y - ay) * dy) / l2, 0.0, 1.0) : 0.0;
        if (sq(r - ar - f * dr) + sq(y - ay - f * dy) < margin * margin) return false;
    }
    return true;
}

double VesselGeometry::cavityVolume(double h) const {
    // Clip the meridian polygon to y <= h (Sutherland-Hodgman), then Pappus:
    // V = 2 pi * integral of r dA = 2 pi * (1/6) sum (r_i + r_j)(r_i y_j - r_j y_i).
    std::vector<double> cr, cy;
    const std::size_t n = polyR_.size();
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        const bool inI = polyY_[i] <= h, inJ = polyY_[j] <= h;
        if (inI) {
            cr.push_back(polyR_[i]);
            cy.push_back(polyY_[i]);
        }
        if (inI != inJ) {
            const double f = (h - polyY_[i]) / (polyY_[j] - polyY_[i]);
            cr.push_back(polyR_[i] + f * (polyR_[j] - polyR_[i]));
            cy.push_back(h);
        }
    }
    double sum = 0;
    for (std::size_t i = 0; i < cr.size(); ++i) {
        const std::size_t j = (i + 1) % cr.size();
        sum += (cr[i] + cr[j]) * (cr[i] * cy[j] - cr[j] * cy[i]);
    }
    return 2 * kPi * std::abs(sum) / 6;
}

double VesselGeometry::fillHeight(double volume) const {
    double lo = *std::min_element(polyY_.begin(), polyY_.end()), hi = rimY_;
    if (volume >= cavityVolume(hi)) return hi;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        (cavityVolume(mid) < volume ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

// ---- Mass ---------------------------------------------------------------------

void VesselGeometry::massProperties() {
    const double rho = def_->glassDensity;
    const int n = std::max(400, static_cast<int>(s_.back() / 2.5e-4));
    const double ds = s_.back() / n;
    double M = 0, My = 0, Iyy = 0, Ixx = 0;
    for (int i = 0; i < n; ++i) {
        double r, y, nr, ny, t;
        at((i + 0.5) * ds, r, y, nr, ny, t);
        const double tx = ny;  // tangent = (ny, -nr)
        const double e = 0.5 * t * std::abs(nr) + 0.5 * ds * std::abs(tx);
        double dV, rc;
        wallElement(r, t, ds, e, dV, rc);
        const double dm = rho * dV;
        M += dm;
        My += dm * y;
        Iyy += dm * rc * rc;
        Ixx += dm * (0.5 * rc * rc + y * y);
    }
    mass_ = M;
    comY_ = M > 0 ? My / M : 0;
    const double Ic = Ixx - M * comY_ * comY_;
    inertia_ = Mat3::diag(Ic, Iyy, Ic);
}

void VesselGeometry::cellCentre(int cell, double& s, double& theta) const {
    const int i = cell / columns_, j = cell % columns_;
    s = 0.5 * (rowS_[i] + rowS_[i + 1]);
    theta = 2 * kPi * (j + 0.5) / columns_;
}

double VesselGeometry::cellVolume(int cell) const {
    const int i = cell / columns_;
    const double ds = rowS_[i + 1] - rowS_[i];
    double r, y, nr, ny, t;
    at(0.5 * (rowS_[i] + rowS_[i + 1]), r, y, nr, ny, t);
    const double e = 0.5 * t * std::abs(nr) + 0.5 * ds * std::abs(ny);
    double dV, rc;
    wallElement(r, t, ds, e, dV, rc);
    return dV / columns_;
}

Vec3 VesselGeometry::cellPoint(int cell) const {
    double s, th, r, y, nr, ny, t;
    cellCentre(cell, s, th);
    at(s, r, y, nr, ny, t);
    return {r * std::cos(th), y, r * std::sin(th)};
}

// ---- Contact samples ------------------------------------------------------------

void VesselGeometry::buildSamples() {
    samples_.clear();
    const int m = rows();
    const int stride = std::max(1, m / 48);
    const int ncol = 16;
    auto curvature = [&](int i) {
        // Meridional radius of curvature through rows i-1, i, i+1 (circumcircle).
        const int a = std::max(0, i - 8), b = std::min(m - 1, i + 8);
        if (a == i || b == i) return 1.0;
        const double x1 = rowR_[a], y1 = rowY_[a], x2 = rowR_[i], y2 = rowY_[i], x3 = rowR_[b], y3 = rowY_[b];
        const double d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
        if (std::abs(d) < 1e-12) return 1.0;
        const double a2 = std::hypot(x2 - x1, y2 - y1), b2 = std::hypot(x3 - x2, y3 - y2), c2 = std::hypot(x3 - x1, y3 - y1);
        return std::min(1.0, a2 * b2 * c2 / std::abs(d));
    };
    auto add = [&](int i, double sign, bool edge, double capDir) {
        const double r = rowR_[i], y = rowY_[i], nr = rowNr_[i], ny = rowNy_[i], t = rowT_[i];
        double pr = r + sign * nr * 0.5 * t, py = y + sign * ny * 0.5 * t;
        if (edge) {  // the rounded end of the wall, half a thickness beyond the last row
            pr = r + capDir * ny * 0.5 * t;
            py = y - capDir * nr * 0.5 * t;
        }
        pr = std::max(pr, 0.0);
        const double R1 = curvature(i);
        const double R2 = std::abs(nr) > 0.05 ? std::max(r, 1e-4) / std::abs(nr) : 1.0;
        const double Re = edge ? 0.5 * t : clamp(std::sqrt(R1 * R2), 0.5 * t, 0.5);
        const double Rs = clamp(std::min(R1, R2), t, 0.2);
        const double Rshell = clamp(std::max(R1, R2), t, 0.2);
        const int nc = pr < 1e-4 ? 1 : ncol;
        for (int j = 0; j < nc; ++j) {
            const double th = 2 * kPi * (j + 0.5 * (i % 2)) / ncol;
            WallSample w;
            w.p = {pr * std::cos(th), py, pr * std::sin(th)};
            w.s = rowS_[i];
            w.theta = th;
            w.radius = Re;
            w.reach = 0.76 * std::sqrt(Rs * t);
            w.t = t;
            w.shellRadius = Rshell;
            w.wallNormal = Vec3{nr * std::cos(th), ny, nr * std::sin(th)} * (sign < 0 ? -1.0 : 1.0);
            w.edge = edge;
            samples_.push_back(w);
        }
    };
    const int ce = static_cast<int>(std::lower_bound(rowS_.begin(), rowS_.end(), s_[def_->cavityEnd] - 1e-9) - rowS_.begin());
    for (int i = 0; i < m; i += stride) {
        add(i, +1, false, 0);
        if (i > ce) add(i, -1, false, 0);
    }
    add(m - 1, +1, false, 0);
    add(m - 1, -1, false, 0);
    add(0, 0, true, -1);
    add(m - 1, 0, true, +1);
}

// ---- Mesh -----------------------------------------------------------------------

VesselMesh VesselGeometry::buildMesh(const std::vector<char>& kept) const {
    VesselMesh out;
    const int m = rows(), C = columns_;
    std::vector<int> outerIdx(static_cast<std::size_t>(m) * C, -1), innerIdx(static_cast<std::size_t>(m) * C, -1);
    auto vert = [&](int i, int j, bool outer) {
        j = ((j % C) + C) % C;
        int& slot = (outer ? outerIdx : innerIdx)[static_cast<std::size_t>(i) * C + j];
        if (slot >= 0) return static_cast<std::uint32_t>(slot);
        const double sgn = outer ? 1 : -1;
        const double th = 2 * kPi * j / C, ct = std::cos(th), st = std::sin(th);
        const double r = std::max(0.0, rowR_[i] + sgn * rowNr_[i] * 0.5 * rowT_[i]);
        const double y = rowY_[i] + sgn * rowNy_[i] * 0.5 * rowT_[i];
        out.pos.push_back({r * ct, y, r * st});
        out.nrm.push_back(Vec3{rowNr_[i] * ct, rowNy_[i], rowNr_[i] * st} * sgn);
        slot = static_cast<int>(out.pos.size()) - 1;
        return static_cast<std::uint32_t>(slot);
    };
    auto keptCell = [&](int i, int j) {
        if (i < 0 || i >= m - 1) return false;
        return kept[cellOf(i, j)] != 0;
    };
    // A flat quad with its own vertices and a face normal pointing away from `inside`.
    auto side = [&](const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, const Vec3& inside) {
        Vec3 n = cross(b - a, d - a);
        if (length2(n) < 1e-24) n = cross(c - b, d - b);
        n = normalize(n);
        const Vec3 mid = (a + b + c + d) * 0.25;
        if (dot(n, mid - inside) < 0) n = -n;
        const std::uint32_t base = static_cast<std::uint32_t>(out.pos.size());
        for (const Vec3& v : {a, b, c, d}) {
            out.pos.push_back(v);
            out.nrm.push_back(n);
        }
        out.tri.insert(out.tri.end(), {base, base + 1, base + 2, base, base + 2, base + 3});
    };
    for (int i = 0; i + 1 < m; ++i)
        for (int j = 0; j < C; ++j) {
            if (!kept[cellOf(i, j)]) continue;
            const std::uint32_t o00 = vert(i, j, true), o10 = vert(i + 1, j, true), o11 = vert(i + 1, j + 1, true),
                                o01 = vert(i, j + 1, true);
            out.tri.insert(out.tri.end(), {o00, o11, o10, o00, o01, o11});
            const std::uint32_t i00 = vert(i, j, false), i10 = vert(i + 1, j, false), i11 = vert(i + 1, j + 1, false),
                                i01 = vert(i, j + 1, false);
            out.tri.insert(out.tri.end(), {i00, i10, i11, i00, i11, i01});
            const Vec3 centre = cellPoint(cellOf(i, j));
            auto P = [&](std::uint32_t k) { return out.pos[k]; };
            if (!keptCell(i - 1, j)) side(P(o00), P(o01), P(i01), P(i00), centre);
            if (!keptCell(i + 1, j)) side(P(o10), P(o11), P(i11), P(i10), centre);
            if (!keptCell(i, j - 1)) side(P(o00), P(o10), P(i10), P(i00), centre);
            if (!keptCell(i, j + 1)) side(P(o01), P(o11), P(i11), P(i01), centre);
        }
    return out;
}

// ---- Shards ---------------------------------------------------------------------

int ShardPattern::shardAt(const VesselGeometry& g, double s, double theta) const {
    if (cellShard.empty()) return 0;
    const int m = g.rows();
    int i = 0, lo = 0, hi = m - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        (g.rowS(mid) <= s ? lo : hi) = mid;
    }
    i = std::clamp(lo, 0, m - 2);
    double th = theta;
    while (th < 0) th += 2 * kPi;
    int j = static_cast<int>(th / (2 * kPi) * g.columns());
    j = std::clamp(j, 0, g.columns() - 1);
    return cellShard[g.cellOf(i, j)];
}

ShardPattern makeShards(const VesselGeometry& g, double s0, double theta0, int count, double spread, Rng& rng) {
    ShardPattern P;
    count = std::max(2, count);
    const double L = g.length();
    auto radiusAt = [&](double s) {
        double r, y, nr, ny, t;
        g.at(s, r, y, nr, ny, t);
        return std::max(r, 0.002);
    };
    for (int k = 0; k < count; ++k) {
        double s, th;
        if (k == 0) {
            s = s0;
            th = theta0;
        } else if (rng.uniform() < 0.75) {
            // Crowded round the impact: exponential distance in the developed wall.
            const double rho = -spread * std::log(std::max(rng.uniform(), 1e-9)), phi = 2 * kPi * rng.uniform();
            s = clamp(s0 + rho * std::sin(phi), 0.0, L);
            th = theta0 + rho * std::cos(phi) / radiusAt(s0);
        } else {
            // Anywhere, weighted by wall area (rejection on the radius).
            const double rmax = g.maxRadius();
            do {
                s = rng.uniform() * L;
            } while (rng.uniform() * rmax > radiusAt(s));
            th = 2 * kPi * rng.uniform();
        }
        P.seedS.push_back(s);
        P.seedTheta.push_back(wrapAngle(th));
    }
    std::vector<double> seedR(count);
    for (int k = 0; k < count; ++k) seedR[k] = radiusAt(P.seedS[k]);
    P.cellShard.assign(g.cells(), 0);
    std::vector<int> used(count, 0);
    for (int c = 0; c < g.cells(); ++c) {
        double s, th;
        g.cellCentre(c, s, th);
        const double rc = radiusAt(s);
        double best = std::numeric_limits<double>::max();
        int bk = 0;
        for (int k = 0; k < count; ++k) {
            const double dx = wrapAngle(th - P.seedTheta[k]) * 0.5 * (rc + seedR[k]), dy = s - P.seedS[k];
            const double d2 = dx * dx + dy * dy;
            if (d2 < best) {
                best = d2;
                bk = k;
            }
        }
        P.cellShard[c] = bk;
        used[bk] = 1;
    }
    // Renumber so that every shard id has cells.
    std::vector<int> remap(count, -1);
    int n = 0;
    for (int k = 0; k < count; ++k)
        if (used[k]) remap[k] = n++;
    for (auto& c : P.cellShard) c = remap[c];
    std::vector<double> ss, st;
    for (int k = 0; k < count; ++k)
        if (used[k]) {
            ss.push_back(P.seedS[k]);
            st.push_back(P.seedTheta[k]);
        }
    P.seedS = ss;
    P.seedTheta = st;
    P.count = n;
    return P;
}

}  // namespace lqd
