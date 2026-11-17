// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Small vector, quaternion and 3x3 matrix types. The world frame is glTF's:
// metres, right-handed, +y up. Hosts convert at their boundary (see
// docs/integration.md); nothing inside the core depends on the host's axes.
#pragma once

#include <algorithm>
#include <cmath>

namespace lqd {

inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kGravity = 9.80665;  // m/s^2, standard gravity (published, CGPM 1901)

struct Vec3 {
    double x = 0, y = 0, z = 0;
    constexpr Vec3() = default;
    constexpr Vec3(double ax, double ay, double az) : x(ax), y(ay), z(az) {}
    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator-() const { return {-x, -y, -z}; }
    constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
    double operator[](int i) const { return i == 0 ? x : i == 1 ? y : z; }
    double& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }
};
inline constexpr Vec3 operator*(double s, const Vec3& v) { return v * s; }
inline constexpr double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double length(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline constexpr double length2(const Vec3& v) { return dot(v, v); }
inline Vec3 normalize(const Vec3& v, const Vec3& fallback = {0, 1, 0}) {
    const double l = length(v);
    return l > 1e-300 ? v / l : fallback;
}
// Any unit vector perpendicular to the unit vector n.
inline Vec3 perpendicular(const Vec3& n) {
    const Vec3 a = std::abs(n.x) < 0.8 ? Vec3{1, 0, 0} : Vec3{0, 0, 1};
    return normalize(cross(n, a));
}

struct Quat {
    double w = 1, x = 0, y = 0, z = 0;
    static Quat axisAngle(const Vec3& axis, double angle) {
        const Vec3 a = normalize(axis);
        const double s = std::sin(0.5 * angle);
        return {std::cos(0.5 * angle), a.x * s, a.y * s, a.z * s};
    }
    Quat operator*(const Quat& o) const {
        return {w * o.w - x * o.x - y * o.y - z * o.z, w * o.x + x * o.w + y * o.z - z * o.y,
                w * o.y - x * o.z + y * o.w + z * o.x, w * o.z + x * o.y - y * o.x + z * o.w};
    }
    Quat conj() const { return {w, -x, -y, -z}; }
    Quat normalized() const {
        const double n = std::sqrt(w * w + x * x + y * y + z * z);
        return n > 0 ? Quat{w / n, x / n, y / n, z / n} : Quat{};
    }
    Vec3 rotate(const Vec3& v) const {
        const Vec3 q{x, y, z};
        const Vec3 t = 2.0 * cross(q, v);
        return v + w * t + cross(q, t);
    }
    Vec3 unrotate(const Vec3& v) const { return conj().rotate(v); }
    // Advances the orientation by angular velocity omega (world frame) over dt.
    Quat integrated(const Vec3& omega, double dt) const {
        const double a = length(omega) * dt;
        if (a < 1e-15) return *this;
        return (axisAngle(omega, a) * *this).normalized();
    }
};

struct Mat3 {
    double m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    static Mat3 identity() {
        Mat3 r;
        r.m[0][0] = r.m[1][1] = r.m[2][2] = 1;
        return r;
    }
    static Mat3 diag(double a, double b, double c) {
        Mat3 r;
        r.m[0][0] = a;
        r.m[1][1] = b;
        r.m[2][2] = c;
        return r;
    }
    static Mat3 fromQuat(const Quat& q) {
        Mat3 r;
        const Vec3 c0 = q.rotate({1, 0, 0}), c1 = q.rotate({0, 1, 0}), c2 = q.rotate({0, 0, 1});
        for (int i = 0; i < 3; ++i) {
            r.m[i][0] = c0[i];
            r.m[i][1] = c1[i];
            r.m[i][2] = c2[i];
        }
        return r;
    }
    Vec3 operator*(const Vec3& v) const {
        return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z, m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
    }
    Mat3 operator*(const Mat3& o) const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) r.m[i][j] += m[i][k] * o.m[k][j];
        return r;
    }
    Mat3 operator+(const Mat3& o) const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[i][j] + o.m[i][j];
        return r;
    }
    Mat3 operator-(const Mat3& o) const { return *this + o * -1.0; }
    Mat3 operator*(double s) const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[i][j] * s;
        return r;
    }
    Mat3 transposed() const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[j][i];
        return r;
    }
    double det() const {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }
    Mat3 inverse() const {
        const double d = det();
        Mat3 r;
        if (std::abs(d) < 1e-300) return r;
        const double i = 1.0 / d;
        r.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * i;
        r.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * i;
        r.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * i;
        r.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * i;
        r.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * i;
        r.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * i;
        r.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * i;
        r.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * i;
        r.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * i;
        return r;
    }
};
// Outer product a b^T.
inline Mat3 outer(const Vec3& a, const Vec3& b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r.m[i][j] = a[i] * b[j];
    return r;
}

inline constexpr double clamp(double v, double lo, double hi) { return v < lo ? lo : v > hi ? hi : v; }
inline constexpr double sq(double v) { return v * v; }

}  // namespace lqd
