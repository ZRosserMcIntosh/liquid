// Liquid Dynamics — minimal dependency-free test harness (same style as
// OpenRoad Dynamics' ord_test.hpp and Aviation Dynamics' avd_test.hpp and Glass Dynamics' lqd_test.hpp).
// SPDX-License-Identifier: MIT
#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace lqdtest {

struct Case {
    std::string name;
    std::function<void()> fn;
};

inline std::vector<Case>& registry() {
    static std::vector<Case> r;
    return r;
}
inline int& failures() {
    static int f = 0;
    return f;
}
struct Registrar {
    Registrar(const char* n, std::function<void()> f) { registry().push_back({n, std::move(f)}); }
};

// Long scenarios run at full length in optimised builds and shortened in
// Debug/sanitizer builds (same invariants, fewer simulated seconds).
inline constexpr bool kOptimisedBuild =
#ifdef NDEBUG
    true;
#else
    false;
#endif

}  // namespace lqdtest

#define LQD_CAT2(a, b) a##b
#define LQD_CAT(a, b) LQD_CAT2(a, b)
#define TEST_CASE(name)                                                                                 \
    static void LQD_CAT(lqd_test_fn_, __LINE__)();                                                      \
    static lqdtest::Registrar LQD_CAT(lqd_test_reg_, __LINE__)(name, &LQD_CAT(lqd_test_fn_, __LINE__)); \
    static void LQD_CAT(lqd_test_fn_, __LINE__)()

#define CHECK(cond)                                                                        \
    do {                                                                                   \
        if (!(cond)) {                                                                     \
            ++lqdtest::failures();                                                         \
            std::cerr << "  FAILED " << __FILE__ << ":" << __LINE__ << ": " #cond << "\n"; \
        }                                                                                  \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                                                  \
    do {                                                                                                       \
        const double lqd_a_ = (a), lqd_b_ = (b), lqd_t_ = (tol);                                               \
        if (!(std::abs(lqd_a_ - lqd_b_) <= lqd_t_)) {                                                          \
            ++lqdtest::failures();                                                                             \
            std::cerr << "  FAILED " << __FILE__ << ":" << __LINE__ << ": |" #a " - " #b "| <= " #tol << "  (" \
                      << lqd_a_ << " vs " << lqd_b_ << ")\n";                                                  \
        }                                                                                                      \
    } while (0)

#define CHECK_RANGE(v, lo, hi)                                                                                     \
    do {                                                                                                           \
        const double lqd_v_ = (v);                                                                                 \
        if (!(lqd_v_ >= (lo) && lqd_v_ <= (hi))) {                                                                 \
            ++lqdtest::failures();                                                                                 \
            std::cerr << "  FAILED " << __FILE__ << ":" << __LINE__ << ": " #v " = " << lqd_v_ << " not in [" << (lo) \
                      << ", " << (hi) << "]\n";                                                                    \
        }                                                                                                          \
    } while (0)

