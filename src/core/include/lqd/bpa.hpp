// Liquid Dynamics — engine-independent fluid core.
// SPDX-License-Identifier: MIT
//
// Bloodstain pattern analysis on the stains of one surface, done the way an
// analyst does it at a scene, then scored against the truth the simulation
// knows:
//   1. impact angle of each elongated stain from its shape, alpha = asin(W/L);
//   2. directionality from the tail;
//   3. area of convergence: the point in the surface's plane closest, in least
//      squares, to the lines drawn back through every stain;
//   4. area of origin by the tangent method: the height above the surface of
//      each stain's line is its distance to the convergence times tan(alpha).
// Straight "strings" ignore gravity and drag, so on a floor they put the
// origin too high; the report says by how much (Carter 2001; Illes & Boue 2011
// quantify the same bias at real scenes).
#pragma once

#include <string>
#include <vector>

#include "lqd/collider.hpp"
#include "lqd/math.hpp"
#include "lqd/spray.hpp"

namespace lqd {

struct BpaOptions {
    int collider = -1;   // the surface to analyse
    int source = -1;     // only stains from this source (-1: all)
    double minLength = 1.0e-3;   // m, stains smaller than this are too small to measure
    double maxRatio = 0.94;      // W/L above this: too round to give a direction
    double minAlpha = 0.17;      // rad (10 deg): grazing stains are excluded
    bool satellites = false;     // include satellite stains
    // Measurement error of a ruler at the scene: W and L read to this (m).
    double readError = 0;
};

struct BpaLine {
    int stain = -1;
    Vec3 pos;       // the stain
    Vec3 back;      // unit, in the surface, back towards the convergence
    double ratio = 0, alpha = 0, alphaTrue = 0;  // W/L, estimated and true impact angle
    double distance = 0, height = 0;             // to the convergence, tangent-method height
};

struct BpaResult {
    bool ok = false;
    int considered = 0, used = 0;
    Vec3 planePoint, planeNormal;
    Vec3 convergence;          // world point in the surface
    double convergenceRms = 0; // rms distance of the lines from it, m
    double height = 0;         // median tangent-method height above the surface, m
    double heightLow = 0, heightHigh = 0;  // interquartile range
    Vec3 origin;               // convergence + height along the normal
    double meanAngleError = 0; // mean |alpha_est - alpha_true|, rad
    bool truth = false;
    Vec3 trueOrigin;
    double convergenceError = 0;  // m, in the surface plane
    double heightError = 0;       // m, estimated minus true
    std::vector<BpaLine> lines;
    std::string summary;
};

BpaResult analyseStains(const Spray& spray, const std::vector<Collider>& colliders, const BpaOptions& options);

}  // namespace lqd
