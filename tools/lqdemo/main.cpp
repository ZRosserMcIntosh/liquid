// Liquid Dynamics — lqdemo: headless runs, exports and the validation table.
// SPDX-License-Identifier: MIT
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "lqd/bpa.hpp"
#include "lqd/export.hpp"
#include "lqd/library.hpp"
#include "lqd/scenarios.hpp"
#include "lqd/world.hpp"

using namespace lqd;

int main(int argc, char** argv) {
    Library lib;
    std::vector<std::string> errors;
    if (!loadDefaultLibrary(lib, &errors)) {
        for (const auto& e : errors) std::cerr << e << "\n";
        return 1;
    }
    const std::string cmd = argc > 1 ? argv[1] : "help";
    if (cmd == "run") {
        std::string json = argc > 2 ? argv[2] : "{}";
        double until = argc > 3 ? std::atof(argv[3]) : 2.0;
        SceneConfig c;
        std::string err;
        if (!parseSceneConfig(json, c, &err)) {
            std::cerr << err << "\n";
            return 1;
        }
        World w(lib, budgetForTier(c.tier), c.seed);
        if (!buildScene(w, c, &err)) {
            std::cerr << err << "\n";
            return 1;
        }
        const auto t0 = std::chrono::steady_clock::now();
        std::size_t shown = 0;
        for (double t = 0; t < until; t += 0.05) {
            w.advance(t + 0.05);
            for (; shown < w.events.size(); ++shown) std::printf("%8.4f s  %s\n", w.events[shown].time, w.events[shown].text.c_str());
            const VolumeBook v = w.volumes();
            std::printf("%8.3f s  steps %6d  bulk %5zu (%.1f ml)  drops %5zu  stains %5zu  wet %.3f m2  err %.3f  vmax %.2f  dt %.1e\n", w.time(),
                        w.steps(), w.fluid.size(), v.bulk * 1e6, w.spray.drops.size(), w.spray.stains.size(),
                        w.wettedArea(), w.fluid.maxCompression(), w.fluid.maxSpeed(), w.nextStep());
            if (w.quiet() && t > 0.2) break;
        }
        for (; shown < w.events.size(); ++shown) std::printf("%8.4f s  %s\n", w.events[shown].time, w.events[shown].text.c_str());
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const VolumeBook v = w.volumes();
        std::printf("volume: in %.3f ml = bulk %.3f + air %.3f + stains %.3f + absorbed %.3f + lost %.3f (%.3f)\n",
                    v.initial * 1e6, v.bulk * 1e6, v.air * 1e6, v.stained * 1e6, v.absorbed * 1e6, v.lost * 1e6,
                    v.total() * 1e6);
        std::printf("%.2f s of simulation in %.2f s\n", w.time(), secs);
        for (const auto& col : w.colliders) {
            int n = 0;
            for (const auto& st : w.spray.stains) n += st.collider == col.id ? 1 : 0;
            if (n < 5) continue;
            BpaOptions o;
            o.collider = col.id;
            const BpaResult r = analyseStains(w.spray, w.colliders, o);
            std::printf("stain analysis, %s (%d stains): %s\n", col.name.c_str(), n, r.summary.c_str());
        }
        return 0;
    }
    std::printf("lqdemo run '<scene json>' [seconds]\n");
    return 0;
}
