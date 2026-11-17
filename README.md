# Liquid Dynamics

A physically grounded fluid engine in C++20 for the São Paulo open world. It
spills, splashes and sprays wine, water, oil and (green) stage blood: a wine
glass knocked off a table, a bottle bursting on the floor, a crash dummy shot
in front of or behind a pane of glass, struck by a car or a wrecking ball, or
dropped from a height. It hands back the liquid, the drops,
the stains and the broken glass, and reads the stains the way a bloodstain
analyst does. It has no renderer and no game-engine dependency; liquids,
surfaces, vessels, rounds, dummies and impactors are JSON with a provenance note on every
value. The same core runs natively and in the browser as WebAssembly, where it
drives the lab at
[zrossermcintosh.com/sao-paulo/liquids](https://zrossermcintosh.com/sao-paulo/liquids/),
and exports to glTF for Unreal and Unity.

It is the liquid partner of the project's other engines: **Glass Dynamics**
breaks the panes the rounds go through (the lab runs both side by side),
**OpenRoad Dynamics** and **Aviation Dynamics** are the vehicles and aircraft.
[docs/integration.md](docs/integration.md) shows the couplings and the export
pipeline.

## What it does

| Area | What is implemented | Evidence |
|---|---|---|
| Bulk liquid | Smoothed-particle hydrodynamics kept incompressible by position-based density constraints (Macklin & Müller 2013), cubic-spline kernel, spatial-hash neighbours. Walls are analytic: the kernel volume a half-space fills, combined by inclusion-exclusion at corners, so floors, tables, glass walls and shards need no wall particles | `tests/test_scenes.cpp` (a resting tank: level surface, <3% compression; Martin & Moyce's dam break) |
| Liquid properties | The liquid's real viscosity as the SPH Laplacian (Morris et al. 1997), shear-thinning where it is (Carreau-Yasuda on the local strain rate: whole blood after Cho & Kensey 1991); surface tension as the continuum surface force with the real σ (Brackbill et al. 1992; Morris 2000); contact-angle hysteresis per surface, which stops a spill spreading | `tests/test_foundations.cpp` |
| Two-way coupling | Every push the liquid takes from a moving glass wall goes back to the glass body, relaxed by the added mass of the liquid in contact; pieces below the liquid's resolution blend to one-way | `tests/test_scenes.cpp` |
| Drops in air | Spheres with the standard drag curve (Schiller & Naumann; Newton above Re 1000), semi-implicit so the finest mist is stable; aerodynamic breakup above We 12 after the Pilch & Erdman (1987) time. Drops at their terminal speed land at once where they would land. Liquid thrown off faster than 2 m/s relative to the rest tears into millimetre drops; slower blobs rejoin the bulk; a surface holds a thin film where drops land and what lands beyond it pools back into liquid (a puddle, a run), volume exact | Terminal speeds against Gunn & Kinzer (1949), `tests/test_foundations.cpp` |
| Stains | Width from the maximum spreading of Laan et al. (2014) (fitted on blood and water), length from the impact angle (W/L = sin α), tail along the travel, spines and satellite drops above Mundo et al.'s (1995) splash threshold, lowered on rough surfaces; absorbent surfaces wick | `tests/test_foundations.cpp`, the drip and impact-angle cards in `tests/test_scenes.cpp` |
| Glassware | A wine glass and a 750 ml bottle as revolved wall profiles: exact signed distance, cavity and fill level, mass and inertia by Pappus, a render mesh with real wall thickness | Masses and capacities against stemware and bottle norms, `tests/test_foundations.cpp` |
| Breaking glass | A strike's peak force from Hertz contact (or a carpet's pile as a Winkler foundation) in series with the wall's shell stiffness, its bending and membrane stress, and a Weibull probability of failure; broken glass becomes Voronoi shards crowded round the strike, flown as rigid bodies while the liquid escapes between them | `tests/test_foundations.cpp`, `tests/test_scenes.cpp` (a knock holds, tile breaks it, carpet saves it; a bottle from 1 m onto concrete) |
| Ragdoll | Optional (off by default): the dummy as a 78 kg ragdoll (Hybrid III mass): capsule ends as point masses; a flexing spine and neck, shoulders and hips within human ranges of motion, knees and elbows as one-way hinges, rocking ankles, limbs that collide with the body; position-based dynamics with friction against the scene, energy-bounded; hits push it at the point they strike (real momentum, or a fraction of the energy deposited: `knockback`); it stands on a release stand until hit and falls with the momentum the rounds leave; wounds, drips and stains on it move with it | `tests/test_scenes.cpp` |
| Rounds and the dummy | A capsule crash dummy; rounds cross it with Poncelet's resistance (fitted to a 9 mm FMJ penetrating 0.65 m of gelatin), throw back spatter from the entry wound and forward spatter from the exit, and the wounds drip (Tate's law). .50 BMG among the rounds. Shot loads (12-gauge 00 buck: nine pellets, cylinder-bore spread, shot string) fire each pellet as its own round; spatter too fine for the drop budget is sampled, each drop standing for several (super-droplets), so volumes stay exact. A pane in front, behind or both is resolved by Glass Dynamics through a hand-off, or by a built-in estimate | `tests/test_scenes.cpp` |
| Impacts and energy absorption | Heavy blunt impactors with their real mass (a sledgehammer, a wrecking ball, a car at 50 km/h, a truck at 80 km/h) hit the ragdoll through inelastic, momentum-conserving contacts and brake; stunt falls from a height. The dummy keeps an energy book (given, absorbed in the blow and afterwards, still moving it) that balances to rounding; a segment struck fast and hard enough bursts its skin, scaled by `spatter_boost`. A car at 50 km/h hands it the inelastic M m v²/(M+m) and throws it 15 m | `tests/test_scenes.cpp`, `docs/validation.md` |
| Glass in pieces | Glass Dynamics hands over a broken pane's pieces (its fragment map, outlines and poses): the pane opens where its pieces have gone, loose pieces become slab colliders, and spatter lands on the pieces and rides with them (back spatter on the inside of a windscreen, forward spatter on a window falling out behind); the glTF export animates every piece with its stains | `tests/test_scenes.cpp`, `tests/test_foundations.cpp` |
| Stain analysis | Bloodstain pattern analysis on any surface: impact angle from shape, directionality, area of convergence by least squares, area of origin by the tangent method, scored against the true source | `tests/test_foundations.cpp` (exact on straight lines), `tests/test_scenes.cpp` (a pane behind the dummy points back to the exit wound) |
| Export | Binary glTF 2.0: the liquid's surface as a flipbook, drops, every glass body and shard animated, the pieces of broken panes animated with their stains, stains as decals, the static scene. PLY point caches, stains as JSON. Files pass the Khronos glTF validator with no errors or warnings | `tests/test_scenes.cpp` |
| Any device | Five tiers (potato to ultra) change the particle size (8 to 3 mm) and the drop budget, not the physics. Steps are chosen from the state alone, so results do not depend on how the host chops its frames | `tests/test_scenes.cpp` (same seed, same stains at 60 and 240 Hz) |
| Browser | WebAssembly build (wasi-sdk, no exceptions or RTTI, ~450 KB) with a C ABI; Node smoke test | `web/`, `web/smoke_test.mjs` |

45 test cases, on Linux (gcc, clang, sanitizers), macOS, Windows and WASM in
CI. How each scenario compares with published data is in
[docs/validation.md](docs/validation.md) (regenerated by `lqdemo validate`);
how far each model goes, and its main gap, in [docs/physics.md](docs/physics.md).

## Build and run

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/tests/lqd_tests                 # every test case
./build/tests/lqd_tests "gunshot"       # a filter
./build/lqdemo validate                 # the comparison table with published data
./build/lqdemo run '{"scenario":"wine_glass","floor":"oak_floor"}' 2
./build/lqdemo run '{"scenario":"gunshot","pane":"behind","glazing":"tempered_6"}' 3
./build/lqdemo run '{"scenario":"impact","impactor":"car_50","spatter_boost":3}' 4   # prints the dummy's energy book
./build/lqdemo run '{"scenario":"fall","fall_m":10}' 3
./build/lqdemo export '{"scenario":"bottle_drop","tier":3}' 2.0 60 out/ --ply   # glTF + stains + PLY for a game engine
./build/lqdemo list
```

Options: `LQD_BUILD_TESTS`, `LQD_BUILD_TOOLS`, `LQD_WARNINGS_AS_ERRORS`.

The browser build needs [wasi-sdk](https://github.com/WebAssembly/wasi-sdk);
version 25 is tested.

```sh
WASI_SDK=/path/to/wasi-sdk web/build_wasi.sh   # web/dist: liquid.wasm, liquid.mjs, library.json
node web/smoke_test.mjs
```

## Using it

```cpp
#include "lqd/scenarios.hpp"
#include "lqd/export.hpp"

lqd::Library lib;
lqd::loadDefaultLibrary(lib);                         // data/*.json

lqd::SceneConfig c;                                   // or parseSceneConfig(json, c, &err)
c.scenario = "wine_glass";
c.floor = "ceramic_tile";
lqd::World w(lib, lqd::budgetForTier(2), /*seed*/ 1);
lqd::buildScene(w, c, &err);

lqd::Recorder rec(60);
while (w.time() < 2.0) { w.advance(w.time() + 1.0 / 60); rec.capture(w); }   // whole internal steps

for (const auto& p : w.fluid.x) { /* bulk particles */ }
for (const auto& d : w.spray.drops) { /* x, v, diameter */ }
for (const auto& s : w.spray.stains) { /* pos, normal, tail dir, width, length, impact angle */ }
for (const auto& b : w.bodies) { /* the glass and its shards: x, q; lqd::bodyMesh(w, i) */ }
const auto glb = lqd::exportGlb(w, rec);             // for Unreal / Unity / Blender

lqd::BpaOptions o; o.collider = w.floorCollider();
const lqd::BpaResult r = lqd::analyseStains(w.spray, w.colliders, o);   // area of convergence and origin
```

## Repository layout

```
src/core/      lqd_core: the engine-independent fluid core (namespace lqd)
tools/lqdemo/  headless runs, exports for game engines, the validation table
tests/         dependency-free test harness and suites
data/          liquids, surfaces, vessels, projectiles, dummies and impactors (JSON, with provenance)
web/           C ABI, WASM build, JS loader, smoke test
docs/          physics, validation, integration
```

## Principles

Shared with the vehicle, flight and glass cores:

- **Nothing canned.** Every splash, drop, stain and shard comes out of the
  mechanics; the renderer only draws.
- **Lower tiers lower resolution, not physics.** A weak device gets bigger
  particles and fewer drops, and the clock slows down rather than skipping physics.
- **Scientific honesty.** Every value is tagged `published`, `calibrated`,
  `estimated` or `gameplay` with its source. The calibration constants are
  listed with what they were fitted to, and the model's known biases (the
  tangent method's overshoot on a floor, the coarse films below the particle
  size) are reported, not tuned away.
- **G-rated.** The dummy is a crash-test target, not an injury model, and its
  stage blood is green.

## Licence

MIT. Forensic method names (impact angle, area of convergence, area of origin,
tangent method) describe the published bloodstain pattern analysis practice
the analysis follows; nothing here is evidence or a substitute for a trained
analyst.
