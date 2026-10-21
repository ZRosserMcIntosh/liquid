# Browser runtime

The fluid core compiles to WebAssembly with wasi-sdk. The build has no
exceptions or RTTI, and the WASI imports are stubs because the core does no
I/O. The C ABI is in `lqd_web.cpp`; one world at a time.

```sh
WASI_SDK=/path/to/wasi-sdk web/build_wasi.sh   # → web/dist/liquid.wasm, liquid.mjs, library.json
WASI_SDK=/path/to/wasi-sdk web/build_wasi.sh ../Sao-Paulo/sao-paulo/liquids   # straight into the site
node web/smoke_test.mjs                        # wine glass, gunshot with a host-resolved pane, export, bench
```

`liquid.mjs` exports `createLiquidModule(bytes?)`, which returns the raw
`lqd_*` exports plus these helpers (views point into the module's memory:
copy what you keep):

| Helper | Returns |
|---|---|
| `load(json)` | errors (strings) after loading a library document |
| `catalog()` | liquids, surfaces, vessels, projectiles, dummies, tiers |
| `build(config)` | builds a scene from a `SceneConfig` object (the JSON keys of `lqdemo`), returns `scene()` |
| `scene()` | the colliders (kind, role, pose, size, colour, surface, glazing), liquid, particle spacing, muzzle and aim |
| `state()` | 32 doubles: [0] time, [1] quiet, [2] particles, [3] drops, [4] stains, [5] bodies, [6] events, [7] steps, [8] spacing, [9] wetted m², [10] covered m², [11] poured m³, [12] bulk, [13] in the air, [14] stained, [15] absorbed, [16] lost, [17] breaks, [18] fastest liquid m/s, [19] rounds, [20] pending pane strikes, [21] sources, [22] next step s, [23] wounds, [24] breakups, [25] splashes, [26] colliders, [27] bled m³ |
| `particles()` | x y z speed per bulk particle |
| `drops()` | x y z diameter per drop |
| `stains(from)` | 16 per stain: pos, normal, tail (xyz each), width, length, impact angle, time, collider, satellite, source |
| `bodies()` | 12 per glass body: x y z, qx qy qz qw, shown, shard, mass, born, broken at |
| `bodyMesh(i)`, `surface(cell)` | `{positions, normals, indices}` of a body (own frame) or of the liquid's surface |
| `projectiles()`, `path(k)`, `wounds()` | rounds (pos, dir, speed, state, launch, corners), a round's path, wounds (pos, normal, exit, time) |
| `events(from)`, `sources()` | the engine log `{time, text}`, the spatter sources |
| `strikes()` | pending pane strikes in Glass Dynamics' pane frame, with the round's Glass Dynamics impactor id |
| `bpa(collider, source, readErrorMm, satellites)` | the stain analysis of one surface |
| `exportGlb(scene, stride)`, `exportPly()`, `exportStains()` | the recording as .glb bytes, this moment as .ply bytes, the stains as JSON |

Raw calls the helpers don't wrap: `lqd_advance(t, maxSteps)` (whole steps to
time t), `lqd_time()`, `lqd_record(fps, maxFrames)` and `lqd_capture()` (call
after each advance), `lqd_resolve_strike(i, penetrated, speedOut, holeRadius,
shattered)`, `lqd_estimate_strike(i)`, `lqd_collider_enable(id, on)`,
`lqd_hole(id, x, y, r)`, `lqd_shoot(projectile, ox, oy, oz, dx, dy, dz, speed,
delay)`, and `lqd_bench()` (ms per step of a 2,000-particle splash, for picking
a tier).
