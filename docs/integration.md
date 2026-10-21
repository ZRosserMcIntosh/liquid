# Integration

## In a host

The core is a library (`lqd_core`) with no I/O of its own beyond loading the
JSON library from disk when asked. A host builds a `World`, adds its surfaces
and objects (or uses `buildScene` with a `SceneConfig`), advances it with its
own clock and reads the state back:

```cpp
lqd::World w(lib, lqd::budgetForTier(2), seed);
w.addPlane(lqd::Role::Floor, "floor", {0, 0, 0}, {0, 1, 0}, lib.surface("ceramic_tile"));
const int glass = w.addVessel(*lib.vessel("bordeaux_glass"), {0.4, 0.75, 0}, lqd::Quat{});
w.fillVessel(glass, *lib.liquid("red_wine"), 150e-6);
w.bodies[glass].applyImpulse({0.1, 0, 0}, w.bodies[glass].fromVessel({-0.045, 0.165, 0}));
w.advance(hostTime);   // whole internal steps: the result does not depend on the host's frame rate
```

Frames: the core works in glTF's frame (metres, right-handed, +y up). A host
converts at its boundary: OpenRoad Dynamics uses ISO 8855 (x forward, z up),
Aviation Dynamics FRD; a rotation and, for left-handed engines, one mirrored
axis. Static scenery the host already has can be given as planes, boxes and
capsules; moving scenery the host simulates itself can be mirrored by moving
static colliders each step (the liquid sees them; they do not feel it).

## Panes: the Glass Dynamics hand-off

When a round reaches a collider with `Role::Pane`, the world stops it there
and asks for the outcome. Three ways to answer:

1. `world.paneResolver = [&](const lqd::PaneStrike& s) { ... return result; };`
   answers synchronously (a native host running Glass Dynamics in-process);
2. `world.deferPanes = true;` parks the round (`Projectile::State::AtPane`)
   and lists the strike in `pendingStrikes()`; the host resolves it with
   `resolveStrike(i, result)` before the next `advance` (the browser lab does
   this with the Glass lab's WebAssembly module);
3. neither: the built-in estimate (18 J per mm of glass, a hole of 1.5 calibres).

`PaneStrike` is already in Glass Dynamics' pane frame: `x`, `y` from the
pane's lower left corner, `yaw` and `pitch` of the path from the normal,
`fromBack` for the −z face, the round's `speed`, and `glazing` (the pane's
Glass Dynamics glazing id). The round's Glass Dynamics impactor id is its
`ProjectileDef::glassImpactor`. With Glass Dynamics in the same process:

```cpp
gls::Pane pane(*glassLib.glazing(s.glazing), gls::budgetForTier(2), seed);
const int k = pane.fire(*glassLib.impactor(round.glassImpactor), {s.x, s.y}, s.speed, s.yaw, s.pitch, s.fromBack, 0.01);
while (pane.impactors()[k].phase != gls::ImpactorPhase::Done && pane.time() < 0.01) pane.advance(pane.time() + 2e-5);
const gls::ImpactReport& r = pane.impactors()[k].report;
lqd::PaneResult out;
out.penetrated = r.penetrated;
out.speedOut = r.speedOut;
out.holeRadius = std::sqrt(pane.status().openArea / lqd::kPi);
world.resolveStrike(0, out);
```

Keep advancing the pane on the liquid world's clock (offset by the strike
time) to grow its cracks; when a tempered pane comes down, disable its
collider (`colliders[id].enabled = false`) so the spatter flies through.

## Into Unreal and Unity

`lqdemo export '<scene json>' <seconds> <fps> <dir> [--ply]` (or the lab's
Export buttons) writes:

- `<scenario>.glb`: binary glTF 2.0. The liquid's surface is a flipbook: a node
  per recorded frame, each shown during its frame alone by a STEP animation of
  its scale, with the drops in the same mesh. Each glass body and shard is a
  mesh with its translation and rotation keyed every frame, hidden before it
  exists and after it breaks. The stains are flat ellipses on their surfaces
  that appear when they land. The static scene (floor, walls, table, pane,
  dummy) is under its own node, to keep or delete. One animation,
  `simulation`, holds every track. The files pass the Khronos glTF validator
  with no errors or warnings.
- `<scenario>_stains.json`: every stain with its centre `p`, surface normal
  `n`, tail direction `tail`, `width` and `length` in metres, impact angle
  `alpha_deg`, drop size, speed, time `t`, surface name and satellite flag, in
  the glTF frame. Use it to place decals procedurally (a decal box of
  length × width oriented by `n` and `tail`).
- `points_NNNNN.ply` (with `--ply`): one binary PLY per frame, position,
  `pscale` (radius) and `kind` (0 bulk, 1 drop), for Houdini, Blender or a
  point-cache importer (Niagara, VFX Graph) to remesh or re-render at film
  quality.

Unity: with glTFast (or UnityGLTF) installed, drop the .glb into Assets; enable
animation import in the importer settings and play the `simulation` clip.
Unreal: import the .glb with the Interchange glTF importer (import into the
level to keep the node hierarchy and its transform animation). Both importers
convert glTF's frame (metres, right-handed, +y up) to the engine's own; the
stains JSON is in glTF's frame, so convert it with the same mapping your
importer applies (or parent the decals under the imported glTF root).

Sizes: the flipbook stores a full surface mesh per frame; 2 s at 60 fps at
the medium tier is about 40-70 MB. Record at a lower fps, export with a
stride, or use the PLY caches and remesh in the DCC tool for long shots.

## In the browser

`web/build_wasi.sh` builds `liquid.wasm` with a C ABI (`web/lqd_web.cpp`) and
copies the loader `liquid.mjs`; `web/README.md` lists the calls. The lab at
`sao-paulo/liquids/` in the São Paulo site repository is the reference host:
it loads both engines, hands pane strikes to Glass Dynamics, draws the
liquid's surface mesh every frame and exports the recording.
