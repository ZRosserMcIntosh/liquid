// Liquid Dynamics — Node smoke test of the WebAssembly build (web/dist).
// SPDX-License-Identifier: MIT
import { readFile } from 'node:fs/promises';
import createLiquidModule from './dist/liquid.mjs';

const fail = (m) => { console.error('FAIL ' + m); process.exitCode = 1; };
const L = await createLiquidModule();
const lib = await readFile(new URL('./dist/library.json', import.meta.url), 'utf8');
const errors = L.load(lib);
if (errors.length) fail('library: ' + errors.join('; '));
const cat = L.catalog();
console.log(`catalogue: ${cat.liquids.length} liquids, ${cat.surfaces.length} surfaces, ${cat.vessels.length} vessels, ${cat.projectiles.length} rounds, ${cat.impactors.length} impactors`);

// Wine glass knocked off the table onto tile.
let scene = L.build({ scenario: 'wine_glass', tier: 0 });
let t0 = performance.now();
for (let t = 0; t < 1.4; t += 1 / 30) { L.lqd_advance(t, 100000); }
let s = L.state();
console.log(`wine glass: ${s[2]} particles, ${s[3]} drops, ${s[4]} stains, ${s[17]} break(s), wet ${s[9].toFixed(3)} m2 in ${((performance.now() - t0) / 1000).toFixed(1)} s`);
if (s[17] !== 1) fail('the glass should break on tile');
if (Math.abs(s[12] + s[13] + s[14] + s[15] + s[16] - s[11]) > 1e-12) fail('volume not conserved');
const surf = L.surface();
if (surf.indices.length < 30) fail('no liquid surface mesh');
const b = L.bodies();
if (b.length < 12 * 5) fail('no shards');
const m = L.bodyMesh(1);
if (m.indices.length < 30) fail('no shard mesh');
const glb = L.exportGlb();
if (new DataView(glb.buffer).getUint32(0, true) !== 0x46546c67) fail('glb magic');

// Gunshot with a pane in front, resolved by the host.
scene = L.build({ scenario: 'gunshot', tier: 1, pane: 'front', defer_panes: true });
L.lqd_advance(0.03, 100000);
const strikes = L.strikes();
if (strikes.length !== 1) fail('expected one pane strike, got ' + strikes.length);
else {
  console.log(`pane strike at (${strikes[0].x.toFixed(2)}, ${strikes[0].y.toFixed(2)}) m, ${strikes[0].speed.toFixed(0)} m/s, glass impactor ${strikes[0].impactor}`);
  L.lqd_resolve_strike(0, 1, 330, 0.006, 0);
}
L.lqd_advance(2.0, 1000000);
s = L.state();
const ev = L.events(0).map((e) => e.text);
console.log(ev.filter((e) => /wound|pane|struck/.test(e)).join('\n'));
const wall = scene.colliders.find((c) => c.name === 'back wall');
const pane = scene.colliders.find((c) => c.role === 'pane');
const r = L.bpa(wall.id);
console.log('back wall: ' + r.summary);
if (!ev.some((e) => e.includes('exit wound'))) fail('no exit wound');
if (s[4] < 50) fail('too few stains');
// A car into the dummy: it is thrown, absorbs the blow, and the energy balances.
scene = L.build({ scenario: 'impact', tier: 0, impactor: 'car_50' });
L.lqd_advance(1.0, 1000000);
const imp = L.impactors();
const en = L.energy();
console.log(`car: ${(en.fromImpactors / 1e3).toFixed(1)} kJ into the dummy, ${(en.absorbed / 1e3).toFixed(1)} kJ absorbed, ${en.bursts} bursts; car at ${Math.hypot(imp[4], imp[5], imp[6]).toFixed(1)} m/s`);
if (imp.length !== 8) fail('expected one impactor');
if (en.fromImpactors < 5e3) fail('the car should hit the dummy');
if (Math.abs(en.start + en.given - en.kinetic - en.potential - en.absorbed) > 1e-6 * en.given) fail('energy book does not balance');
const bench = L.lqd_bench();
console.log(`bench: ${bench.toFixed(2)} ms per step`);
console.log(process.exitCode ? 'smoke test FAILED' : 'smoke test passed');
