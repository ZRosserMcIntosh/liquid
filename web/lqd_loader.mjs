// Liquid Dynamics — loader for the wasi-sdk build (web/build_wasi.sh).
// SPDX-License-Identifier: MIT
// `const L = await createLiquidModule(bytes?)` then L.lqd_*() and helpers. The
// core does no I/O in the browser, so the WASI imports are minimal stubs.
export default async function createLiquidModule(wasmBytes) {
  let bytes = wasmBytes;
  if (!bytes) {
    const url = new URL('./liquid.wasm', import.meta.url);
    if (typeof process !== 'undefined' && process.versions?.node && url.protocol === 'file:') {
      const { readFile } = await import('node:fs/promises');
      bytes = await readFile(url);
    } else {
      bytes = await (await fetch(url)).arrayBuffer();
    }
  }
  let memory = null;
  const ENOSYS = 52;
  const wasi = new Proxy({
    proc_exit: (c) => { throw new Error('liquid.wasm exited with code ' + c); },
    fd_write: (fd, iovs, n, nwritten) => {
      const v = new DataView(memory.buffer); let s = '', total = 0;
      for (let i = 0; i < n; ++i) { const p = v.getUint32(iovs + i * 8, true), l = v.getUint32(iovs + i * 8 + 4, true); s += new TextDecoder().decode(new Uint8Array(memory.buffer, p, l)); total += l; }
      v.setUint32(nwritten, total, true); if (s.trim()) (fd === 2 ? console.warn : console.log)(s.trimEnd()); return 0;
    },
    environ_sizes_get: (c, s) => { const v = new DataView(memory.buffer); v.setUint32(c, 0, true); v.setUint32(s, 0, true); return 0; },
    environ_get: () => 0,
    clock_time_get: (id, prec, out) => { new DataView(memory.buffer).setBigUint64(out, BigInt(Math.round(performance.now() * 1e6)), true); return 0; },
    random_get: (p, n) => { crypto.getRandomValues(new Uint8Array(memory.buffer, p, n)); return 0; },
  }, { get: (t, k) => t[k] ?? (() => ENOSYS) });
  const { instance } = await WebAssembly.instantiate(bytes, { wasi_snapshot_preview1: wasi });
  const ex = instance.exports;
  memory = ex.memory;
  ex._initialize?.();
  const L = { ...Object.fromEntries(Object.entries(ex).filter(([k]) => k.startsWith('lqd_'))) };
  const enc = new TextEncoder(), dec = new TextDecoder();
  L.text = (len) => dec.decode(new Uint8Array(memory.buffer, ex.lqd_text(0), len ?? ex.lqd_text_length()));
  const put = (s) => { const b = enc.encode(s); new Uint8Array(memory.buffer, ex.lqd_text(b.length + 1), b.length).set(b); return b.length; };
  // Loads a JSON library document; returns the error strings (empty when fine).
  L.load = (json) => {
    let n = ex.lqd_load(put(json));
    n = ex.lqd_resolve();
    const errors = [];
    for (let i = 0; i < n; ++i) errors.push(L.text(ex.lqd_error_text(i)));
    return errors;
  };
  L.catalog = () => JSON.parse(L.text(ex.lqd_catalog()));
  // Builds a scene from a config object; throws with the engine's message.
  L.build = (config) => { if (!ex.lqd_build(put(JSON.stringify(config)))) throw new Error(L.text()); return L.scene(); };
  L.scene = () => JSON.parse(L.text(ex.lqd_scene()));
  L.state = () => new Float64Array(memory.buffer, ex.lqd_state(), 32).slice();
  // Views below point into the module's memory: copy what you keep.
  const counted = (p, stride) => { const n = new Float32Array(memory.buffer, p, 1)[0]; return new Float32Array(memory.buffer, p + 4, n * stride); };
  L.particles = () => counted(ex.lqd_particles(), 4);
  L.drops = () => counted(ex.lqd_drops(), 8);
  L.dummy = () => counted(ex.lqd_dummy(), 7);
  L.dummyStains = () => counted(ex.lqd_dummy_stains(), 10);
  L.impactors = () => counted(ex.lqd_impactors(), 8);
  L.energy = () => {
    const e = new Float64Array(memory.buffer, ex.lqd_energy(), 13).slice();
    return { start: e[0], fromHits: e[1], fromImpactors: e[2], inBlows: e[3], afterwards: e[4], inTissue: e[5],
             kinetic: e[6], potential: e[7], peakKinetic: e[8], bursts: e[9], burstVolume: e[10], given: e[11], absorbed: e[12] };
  };
  L.stains = (from = 0) => counted(ex.lqd_stains(from), 16);
  L.movedStains = () => counted(ex.lqd_moved_stains(), 11);
  // A broken pane from Glass Dynamics (see docs/integration.md): its fragment
  // map ({nx, ny, label} for the front ply, and optionally the back ply's), its
  // pieces (G.fragmentShapes records) and every piece's pose (G.fragmentPoses).
  const input = (n) => ex.lqd_input(n);
  L.paneMap = (pane, cell, front, back = null, laminated = false) => {
    const n = front.nx * front.ny;
    new Int32Array(memory.buffer, input(2 * n), 2 * n).set(back ? [...front.label, ...back.label] : front.label);
    ex.lqd_pane_map(pane, front.nx, front.ny, cell, back ? 1 : 0, laminated ? 1 : 0);
  };
  L.panePiece = (pane, shape) => {
    const n = shape.outline.length / 2, holes = shape.holes || [];
    const data = [...shape.outline, holes.length];
    for (const h of holes) data.push(h.length / 2, ...h);
    new Float32Array(memory.buffer, input(data.length), data.length).set(data);
    ex.lqd_pane_piece(pane, shape.index, n, shape.thickness, shape.chip ? 1 : 0);
  };
  L.panePoses = (pane, poses, n) => {
    new Float32Array(memory.buffer, input(Math.max(12 * n, 1)), 12 * n).set(poses.subarray(0, 12 * n));
    ex.lqd_pane_poses(pane, n);
  };
  L.bodies = () => counted(ex.lqd_bodies(), 12);
  L.wounds = () => counted(ex.lqd_wounds(), 8);
  L.projectiles = () => counted(ex.lqd_projectiles(), 10);
  L.path = (k) => counted(ex.lqd_projectile_path(k), 3);
  const mesh = (p) => {
    const h = new Float32Array(memory.buffer, p, 2), nv = h[0], nt = h[1];
    const f = new Float32Array(memory.buffer, p + 8, nv * 6 + nt * 3);
    return { positions: f.slice(0, nv * 3), normals: f.slice(nv * 3, nv * 6), indices: Uint32Array.from(f.subarray(nv * 6)) };
  };
  L.bodyMesh = (i) => mesh(ex.lqd_body_mesh(i));
  L.surface = (cellFactor = 0) => mesh(ex.lqd_surface(cellFactor));
  L.events = (from) => { const t = L.text(ex.lqd_events(from)); return t ? t.split('\n').filter(Boolean).map((l) => { const i = l.indexOf('\t'); return { time: +l.slice(0, i), text: l.slice(i + 1) }; }) : []; };
  L.sources = () => JSON.parse(L.text(ex.lqd_sources()));
  L.strikes = () => JSON.parse(L.text(ex.lqd_strikes()));
  L.bpa = (collider, source = -1, readErrorMm = 0.25, satellites = false) => JSON.parse(L.text(ex.lqd_bpa(collider, source, readErrorMm, satellites ? 1 : 0)));
  const blob = (n) => new Uint8Array(memory.buffer, ex.lqd_blob(), n).slice();
  L.exportGlb = (scene = true, stride = 1) => blob(ex.lqd_export_glb(scene ? 1 : 0, stride));
  L.exportPly = () => blob(ex.lqd_export_ply());
  L.exportStains = () => L.text(ex.lqd_export_stains());
  return L;
}
