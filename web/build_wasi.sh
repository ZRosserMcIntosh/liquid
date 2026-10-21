#!/bin/sh
# Liquid Dynamics — build the browser lab's WebAssembly module with wasi-sdk.
# SPDX-License-Identifier: MIT
# Usage: WASI_SDK=/path/to/wasi-sdk web/build_wasi.sh [output-dir]   (default: web/dist)
# Tested with wasi-sdk 25. The core builds without exceptions or RTTI here.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/.." && pwd)
OUT=${1:-"$HERE/dist"}
: "${WASI_SDK:?set WASI_SDK to the wasi-sdk directory}"
mkdir -p "$OUT"
EXPORTS=""
for f in lqd_text lqd_text_length lqd_load lqd_resolve lqd_error_text lqd_catalog lqd_build lqd_scene \
         lqd_advance lqd_time lqd_next_step lqd_state lqd_particles lqd_drops lqd_stains lqd_bodies lqd_body_mesh \
         lqd_surface lqd_events lqd_projectiles lqd_projectile_path lqd_shoot lqd_wounds lqd_sources lqd_strikes \
         lqd_resolve_strike lqd_estimate_strike lqd_collider_enable lqd_hole lqd_bpa lqd_record lqd_capture lqd_blob \
         lqd_export_glb lqd_export_ply lqd_export_stains lqd_bench; do
  EXPORTS="$EXPORTS -Wl,--export=$f"
done
"$WASI_SDK/bin/clang++" --target=wasm32-wasip1 --sysroot="$WASI_SDK/share/wasi-sysroot" \
  -std=c++20 -O2 -fno-exceptions -fno-rtti -Wall -Wextra -mexec-model=reactor \
  -I"$REPO/src/core/include" \
  "$HERE/lqd_web.cpp" "$REPO"/src/core/src/*.cpp \
  $EXPORTS -Wl,--initial-memory=67108864 -Wl,--max-memory=2147483648 -Wl,--strip-all \
  -o "$OUT/liquid.wasm"
cp "$HERE/lqd_loader.mjs" "$OUT/liquid.mjs"
# One library file for the page.
node -e '
const fs = require("fs"), d = process.argv[1];
const out = { liquids: [], surfaces: [], vessels: [], projectiles: [], dummies: [] };
for (const f of ["liquids", "surfaces", "vessels", "projectiles", "dummies"]) {
  const j = JSON.parse(fs.readFileSync(`${d}/${f}.json`, "utf8"));
  for (const k of Object.keys(out)) if (j[k]) out[k].push(...j[k]);
}
fs.writeFileSync(process.argv[2], JSON.stringify(out));
' "$REPO/data" "$OUT/library.json"
ls -l "$OUT"
