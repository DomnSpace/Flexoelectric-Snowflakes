# PASS 7 — Versioned physical state and five-sheet optical preview

Pass 7 is deliberately **not** a new growth law or a solved spatial-height field. It freezes the existing 2.5-D vapor/heat/ice growth loop as a comparison baseline and adds a correctly identified view of the one currently solved **uniform basal thickness**.

## Ownership and provenance

- `app/solver/hex25d.js`, `growth25d.js` continue to own `phi`, `sigma`, `dT`, ice masks, total deposited mass, time and the single `basalHalfThickness` scalar. This PR changes no physics equations, numerical timesteps, seed rules, grid limits or growth boundaries.
- `app/solver/ice-snapshot-v1.js` owns the versioned replay boundary. It verifies the shape of full-grid arrays, units, phase contour and thickness; reports unavailable fields explicitly.
- `app/render/uniform-five-layers.js` owns a pure five-slice plan plus a Three.js adapter. It cannot advance simulation time, change phase fractions or transfer ice mass.
- `app/observatory.js` orchestrates render mode, snapshots and replay; it must not become a second physical solver.

Snapshot schema:
```js
{
  snapshotSchema: {
    name: "snowflake-observatory/snapshot",
    version: 1,
    topology: "triangular" | "square",
    unitSystem: "SI",
    fields: {
      phase: "phi",
      vapor: "u",
      temperature: "dT",
      frozenMask: "ice",
      localThickness: null,        // not solved yet
      firstFreezeTime: null         // not solved yet
    }
  },
  iceStructure: {
    representation: "legacy-uniform-slab",
    physicalHalfThicknessM: Number,
    physicalTotalThicknessM: Number,
    localHeightSolved: false,
    verticalVaporTransportSolved: false,
    basalMassCoupled: false,
    formationTimeSolved: false
  }
}
```

Every captured frame is **independently copied** from the existing typed solver arrays. The height model is a compact uniform descriptor, not another 401² height field allocated on every frame. A compatibility accessor reports `top=+h`, `bottom=-h` in occupied cells and zero outside physical ice support, but **does not claim** an independently evolved height.

`upgradeLegacyIceSnapshot(frame)` upgrades an existing unversioned *in-memory* Observatory snapshot, not a generic externally persisted file format. Unrecognized versions fail explicitly.

## Five transparent planes

The viewer adds a select under `ICE VIEW`:

- **Solid · legacy slab** — original optical extrusion, unchanged.
- **5 translucent height slices · preview** — at physical elevations `z = [-h,-h/2,0,+h/2,+h]`, where `h=basalHalfThickness` in SI meters. Each plane uses the **same solved 2-D basal contour**, each with 24.5% optical opacity. A fine contour outline is drawn on each plane. At oblique camera angles the five elevations produce true physical parallax; TOP gives an intentionally overlapped semitransparent silhouette. No invented ridge, concentric band, virtual branch or artificial taper is introduced.

Transparency is render-only: `depthWrite=false`; all GPU geometries and materials are disposed on switching. The old mesh is likewise disposed exactly once. `data-ice-slices` on the WebGL viewport reports the number of actual live sheet meshes for browser verification. Mode can be switched while paused, live, or scrubbing replay without resetting the growth model.

**Known limitation:** These are five flat *optical planes*, not five independently simulated z slices, an optical refraction ray-tracer, a faceted upper surface, or a local `h(x,y,t)` field. The old scalar thickness is still out-of-plane growth not consistently included in in-plane mass/latent-heat budgets. Pass 8 will own that correction, while Pass 9 will contour a physical local-height model. Do not infer central ribs or shell terraces from this preview.

## Acceptance tests

```sh
node app/tests/pass7-snapshot-layers.test.mjs
node app/tests/run-all.mjs
```

`pass7-snapshot-layers.test.mjs` compares **12 identical physical steps** of both hex and Cartesian solvers, with snapshot/preview computations inserted in only one run; vapor/temperature/phase/heat/mass/time must remain identical. It tests version metadata, cloned snapshot arrays, five ordered z elevations with identical solved contours, legacy adapter and rejection of unsupported local-height claims.

Playwright `observatory-browser.test.mjs` additionally checks: live snapshot v1; five actual Three.js mesh objects; three repeated render-mode switches; unchanged elapsed time/deposition/radius; replay to t=0 and return to LIVE; legacy solid extrusion still working; no WebGL/console errors; and the full-resolution screenshot in the five-layer mode.

## Next physically meaningful step (Pass 8)

Introduce local top/bottom height fields `hTop[i],hBottom[i]` as solver-owned quantities and define ice mass consistently. A 2.5-D plate with thickness gradients needs an explicit reduced-order vertical vapor supply/kinetics closure; it must not silently create mass at the basal surfaces. Only once those local height fields pass a documented conservation test should the five slices become different *physical* contours or show raised central ribs.
