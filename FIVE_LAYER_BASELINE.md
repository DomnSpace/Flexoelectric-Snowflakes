# Step 1 / 15 — frozen pre-coupling numerical reference

**Status:** Step 1 baseline captured and verified: full physics suite and Chromium browser smoke passed on PR run [38733133826](https://github.com/DomnSpace/Flexoelectric-Snowflakes/actions/runs/38733133826). No production solver or viewer source changed.  
**Frozen commit (deployed working Pass 7):** [77f118c091a2fd1ae278c899e7d5f6a1715ebc54](https://github.com/DomnSpace/Flexoelectric-Snowflakes/commit/77f118c091a2fd1ae278c899e7d5f6a1715ebc54)  
**Frozen reference branch:** [baseline/five-layer-pre-coupling-20261015](https://github.com/DomnSpace/Flexoelectric-Snowflakes/tree/baseline/five-layer-pre-coupling-20261015)  
**Development branch:** `feat/five-layer-coupling-foundation-v1`  
**Provisional physics model decision:** [FIVE_LAYER_RESEARCH_DECISION.md](FIVE_LAYER_RESEARCH_DECISION.md)

## Why this checkpoint matters

Pass 7 **already exists** in the frozen version: `app/solver/ice-snapshot-v1.js` stores honest `legacy-uniform-slab` version-1 snapshots; `app/render/uniform-five-layers.js` draws 5 optically separated sheets with the **same** physical 2-D ice contour. Neither the old scalar basal thickness nor the optical sheets contain a new solved z-field. The historic `pass7-snapshot-layers.test.mjs` already checks that toggling this viewer does not alter either hexagonal or Cartesian growth.

The requested five-way competition for the same vapor inventory does **not** exist in that solver. The five identical slices in the screenshots are thus an expected correct preview of the current model, not evidence that the physical layers have independent supply. Creating five independent copies of the legacy 2.5-D engine would **double-count** the surrounding vapor supply; do not do it.

Step 1 pins exact reference measurements before any competing-growth law is added.

## Files that MUST remain unchanged through this step

These were read from the source at frozen commit 77f118c091...:

- `app/solver/hex25d.js` — physical triangular lattice growth with scalar basal thickness. Git blob `c4e53850aec6488a5e7df4b81b15a71afb408303`.
- `app/solver/growth25d.js` — independent Cartesian 2.5-D baseline.
- `app/solver/ice-snapshot-v1.js` — replay schema v1, localHeightSolved=false. Git blob `489882d460c553b8e8e6c3a28a690250113f60a3`.
- `app/render/uniform-five-layers.js` — five translucent **render-only** planes. Git blob `f58ebd94a67aeb12fb83c3c0897359a771c82bfb`.
- `app/observatory.js` — viewer orchestration; the existing 5-layer toggle stays unchanged.

The work branch modifies only **new docs/new tests/a CI run step** and the registration of a new regression case in the test runner. It does not touch the physics solver or viewer. The root live GitHub Pages branch is not changed until code review and validation.

## Deterministic reference test

    node app/tests/legacy-baseline-capture.mjs
    node app/tests/five-layer-baseline.test.mjs
    node app/tests/run-all.mjs
    node app/tests/pass7-snapshot-layers.test.mjs

Capture procedure: Node.js 22 on Ubuntu GitHub Actions, three N=41 domain cases (triangular 0°, triangular 30° with kinetic amplitude 0.5, Cartesian), spacing 0.75 µm, temperature 258.15 K, ice supersaturation 7.5%, initial circular seed 1.5 µm, 12 real steps of 0.5 s; 14 vapor and 18 heat relaxation iterations per step. Snapshots at steps 1, 6 and 12. Measures in SI: elapsed time, **legacy** lateral deposited mass, legacy latent heat, shared basal half-thickness, relative sixth angular harmonic, and interface-cell count. For every checkpoint computes SHA256 over `ice`, `phi`, `sigma`, `dT`, `heatSource` typed-array binary bytes in fixed order.

The frozen **external reference** comes from successful CI run [38732698511](https://github.com/DomnSpace/Flexoelectric-Snowflakes/actions/runs/38732698511) executing the working legacy solver, *before* adding the fixture-based regression. The expected SI metrics and all 9 SHA256 field fingerprints are committed in [app/tests/fixtures/five-layer-legacy-v1.json](app/tests/fixtures/five-layer-legacy-v1.json), and checked by [five-layer-baseline.test.mjs](app/tests/five-layer-baseline.test.mjs). This is not a self-generated expected value during the test run.

| Case | t [s] | Legacy deposited mass [kg] | Legacy total latent heat [J] | Legacy shared half-thickness [m] |
|---|---:|---:|---:|---:|
| Triangular 0° | 6 | 4.367291949234789e-12 | 1.2536529904875921e-5 | 1.326911964635121e-5 |
| Triangular 30° | 6 | 3.837836659848035e-12 | 1.1016702023926778e-5 | 1.3138286221884806e-5 |
| Cartesian | 6 | 1.1907328923225283e-11 | 3.4180583040564336e-5 | 2.5627153036626944e-5 |

**Important:** N=41 is a *short, deliberately tiny numerical fingerprint*, not a claim of converged/physically valid snow-crystal morphology. The existing legacy basal height is not included in the deposited mass and latent-energy tally. The frozen test intentionally records even these deficiencies so the next model cannot silently inherit them while calling its accounting “conservative.”

The SHA256 comparison is intentionally strict for the documented Node-22/little-endian CI environment; nonidentical platform math may require comparisons of numerical summaries at tolerances rather than raw binary hashes. The regression test also separately checks physical scalars using declared tolerances.

## Legacy behavior checks

- The existing Pass 7 parity test advances 12 steps on both hex and Cartesian solvers with snapshot creation and five-layer plan generation in only one branch; it checks identical arrays, mass and time against the no-render control.
- The existing Playwright smoke test checks the current renderer has exactly five live transparent meshes, the new UI label still says **PREVIEW ONLY**, switching modes does not affect mass/time, and the entire viewer works on Chromium with WebGL.
- CI must pass the full physics suite and browser smoke on this new work branch before merge.

## Next step in the 15-step checklist

Step 2: formalize candidate B's dimensions, face ownership, initial gas/ice inventories, closed-boundary conservation and explicit boundary fluxes. Before committing any viewer-driven ice growth, test a **small finite-volume vapor transport toy model with zero deposition and zero external flux**: gas mass must remain invariant; then add one deposition sink and check ice gain equals gas loss. Compare 1/2/4 gas halo slabs and Δz. Only after these tests should the five optical sheets become independently simulated ice layers.

The snapshot schema must become `v2` for this new simulation, while `v1` remains loadable for unchanged older replay.
