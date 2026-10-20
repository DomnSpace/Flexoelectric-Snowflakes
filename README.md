# Snowflake Observatory

Runnable interactive front end for the repository's **resolved 2.5-D ice growth** implementation. Open the repository root through a local HTTP server. The root legacy 3-D demo has been replaced with a physics-first instrument; the activation, supercooling and field labs in app/ remain intact.

## Windows one-click

Double-click **run-viewer.cmd** from the repository root. Python 3 (via py or python on PATH) and a modern browser are required. The command opens a local HTTP server at http://localhost:8088/.

Alternately, from the repository root:

    python -m http.server 8088

Then visit http://localhost:8088/. Opening index.html via file:// is not supported because the physics is imported as ES modules.

**Online requirement:** The 3-D renderer uses pinned Three.js 0.160.1 modules from jsDelivr. The physics and chart data use local project modules. The page needs internet for the 3-D library, and is not a fully offline bundle.

## Simulator and interaction

- Play, Pause, Step, Reset, and Dendritic preset.
- Change environmental and facet parameters, then Apply / Restart solver. The pending badge distinguishes unsimulated slider changes from the active model configuration.
- Orbit using drag, zoom using wheel, double-click to reset; TOP / ISO view presets. The solid is a 3-D prism extrusion of the actual 2-D phase-field contour, with depth from model basal thickness.
- Spatial heatmap: reference-normalized vapor supersaturation (u = c_v / c_sat(T_inf) - 1), excess temperature, and actual ice fraction. Hover to inspect grid cells.
- Three synchronized plots: six radial vapor traces, radius and basal thickness over time, and m=6 / m=4 angular modes. Additional solver-derived field maps expose the electrostatic potential and |E| with live residual, polarization and equilibrium-correction diagnostics.
- Replay slider inspects stored numerical snapshots; LIVE selects the newest snapshot. Replay does not overwrite solver history.
- Data CSV exports captured model history with units. View PNG saves the WebGL crystal viewport.

## Physical scope and limitations

The solver in app/solver/growth25d.js includes nine-point isotropic 2-D quasi-steady vapor diffusion, latent thermal conduction, 2.5-D interface kinetics, Gibbs–Thomson correction, prism/basal attachment, optional SDAK edge enhancement, and angular modes from growthMorphology(). Three.js only extrudes the numerical contour.

This is **not** a complete 3-D growth model, a calibrated atmospheric forecast, or a prediction of exact natural dendrites. There is no coupled airflow velocity field or full 3-D vapor/heat dynamics. The optional electrostatic module solves ∇·D=0 with spatial dielectric contrast and imposed flexoelectric polarization from a **prescribed**, not mechanically solved, strain gradient. Its dielectric energy correction to interface equilibrium is separately opt-in and not a validated growth law. See [ELECTROSTATICS.md](ELECTROSTATICS.md); the legacy electric-field toy is not reused.

The uploaded SnowflakeWorkSample.cpp ends inside the relax function, and cannot be compiled verbatim; the existing JS physics module is used as the live source of truth.

The fixed 0.5-second model step and bounded 360-frame history make replay deterministic and bounded in memory. A 73 × 73 grid at 0.75 micrometers per cell is used. The UI stops as the growing crystal nears the finite computational boundary.

## Validation

    node --check app/observatory.js
    node --check app/solver/electric25d.js
    node app/tests/electric25d.test.mjs
    node app/tests/observatory-contract.test.mjs
    node app/tests/run-all.mjs

The source-level contract test verifies expected DOM bindings, real solver imports, and numerical snapshot fields. The existing numerical suite covers physical modules. Browser-level WebGL rendering and CDN retrieval still require verification on a modern Chromium browser. GitHub Actions runs syntax and Node tests on the feature branch and pull requests.


## Default stellar ice habit / thermal five-layer coupling (Step 4/15)

The Observatory now **starts with the triangular stellar experiment at a 30° kinetic phase**, 14% ice supersaturation, 0.50 sixfold anisotropy and SDAK on, using its physical 129-index hexagonal domain with automatic growth up to 300 µm. These settings reproduce the explicitly selected stellar research configuration without presuming that the 30° orientation is universally correct for ice Ih. Choosing another habit or using Cartesian/electric controls still works, and changing back to Stellar reloads 30° instead of reverting to 0°. See [FIVE_LAYER_THERMAL.md](FIVE_LAYER_THERMAL.md) for the physics and tests.

The independent [five-layer Model B laboratory](app/layered-lab.html) now runs with **coupled latent-heat conduction by default**, adjustable to an uncoupled no-heat control. The five ice maps are actual separately deposited mass/volume fractions sharing vapor; the optional temperature maps, maximum ΔT and energy residual are measured from an implicit conservative voxel gas/ice heat solver. This remains a **small-grid experimental 3D-transport reduction**, not a converged 300 µm faceted 3D ice solver or molecular snowflake rendering.

## Model B — five genuinely competing ice slabs (Step 3/15)

[Launch the independent five-layer physical growth laboratory](app/layered-lab.html) to inspect **five separate calculated ice-occupation fields**, their layer-specific vapor deposition and a shared finite-volume gas inventory. [FIVE_LAYER_MODEL_B.md](FIVE_LAYER_MODEL_B.md) provides the equations, computational compromises and exact conservation tests.

A backward-Euler implicit solver replaces Step 2's nanosecond-limited explicit transport in this opt-in experiment. Partial ice volume fractions now grow from measured water mass and displace local vapor conservatively. The top/bottom can grow faster than the shielded central slab — the model does not hardcode a center spine. **This remains a coarse 3-D vapor experiment without a thermal PDE or calibrated full-surface kinetics.** The main Observatory and its legacy five-slice preview are unchanged.

## Five-layer competition — step 2/15 conservative gas finite volumes

[FIVE_LAYER_STEP2_FV.md](FIVE_LAYER_STEP2_FV.md) describes the opt-in numerical pilot on triangular 2D × vertical gas slabs with **one shared vapor mass inventory**, genuine lateral/vertical diffusion and attachment sinks on exposed ice faces. Five ice-eligible layers plus gas-only upper/lower halos are independently inventoried, but their fixed geometry does not yet move. Positive flux limiting, sealed/open boundary mass conservation and interface latent-heat **bookkeeping** are tested. The old 2.5-D snowflake solver and existing five-layer optical preview remain unchanged. The explicit time-step restriction is extremely short at micron spacing; this is a validation kernel, **not yet an interactive 300 µm model**.

## Five-layer competition roadmap — baseline step 1/15

[FIVE_LAYER_RESEARCH_DECISION.md](FIVE_LAYER_RESEARCH_DECISION.md) records the evidence and provisional choice of **five ice-eligible layers sharing a vertically resolved finite-volume vapor/heat environment**, with separate gas halo slabs. Crucially, it does **not** prescribe fastest central-layer growth: all exposure and attachment preferences must emerge from computed flux and kinetics. It treats this as a coarse 3D transport discretization, not solved molecular ice or an already calibrated five-layer model.

[FIVE_LAYER_BASELINE.md](FIVE_LAYER_BASELINE.md) pins the exact immutable pre-coupling Pass 7 commit and a frozen branch. The legacy hex/Cartesian solvers stay unchanged. A new independent capture produces nine deterministic SHA256 references for vapor, heat, ice, phase and sources (three physical checkpoint steps × three solver cases); SI mass/latent-heat/basal-height results are committed in `app/tests/fixtures/five-layer-legacy-v1.json`. Run `node app/tests/five-layer-baseline.test.mjs` to verify that new modules have not unintentionally changed the baseline.

**No new vertical ice growth, local solved thickness or vapor competition is claimed by Step 1.** The five transparent layers remain explicitly labeled a uniform-slab preview. Model B's conservative gas/ice face-flux interface is the next implementation gate; the original solver and current GitHub Pages remain intact.

## Pass 10 — Experiment presets, humidity validity and growth histories

[PASS10_GROWTH_REGIMES.md](PASS10_GROWTH_REGIMES.md) documents the new observational candidate presets, temperature-dependent water-saturation reference, the staged **−15 °C: 3 → 13 → 5 → 14% ice supersaturation** continuous growth-history experiment, independent basal σ₀ loading, and versioned environment metadata/CSV. **Stage preset** prepares complete conditions without restarting the active crystal; **Run experiment** applies those conditions. Time-varying environments rebase vapor from preserved absolute gas density before resetting the OPEN Dirichlet far-field boundary. Needle/capped column entries remain unsupported by the 2.5-D solver. The attachment coefficients are explicitly uncalibrated.

## Pass 9 — Real local surface reconstruction and ice ridge research

[PASS9_HEIGHT_AND_RIDGE_RESEARCH.md](PASS9_HEIGHT_AND_RIDGE_RESEARCH.md) reviews observed snowflake **macrosteps, step bunching, and ridge formation**, explicitly rejecting ad hoc van der Waals alignment/rotation as a manufactured ridge source. Pass 9 reconstructs the live opt-in Pass 8 local top/bottom height surface and ice rim on the physical triangular grid and clips **five actual z cross-sections**. Choose **ICE VIEW → Physical local top / bottom / rim** or **5 actual local iso-height sections** after enabling **Resolve local upper/lower ice heights**. The old uniform views remain. Surface details only appear if the solver produced them: no ridge deposition or artificial branches are added. Render LOD on 300 µm domains is explicitly distinguished from the 0.75 µm physical lattice and does not change mass, time, vapor or heat calculations.

## Pass 8 — Local ice thickness and explicit prism/basal ice-mass ledger

[PASS8_LOCAL_HEIGHT.md](PASS8_LOCAL_HEIGHT.md) documents the **opt-in** locally varying ice `hTop`/`hBottom` model on the triangular lattice. It computes actual local mass from the phase fraction and thickness, adds both prism and basal deposits to latent heat, tracks distinct external planar and unresolved z-direction vapor supplies, and records an ice-side mass ledger residual. **This is not a fully conservative 3-D vapor-ice solution**: the vertical reservoir has no resolved concentration/thermal PDE. In the Observatory, enable **Resolve local upper/lower ice heights**, apply the solver, then choose **Spatial Fields → Local ice thickness h**. Five optical sheets remain a *mean-thickness visual proxy* until Pass 9.

## Pass 7 — Versioned ice snapshots and five transparent layers

[PASS7_STATE_AND_LAYERS.md](PASS7_STATE_AND_LAYERS.md) describes the versioned SI snapshot boundary, a render-only **five-sheet preview** at five physically spaced z elevations, and tests that the solver's deposition, heating and fields remain unchanged. In the Observatory, switch **ICE VIEW → 5 translucent height slices · preview** to compare with the legacy solid slab without resetting the simulation. The five sheets **currently share the same measured basal contour** because the original solver has only one uniform height. This is an honest optical preview for Pass 8's locally varying mass-conserving height model, **not** a calculated set of internal crystal ribs.

## Crystal orientation / tip-splitting investigation

[ANISOTROPY_PHASE_STUDY.md](ANISOTROPY_PHASE_STUDY.md) explains the physically separated and falsifiable hex solver experiments: a selectable 0°/15°/30° phase in the **local surface-normal** attachment factor, optional positive anisotropic surface stiffness (gamma + gamma-double-prime), and a real rotated hexagonal ice seed. These are NOT dielectric effects; the electric module remains Cartesian-only. The controls default to the original physics and do not enforce six stems or manufacture morphology. Compare with zero kinetic anisotropy and the phase-sweep test results before claiming the 30° setting suppresses splitting.

## 300 µm live hex-domain / crisp field display

[DOMAIN_300UM.md](DOMAIN_300UM.md) details the physically **300 µm** hexagonal domain (401 axial index width, 120,601 active sites at unchanged 0.75 µm spacing), optimized six-neighbor relaxation, the automatic domain ladder, bounded replay memory and 768/1024/1536-pixel *display-only* field interpolation.

The basal-plane ice contour uses 720 phase-derived angular samples for the experimental hex mode, as opposed to the earlier 180, but the 3D rendering remains a basal contour extrusion, not a resolved 3D ice crystal. Larger arrays do not validate branching or replicate molecular hydrogen-bond geometry.

## Experimental 60° triangular transport grid

[HEX_LATTICE.md](HEX_LATTICE.md) describes the optional six-neighbor triangular transport solver, true hexagonal Voronoi field visualization, automatic 73→97→129→161 hex-domain expansion without restarting the ice, and why neither this lattice nor oxygen/hydrogen molecular bond geometry automatically predicts six-armed snowflakes.

Use the **Numerical transport lattice** selector to compare the existing Cartesian nine-point scheme with the new triangular six-neighbor scheme. Press **Apply / restart**, or load one of the 2.5-D plate/dendrite candidate habits (which now select triangular transport by default). The electric Poisson experiment remains Cartesian-only; it is disabled and explicitly identified as unavailable in triangular mode. For honest comparisons, inspect both m4 and m6 and vary kinetic anisotropy, domain size and site spacing.

The hex solver is an **experimental alternate 2.5-D model**, not calibrated molecular simulation or proof of sixfold morphological instability. Its adaptive boundary still imposes a finite, hexagonal Dirichlet far field.

## Grounded habit scenarios / numerical validity

See [SCENARIO_VALIDITY.md](SCENARIO_VALIDITY.md) for the literature basis, exact temperature/supersaturation presets, adjustable solver-domain geometry, basal-thickness warnings, finite-domain stop conditions and **limitations that prevent any guarantee of dendrites**.

The Habits control loads three exploratory, neutral-field, **−15 °C plate-regime** scenarios using the same solver. The `−5 °C needle/column` reference is explicitly disabled: it requires 3-D vapor/heat and face-specific growth rather than a 2-D outline extruded to a fake needle. Domain size is 73²–161² at constant 0.75 µm spacing; grid-point width is `(N−1)dx`, not `N dx`. The UI tracks minimum remaining halo to the square Dirichlet boundary and actual plate thickness relative to tip radius, and pauses runs that leave stated validity limits. Mean basal attachment prefactor A is now visible instead of being hidden at A = 1, but it remains **uncalibrated**.

The previous very tall-column screenshot is an expected sign that out-of-plane growth was unconstrained by a resolved z transport model. The present guard prevents presenting such objects as validated snowflakes; it does not add that missing physics.

## GitHub Pages — browser-only launch

GitHub Pages is already enabled for this public repository. Its current setting was verified through the Pages API on 2026-10-14:

- **Publishing mode:** Deploy from a branch (legacy).
- **Current published branch:** `feat/flexoelectric-snowflake` at `/(root)` (the earlier viewer).
- **Confirmed Pages address:** https://domnspace.github.io/Flexoelectric-Snowflakes/

To publish the complete current Observatory with electrostatics, change one repository setting:

1. Open [Repository Settings → Pages](https://github.com/DomnSpace/Flexoelectric-Snowflakes/settings/pages).
2. Keep **Build and deployment → Source: Deploy from a branch**.
3. Set **Branch: `feat/snowflake-electrostatics-v1`**, **Folder: `/(root)`**, then click **Save**.

GitHub Pages will publish directly from that feature branch, and future commits pushed there will trigger Pages rebuilds. The app needs no build tooling: root `index.html` imports its local solver modules and Three.js from jsDelivr. All navigation and imports use relative paths compatible with `/Flexoelectric-Snowflakes/`. There is no backend or data storage; computations run in the visitor's browser.

**Why this setting is necessary:** the repository's existing Pages site is configured for branch publishing, and the Github connector does not expose a Pages-settings edit. An attempted custom Actions deployment packaged and passed all 13 physics tests but could not replace the current branch-based publishing configuration. We removed the redundant Actions publishing workflow rather than leaving broken CI or silently overriding the existing site.

For local development, `run-viewer.cmd` or `python -m http.server 8088` remains available. GitHub Pages serves the same source.

