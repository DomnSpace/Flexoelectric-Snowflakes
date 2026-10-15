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

