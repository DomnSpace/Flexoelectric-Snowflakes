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
- Three synchronized plots: six radial vapor traces, radius and basal thickness over time, and m=6 / m=4 angular modes.
- Replay slider inspects stored numerical snapshots; LIVE selects the newest snapshot. Replay does not overwrite solver history.
- Data CSV exports captured model history with units. View PNG saves the WebGL crystal viewport.

## Physical scope and limitations

The solver in app/solver/growth25d.js includes nine-point isotropic 2-D quasi-steady vapor diffusion, latent thermal conduction, 2.5-D interface kinetics, Gibbs–Thomson correction, prism/basal attachment, optional SDAK edge enhancement, and angular modes from growthMorphology(). Three.js only extrudes the numerical contour.

This is **not** a complete 3-D growth model, a calibrated atmospheric forecast, or a prediction of exact natural dendrites. There is no coupled airflow velocity field, full 3-D vapor/heat dynamics, or electrically driven flexoelectric term in this viewer. Parameters from the legacy electric-field toy should not be misrepresented as coupled to this resolved-growth solver.

The uploaded SnowflakeWorkSample.cpp ends inside the relax function, and cannot be compiled verbatim; the existing JS physics module is used as the live source of truth.

The fixed 0.5-second model step and bounded 360-frame history make replay deterministic and bounded in memory. A 73 × 73 grid at 0.75 micrometers per cell is used. The UI stops as the growing crystal nears the finite computational boundary.

## Validation

    node --check app/observatory.js
    node app/tests/observatory-contract.test.mjs
    node app/tests/run-all.mjs

The source-level contract test verifies expected DOM bindings, real solver imports, and numerical snapshot fields. The existing numerical suite covers physical modules. Browser-level WebGL rendering and CDN retrieval still require verification on a modern Chromium browser. GitHub Actions runs syntax and Node tests on the feature branch and pull requests.
