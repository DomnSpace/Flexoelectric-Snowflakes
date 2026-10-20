# 300 µm hexagonal domain / high-resolution display pass

This pass addresses the screenshots of the 120 µm, 161²-node hexagonal ice-growth model where arms already reach a 6-ring Dirichlet halo and the crystal contour and vapor raster look coarse.

## What "300" means

The target is **300 µm across the physical computational hexagon**, not 300 nodes:

    a = 0.75 µm nearest-neighbor site spacing (unchanged)
    N = 401 axial-index rows/columns (not a 401×401 fully occupied square)
    hex shell radius R = (N−1)/2 = 200
    hexagonal domain across vertices = 2 R a = 300 µm
    active sites = 1 + 3R(R+1) = 120,601

The permitted progression is **73, 97, 129, 161, 201, 257, 321, 401** lattice-index widths, corresponding to 54, 72, 96, 120, 150, 192, 240, 300 µm. The expanded lattice retains **the same physical mesh spacing** and therefore does not by itself refine microscale attachment kinetics.

Automatic hex-domain growth remains on by default; it expands at 12-ring clearance (before reaching the 6-ring hard stop), preserving the full phase/ice field, time, numerical flux/temperature/source fields and cumulative model budget. At 401 the solver must still stop near its fixed hexagonal Dirichlet boundary. This is not an infinite unbounded cloud.

The existing Cartesian branch, including the Poisson/flexoelectric experiment, remains selectable and limited to its previously validated ≤161² physical domain; large-N hex transport does **not** automatically support electric physics.

## Numerical optimization

The previous Jacobi six-neighbor stencil invoked a JavaScript callback for each neighboring cell in every sweep, with nested full-square iteration. A 401² array makes that expensive. The hex solver now caches the active interior and Dirichlet-shell indices, and evaluates the same six-site flux stencil by direct typed-array neighbor offsets. It holds both ping-pong buffers to the same Dirichlet far field and updates evolving ice nodes' vapor equilibrium each iteration. This is a computational implementation change, not a new PDE or approximation to the six-neighbor operator.

Test the 300 µm solver independently:

    node app/tests/hex25d.test.mjs
    node app/tests/run-all.mjs

Tests check actual active-site count, 300 µm width, finite deposition and temperature, and conservative continuation from 321 to 401 index-width grid. CI Chromium smoke also checks actual >0 model-time advancement on 401 with the raster drawn.

## Higher fidelity, independent of numerical physics

The old field panel painted 600px pixels by manually drawing a hexagon for each control cell. At 120,601 active sites this would be slow and magnify jagged nearest-cell boundaries. The new field panel samples a **continuous piecewise-linear interpolation of the same three neighboring triangular solver sites** at 768, 1024 or 1536 display pixels. It does not introduce extra physical degrees of freedom or pretend atomic-scale fields. Vapor and temperature are interpolated within the triangular mesh; solid phase is masked using the interpolated physical phase fraction.

The 3D rendered basal outline now uses **720 angular rays** rather than 180 for the hex solver and interpolates physical phase values at subcell positions to locate the 0.15-occupancy iso-contour. This is a more faithful depiction of the stored phase field, not an artist-created or precomputed ice-flake polygon. Reentrant contours and disconnected islands are still not fully represented by a single radial outline; a future marching-triangles multi-loop extractor would be preferable for mature dendrites.

The higher-resolution 1024/1536-pixel raster is independent of the solver's a=0.75 µm; there is still no molecular ice Ih simulation. Larger CPU/GPU costs can be reduced using the 768-pixel selection without changing the physics.

## Important unsolved conditions

- Boundary-induced m6 is not separated from physical anisotropy or ice surface kinetics.
- No full 3D thermal/vapor transport, basal mass conservation, surface diffusion, or experimentally calibrated branching.
- Larger domains require more iteration for an equivalent steady-state convergence level: the optimized Jacobi relaxation is not an adaptive multigrid or converged infinite-domain solution. Growth histories should be compared only with numerical convergence checks.
- The contour is a radial reconstruction with 720 samples, not a genuine 3D ice surface and not guaranteed to capture interior cavities or disconnected protrusions.
- Longer runs can still stop near the 300 µm boundary. This is a guard, not a failure of the branch.

See HEX_LATTICE.md and SCENARIO_VALIDITY.md for the existing 2.5-D and environmental interpretations.
