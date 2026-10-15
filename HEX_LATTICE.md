# Experimental triangular basal-plane lattice (Pass 4)

## What the lattice represents

The classic snowflake has basal-plane **sixfold crystallographic symmetry**, arising from hexagonal ice Ih. Ordinary ice Ih is a hydrogen-disordered, tetrahedrally hydrogen-bonded 3-D oxygen network. Its nearest oxygen–oxygen spacing is about **0.275 nm**, the tetrahedral neighbor angle is about **109.5°**, and individual molecular bond directions are not planar sixfold bonds. See Petrenko & Whitworth, *Physics of Ice*, Ch. 2 (2002), https://doi.org/10.1093/acprof:oso/9780198518945.003.0002 and the 2005–2017 Libbrecht reviews. The isolated water H–O–H angle is ~104.5°; this is **not** the hexagon angle.

**This solver does NOT simulate atoms or hydrogen bonds.** Each point is a **mesoscopic 0.75 µm transport volume**, roughly 2,700 oxygen–oxygen distances in the plane. The sixfold grid represents basal-plane directions for vapor/thermal transport. Crystallographic attachment anisotropy is still an independently selected model parameter. Do not infer the right attachment coefficients from the lattice alone.

## Triangular primal lattice → hexagonal dual cells

On axial coordinates `(q,r)`, with edge-to-edge nearest-site distance `a=0.75 µm`, the physical site centers are

    x = a(q + r/2)
    y = a sqrt(3)r/2

Six equidistant neighboring centers have axial offsets

    (+1,0), (-1,0), (0,+1), (0,-1), (+1,-1), (-1,+1).

Their Voronoi **control volumes are regular hexagons** of area `sqrt(3)/2 a²`; their separating face length is `a/sqrt(3)`. Each site has one vapor, temperature and phase-fraction value. Those values are updated by the solver before visualization. The field map draws the control cells themselves; no cosmetic hexagon filter is applied to the former square-grid data.

The six-neighbor isotropic diffusion stencil is

    Laplacian(u) = [2/(3a²)] (sum_6(u_neighbor) - 6u) + O(a²).

A quadratic consistency check tests `Laplacian(x²+y²)=4` and `Laplacian(linear)=0` across multiple sites.

**Face transport/cell balance:** interface growth uses a combined vapor-diffusion/attachment resistance, local vapor/temperature Gibbs–Thomson equilibrium, and occupied-face area multiplied by the shared modeled ice thickness. Deposited mass increments the local phase fraction and sources the next heat step, with the same *in-plane* mass budget limitations as the older Cartesian solver.

## Why a hex grid is not sufficient to prove hexagonal snowflakes

A triangular mesh is itself sixfold. Even with ice attachment anisotropy set to zero, it may seed an artificial `m6` morphology. The absence of the old square grid does not establish **physical** crystal symmetry. The viewer exposes both `m6` and `m4` but **neither is proof of resolved six arms**. Quantitative acceptance requires **orientation-sweep tests, refinement at constant physical size**, and control runs with anisotropy zero. The natural snowflake's sixfold appearance must come from experimentally grounded attachment/surface kinetics, not simply from the 60° transport directions.

The 2017 quantitative 3-D phase-field snow-crystal study adds anisotropic kinetic, surface diffusion and anisotropic surface-energy terms and succeeds in reproducing much more morphology: https://www.nature.com/articles/s41524-017-0015-1 . Our 2.5-D solver **does not** implement this complete constitutive model.

## Finite hexagonal boundary and automatic domain expansion

The transport domain uses hexagonal shell radius `R=(N-1)/2` and includes sites satisfying `max(|q|,|r|,|q+r|) <= R`. Outer-shell sites are fixed to far-field vapor excess and zero temperature anomaly (Dirichlet). This is a **hexagonal** boundary now: it can introduce m6 bias rather than the former square m4 bias. A finite boundary is not an infinite cloud. The UI always labels the boundary.

The domain audits all occupied phase cells. At **12 remaining neighbor rings** the optional automatic domain extension copies the real ice/phase, vapor, heat and source arrays into the next domain size 73→97→129→161, without resetting the ice, elapsed model time or accumulated deposition/latent heat. Newly added gas is initialized at far-field vapor and zero thermal perturbation, then relaxed; that is a **numerical approximation**, not a dynamically infinite PDE solution. The viewer records the number of expansions and CSV includes both topology and domain size. At **6 rings** or past the maximum grid size the stop remains enforced.

The exact domain's physical horizontal diameter is `(N-1)a` for a hex shell at radius R; it is NOT an atomic length. Increasing N expands the accessible area at constant site spacing and therefore does not improve microscopic resolution.

The automatic extension is only implemented for hex mode. The Cartesian mode retains the earlier square-domain stop and existing separate electrical module.

## Electric physics limitation

`app/solver/electric25d.js` currently uses Cartesian face-dielectric finite differences. It has NOT been ported to the triangular Voronoi cells. While the hex mode is active, electric/field controls are explicitly disabled and the viewer reports "not coupled." Pretending the Cartesian potential array corresponds to hex grid positions would be invalid.

## Verification / numerical acceptance

    node app/tests/hex25d.test.mjs
    node app/tests/run-all.mjs

Regressions check six symmetric equidistant neighbors, expected site count `1+3R(R+1)`, axial↔Cartesian coordinate roundtripping, linear/quadratic Laplacian consistency, thermally sourced positive deposition, finite fields and morphology, tip-to-hex-boundary detection, and conservation of stored phase/time/mass through domain expansion.

Browser CI checks the six-neighbor solver is *actually run*, the hex raster is used, the Cartesian solver continues to work, and electric coupling is not silently reused on an incompatible grid.

## Unsolved

- A genuine molecular-size ice Ih lattice, hydrogen disorder and 3-D tetrahedral geometry.
- A fully conservative basal growth budget and resolved vertical heat/vapor transport.
- An experimentally calibrated temperature-dependent sixfold anisotropy and SDAK.
- Quantitative convergence, surface-diffusion and charge coupling, branching validation, moving adaptive mesh, full 3-D interface.
