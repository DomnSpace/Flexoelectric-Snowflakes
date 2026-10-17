# PASS 9 — Actual surface meshes from local ice height (NO invented ridge force)

## Research: what makes a snow crystal's main axial ribs?

**Do not add an ad hoc “ridge attraction” force or reinterpret ice's H-bond angles as a micron-scale van der Waals torque.** Ice Ih crystallographic basal/prism orientations come from its hydrogen-bonded hexagonal lattice; crystal *rotation* changes which directions align with the external field/lattice but does not by itself create raised ribs. The electric Poisson module remains Cartesian-only, so the live hex mode has no resolved electric ridge effect.

Two targeted results explain why ridge *accretion* is a different, unresolved question:

1. **Gravner & Griffeath**, *Modeling snow-crystal growth: A three-dimensional mesoscopic approach*, Phys. Rev. E **79**, 011601 (2009), DOI https://doi.org/10.1103/PhysRevE.79.011601, model vapor diffusion, anisotropic attachment and a boundary layer. They reproduce snow crystal **ridges, ribs, sandwich plates and hollow columns**. In their case study of plates and ridges, elevations along the main branches arise from **early mesoscopic macrosteps** and continue growing slowly in z, accompanied by flanking grooves. Original manuscript https://www.math.ucdavis.edu/~gravner/papers/h3l.pdf , “Case study I: Ridges and plates”.
2. **Three-dimensional aspects of sidebranch formation during the growth of snow crystals**, *Journal of Crystal Growth* **548**, 125846 (2020), DOI https://doi.org/10.1016/j.jcrysgro.2020.125846, reports interferometric imaging of **step bunching on prism and basal crystal faces** during sidebranch growth. It is observed microstructural evidence, not a calibrated term we can insert without modeling steps.

More broadly, **Libbrecht**, *The physics of snow crystals*, Rep. Prog. Phys. **68**, 855–895 (2005), DOI https://doi.org/10.1088/0034-4885/68/4/R03, emphasizes vapor/heat transport and surface-attachment kinetics. SDAK and premelting-related surface diffusion are serious candidate mechanisms, not an unconditional requirement that every arm has one ridge.

### Decision

Do **not** introduce a numerical ridge source without an independent prognostic model of facet steps / terrace widths, basal/prism attachment, step bunching and vertical mass sourcing. Pass 9 is strictly the previously planned **render-side implementation of physical height geometry**. A natural rib appears only if Pass 8's current hTop/hBottom actually create it. No extra molecule is deposited by rendering.

## Geometry

Source: copied, versioned **snapshot-v2** physical fields:

- the triangular grid: \(x=a(q+r/2), y=\sqrt{3}ar/2\);
- ice phase fraction \(\phi(q,r)\), with physical solid mask;
- local top \(z_t(q,r)\) and bottom \(z_b(q,r)\) from Pass 8.

\`app/render/actual-height-mesh.js\` iterates both *equilateral triangles* per axial index parallelogram, ignoring triangles crossing the hex domain mask. Sutherland–Hodgman clipping of each triangle against \(\phi\ge 0.5\) yields the footprint. It triangulates the upper/lower reconstructed surfaces and generates vertical sidewalls **only** on phase-isoline boundary segments. Multiple disconnected areas and interior holes are retained by construction; there is no radial six-arm outline or Wulff template. Surface normals come from the resulting actual mesh and control the existing Three.js lights.

For each of **five** z planes at relative \([-0.8,-0.4,0,0.4,0.8]\) of the maximum model half-thickness, the renderer intersects actual sampled top/bottom fields:

\[
\phi(x,y)\ge0.5,\quad z_b(x,y)\le z \le z_t(x,y).
\]

Each of the five transparent planes therefore shows its **own clipped physical cross-section**, including gaps/holes where supported. Uniform real local heights deliberately give **identical** five contours. No fake tapered insets or manufactured middle ridges.

### Important quantitative distinction

The solver's ice ledger is finite-volume:

\[
M_{\rm solver}=\rho_{\rm ice}A_{\rm hex}\sum_i\phi_i(z_{t,i}-z_{b,i}).
\]

A display mesh using the **thresholded interpolated** \(\phi=0.5\) footprint approximates the solid geometry at finite resolution; its triangle volume is **not guaranteed equal to this finite-volume mass**. It does not determine the real physical deposition or modify any field. The viewer labels the LOD and phase threshold. Do not use the rendered mesh volume to validate mass conservation. Pass 8's independent ice-side ledger remains the conservation metric. The unresolved vertical vapor reservoir remains *open*, not a solved 3-D PDE.

### Performance / mobile

N≤161 renders all physical triangles at the original a=0.75 µm. Larger N=201–401 renders at **stride 2** (approximately a=1.5 µm mesh sampling) to limit mobile GPU load. This is ONLY render LOD: the solver still tracks all 120,601 hexagonal cells at 401 and still uses a=0.75 µm. A thin sidebranch below the render LOD can disappear visually; the full-resolution 2-D field still contains the numerical data. The UI reports stride and generated triangle counts. The meshes are disposed on mode switches; legacy v1 and Cartesian electric rendering are left intact.

## Feature selection

In ICE VIEW, for a live Pass 8 height run:

- **Physical local top / bottom / rim**: actual v2 upper/lower surfaces, sidewalls, gradient lighting.
- **5 actual local iso-height sections**: actual intersected z-levels; up to five different physical footprints.
- The old **solid legacy slab** and **five uniform optical guide** modes remain available for A/B comparison.

Without Pass 8 local-height enabled, the new physical modes are disabled. Changing the view or scrubbing snapshots NEVER advances growth or changes ice mass. The core growth equations and open prism/basal fluxes are deliberately untouched.

## Validation

\`app/tests/pass9-actual-height.test.mjs\`: synthetic annular crystal with real hole and a known thick ring must give distinct center vs upper/lower slice areas, wall segments and finite surface triangles. A constant-height crystal must yield five equal footprints; changing the renderer must leave time, mass, all phase, vapor, temperature and height fields bitwise unchanged. Real Pass 8 frames, including a large-domain render LOD, are included.

Playwright: actual physical-solid, five real slice objects, v2-only state gating, mesh metadata, no solver changes across view transitions, preserved v1 old mode and electric Cartesian baseline.

## What still needs a separate *physics* pass

- A prognostic density / propagation law for mesoscopic macrosteps on basal surfaces, prism/basal step bunching, and terrace-width SDAK. This needs quantitative experimental parameters; **not an attractive procedural trunk-ridge mask**.
- Spatially conserved attachment/remobilization and measured ridge-growth z-flux; don't double-count basal ice supply. The existing five-slab finite-volume research prototype may serve as a better foundation for step-driven z flux than applying new arbitrary biases to the live 2.5-D code.
- Refinement and transport convergence, direct sidebranch 3-D micrographs and checks of the open z-reservoir's validity.
