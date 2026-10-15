# Decision record FL-0001 — competitive five-layer ice growth

**Status:** proposed reduced model; research-backed mechanism, numerical parameters uncalibrated.  
**Date:** 2026-10-15  
**Baseline commit:** `77f118c091a2fd1ae278c899e7d5f6a1715ebc54`  
**Frozen baseline branch:** `baseline/five-layer-pre-coupling-20261015`  
**Development branch:** `feat/five-layer-coupling-foundation-v1`

## Evidence evaluated (versus open hypotheses)

1. J. Gravner and D. Griffeath, *Modeling snow-crystal growth: A three-dimensional mesoscopic approach*, **Phys. Rev. E 79**, 011601 (2009), DOI [10.1103/PhysRevE.79.011601](https://doi.org/10.1103/PhysRevE.79.011601). A computationally feasible 3-D mesoscale model incorporating vapor diffusion, anisotropic attachment and boundary-layer dynamics reproduces ribs, sandwich plates, faceted structures and branches. This **supports testing genuine vertical gas transport**. It does not establish that a particular five-slab reduction is accurate.
2. K. G. Libbrecht, *The physics of snow crystals*, **Rep. Prog. Phys. 68**, 855–895 (2005), DOI [10.1088/0034-4885/68/4/R03](https://doi.org/10.1088/0034-4885/68/4/R03). Vapor diffusion, heat diffusion, kinetic attachment and morphological instability all matter.
3. K. G. Libbrecht, *Characterizing Structure-Dependent Attachment Kinetics near −14 °C* (2020), [arXiv:2009.08404](https://arxiv.org/abs/2009.08404). Measurements support reduced attachment barriers on **narrow prism facets** near −14 °C (SDAK); naive equal rate on all ice interfaces is not physically adequate.
4. K. G. Libbrecht, *On the Molecular Dynamics of Structure Dependent Attachment Kinetics* (2020), [arXiv:2012.12916](https://arxiv.org/abs/2012.12916). Edge-enhanced surface diffusion is a plausible microscopic mechanism, but its quantitative details are uncertain.

No external user-provided research report has been received. This record is a **provisional model-selection decision**, not a validated set of five-layer kinetic constants. Do not claim that a middle layer must win; exposed outer basal faces may attract more vapor.

## Competing approaches

| Model | Real vertical competition? | Geometry / physics risk | 300 µm browser cost | Decision |
|---|---|---|---|---|
| A. One 2-D vapor field with five sinks | Only global sink competition; vertical shielding unresolved | Distinct z layers share the same c(x,y), can be force-allocated arbitrarily | Low | Keep as null/control model only |
| **B. Five ice-eligible slabs + upper/lower gas slabs** | **Yes, explicit adjacent gas-cell inventory and shared face fluxes** | Coarse z discretization, incomplete 3-D diffusion and unresolved microscopic terraces | Medium/high | **Provisional selected reference** |
| C. Continuous height field + effective vertical diffusion resistance | Partial, closure-dependent | May silently assign central preference; difficult to falsify without B | Low/medium | Potential reduced surrogate after benchmarking against B |

**Selected:** a finite-volume **2D basal plane × few z planes** gas/ice system (a *coarse 3-D transport discretization*, not a general 3-D surface solver). The *rendered* number of ice layers is five. Gas slabs are separate and extend on both sides. Begin with five ice-eligible slabs plus two gas-only slabs above and two below (9 total) as a **test candidate**; vary gas halo length, Δz and slab count to check convergence before accepting them. In the short term the original hex25d.js stays the truth in legacy mode.

## Geometry and conserved inventory

At basal spacing `a=0.75e-6 m` with axial (q,r) sites:

- Hex Voronoi area `Ahex = sqrt(3)/2 * a*a` [m²].
- Shared in-plane face length `Lface = a/sqrt(3)` [m].
- In-plane gas face area `Axy = Lface * Δz_gas` [m²], centre spacing `a` [m].
- Vertical gas face area `Az = Ahex` [m²], centre spacing `Δz_gas` (or half-width sums).
- Gas volume `V_g = Ahex * Δz_gas * openFraction` [m³].
- Ice-volume contribution per slab `V_ice = Ahex * Δz_ice * phiIce` [m³].
- All concentrations `c` are **water-vapor mass concentration in kg/m³**, not current legacy excess `sigma`. Convert only at the solver boundary, using validated ice saturation density.
- Vapour mass in gas CV `m_g = c * V_g` [kg]; ice mass `m_i = rhoIce * V_ice` [kg].

For neighboring gas control volumes i,j use one single **signed face flux** (kg/s):

    F_ij = D_v * A_ij / l_ij * (c_i - c_j).

Apply **-F_ij dt to i, +F_ij dt to j**, once per face. Lateral/vertical open face areas must be computed from gas accessibility / ice obstruction with a documented transmissibility rule. Never use a separate unaccounted sink for each optical layer.

At exposed gas/ice interfaces compute an attachment-limited deposition flux (kg/s) with local equilibrium concentration c_eq(T,κ,normal,facet/SDAK). Remove deposited mass from *adjacent gas inventory* and add it to *ice volume* **exactly once**. If combined candidate sinks exceed gas inventory, cap/redistribute them together, not independently layer by layer.

For latent heat: each kg deposited releases approximately `Ls` joules into the thermal balance (not energy taken from diffusion). Transfer heat through exposed ice/gas faces, account outer heat fluxes. The legacy heat field `dT` is relative to the far-field temperature; do not reinterpret it as joules.

Per-step checks (zero closed-boundary flux):

    Σ gasMass_before + Σ iceMass_before =
    Σ gasMass_after  + Σ iceMass_after  + residual.

For open boundaries report `boundaryVaporMassIn`, `boundaryVaporMassOut` explicitly. Similarly report energy changes, deposited latent-energy redistribution and external heat exchange. **Zero negative inventories**.

The current hex25d.js `totalDepositedMass` counts in-plane deposition only. Its independently updated scalar `basalHalfThickness` **is not** included in that mass or latent heat. Treat that as a known legacy deficiency, not a validated mass conservation result. Do not retroactively compute 3-D mass from old 2.5-D frames and label it conserved.

## Physical alternatives / falsification

The model should **not** hardcode fast central-layer growth. Test all of:

- Symmetric seed/flux and identical coefficients: reflective upper/lower symmetry survives (up to explicitly quantified error).
- If upper/lower basal surfaces are well exposed, they may dominate uptake and cause central **shielding**, rather than central-spine growth. This is a valid result.
- Central slab advantage only occurs if flux geometry and *measured/local* prism facet attachment rates genuinely favor its exposed lateral tips.
- 5-ice + {1,2,4} exterior gas slabs on each side, Δz refinement, basal-domain expansion and transport convergence.
- Turn off SDAK, surface transport, temperature feedback independently.
- Track `massFluxByLayer[5]` and `heatByLayer[5]`, not fabricated 42%-17%-... flux percentages.
- Compare grown-arm spine occupancy and independently evolved slab occupancy rather than transparent duplicate contours.
- For an unsupported case (no gas around a buried ice cell), deposition **must be zero**. Never supply material from inside solid ice.
- Compare against Model A as an explicit controlled approximation; if B offers no identifiable benefit, revisit cost.

## Computational budget (approximate, not benchmarked)

At N=401 the current hex footprint has **120,601** active basal cells. Nine gas/ice-eligible vertical slabs imply ~**1,085,409** active slab-cells. One Float64 field across these is ~**8.7 MB**, before padding and JS object overhead; 10 such arrays are ~87 MB. The naive snapshot history is not acceptable; store sparse, delta-compressed or checkpointed states and budget the complete browser heap/WebGL memory. A fast conservative stencil, independent iteration residual and convergence criterion are mandatory: `45 vaporIterations` alone does not establish convergence on 300 µm.

First benchmark a small 41²/73² hex domain and a short run. Increase N and number of slabs only after stable closed-boundary tests and interactive timing data. Consider workers/WASM and multigrid only after measurement.

## Implementation seam and versioning

Existing owners, **unchanged until step 8 of the 15-point plan**:

- `app/solver/hex25d.js`: legacy six-neighbor phase/vapor/heat solver, basal scalar, original tests.
- `app/solver/growth25d.js`: Cartesian reference; electric Poisson coupling remains Cartesian only.
- `app/solver/ice-snapshot-v1.js`: version 1 explicitly says `legacy-uniform-slab`, `localHeightSolved:false`, `verticalVaporTransportSolved:false`, `basalMassCoupled:false`.
- `app/render/uniform-five-layers.js`: **preview**: 5 identical xy contours at z=-h,-h/2,0,h/2,h.
- `app/observatory.js`: orchestration and render-only toggling, never source of independently computed ice mass.

Future v2 `five-layer-conservative` uses independent arrays on the solver side, not five independent instances of `hex25d.js` consuming the same far-field budget. New modules may include `app/solver/layered-gas.js`, `layered-ice.js`, `layer-fluxes.js`, `layered-temperature.js`; names and field contract are provisional. Version bump only after conservation tests. Legacy v1 remains selectable in viewer for A/B comparisons.

## Milestone 1 of 15: freeze original before new physics

Pin immutable **commit 77f118c091a2fd1ae278c899e7d5f6a1715ebc54** and baseline branch `baseline/five-layer-pre-coupling-20261015`. Preserve the already implemented Pass 7 snapshot/5-layer preview rather than pretending it needs to be redone.

Create an independent deterministic runner measuring the existing solver's exact fields and budgets; record the reference results and a test for them in the development branch. Do not modify `hex25d.js`, `growth25d.js`, the rendered five preview planes, time step, initial conditions, CSV schema, or Pages source branch as part of this milestone.

Next milestone after baseline acceptance: inspect and validate the numerical formulation of a 5-ice-slab + gas-halo finite-volume cell with closed mass/energy test, **before** any viewer-dependent z-layer growth.
