# Snow-crystal habits & numerical-domain validity

This pass is driven by four user observations: the old model could display a **134.52 µm tall extruded column inside a ~55 µm basal-plane domain**, the vapor map displayed strongly square boundary structure, the former "dendrite" preset did not guarantee branches, and the only domain guard used equivalent-area radius.

## A. Experimental basis versus solver implementation

The temperature/supersaturation dependence is empirical: according to the **Nakaya morphology diagram**, around −15 °C plate-like crystals occur at low ice supersaturation, sectored plates at intermediate humidity and stellar dendrites at higher humidity. Around −5 °C higher humidities favor columns/needles instead.

Reference: Kenneth G. Libbrecht, *The physics of snow crystals*, Reports on Progress in Physics **68** (2005), 855–895, especially the habit diagram in Fig. 2. Author PDF: https://www.its.caltech.edu/~atomic/publist/rpp5_4_R03.pdf . Real crystal morphology is stochastic, history- and ventilation-dependent.

The three implemented scenarios in `app/solver/habits25d.js` share the **same resolved growth solver**. They differ only in exposed initial/environmental and facet-kinetic coefficients (including basal `A`). They do **not** load target crystal polygons, overlays, or synthetic measurements. The facet coefficients and SDAK values are *illustrative, uncalibrated parameters*, not inverse-fitted to lab experiments, so a "stellar candidate" is a research hypothesis, not evidence a stellar dendrite emerges.

| Scenario | Temperature | Far-field ice supersaturation | Grid | Status |
|---|---:|---:|---:|---|
| Compact plate | −15 °C | 3% | 97² / 72 µm | Testable 2.5-D plane |
| Sectored plate candidate | −15 °C | 9% | 129² / 96 µm | Testable; sectors NOT guaranteed |
| Stellar dendrite candidate | −15 °C | 14% | 129² / 96 µm | Testable; six arms NOT guaranteed |
| Column / needle reference | −5 °C | 18% | — | Disabled: 3-D model required |

At −15 °C, the ice supersaturation at *liquid-water saturation* is calculated, rather than guessed, using the existing Murphy–Koop vapor-pressure formulas, as

    100 × (e_sat,water(T) / e_sat,ice(T) − 1).

The UI shows this **water-saturation reference** and flags values above it, which may be permissible transiently in the laboratory but are not an ordinary steady mixed-phase cloud prescription. Pressure defaults to about 101 kPa. This is not an atmospheric parcel, cloud-flow, nucleation or weather forecasting model.

## B. Why the 134 µm tall ice column was not a snowflake

The original solver computed out-of-plane basal growth using a fixed basal nucleation coefficient (A = 1) multiplied by local mean supersaturation, whereas prism growth is diffusion-limited in the basal plane. The old viewer extruded that *shared* basal thickness literally. No 3-D vapor/heat transport or basal-surface mass conservation existed to justify such a tall object.

The revised interface **exposes basal attachment prefactor A** explicitly, with 0.009 as the initial explorative value and lower illustrative values in the plate/dendrite scenarios. This parameter *enters the existing basal attachment kinetic equation* rather than visually squashing a valid tall crystal. Thin-plate assumptions are also checked continuously. **Do not call the chosen A measured** or a complete correction to out-of-plane physics.

The current mass and latent-heat totals account for **in-plane prism deposition only**. Evolving the basal thickness does not update that mass/energy source, a known 2.5-D modeling deficiency. This limitation is explicitly labeled in the viewer and must be fixed by a resolved 3-D model or a validated, conserved reduced-order basal exchange scheme before any quantitative cloud-physics claim.

## C. Finite-domain audit

`auditGrowthDomain` looks at **all** occupied phase cells (`ice = 1` or `phi >= 0.15`), rather than equivalent-area radius. On the square N×N grid at spacing dx:

    physical width = (N − 1) dx
    haloCells = min over occupied (i, j, N−1−i, N−1−j)
    tipRadius = max over occupied hypot((i−c)dx, (j−c)dx)
    thickness = 2 * basalHalfThickness
    plateRatio = basalHalfThickness / max(dx, tipRadius).

Warnings appear when halo ≤ **12 grid cells** or plateRatio > **0.65**. Growth pauses at halo ≤ **6 cells**, before the ice reaches the enforced Dirichlet vapor/thermal boundary. Growth also pauses after 5 s if plateRatio > **1.25**, because that tall object is outside the shared-thickness plate approximation. These **numerical heuristic thresholds** are declared and tested, not universal ice physics constants. A seed may initially be approximately isometric; the warnings do not retroactively assign it a plate habit.

**Crucial:** Domain growth only increases the distance to artificial boundary conditions. Larger N with unchanged dx **does not** improve interface resolution or validate branching physics. The square gradients in the field map are a consequence of finite square Dirichlet far-field boundaries; the visualization must display this fact rather than paint the square away.

Grid settings:

- 73² (54 µm), fast;
- 97² (72 µm);
- 129² (96 µm), default for detailed habit scenarios;
- 161² (120 µm), slower.

Replay storage is reduced as N increases to constrain browser memory. The crystal outline is still derived from true phase occupancy; the graphics are an extrusion of one basal-plane contour, **not a 3-D crystal, faceted basal-tip reconstruction, or resolved heat/vapor field in z**.

## D. Validation and open research tests

`node app/tests/habit-domain.test.mjs` checks preset consistency relative to water saturation, rejection of unsupported columns, physical domain geometry, tip-based edge detection even with tiny equivalent radius, tall-column warning/stop and actual basal kinetic sensitivity. The full suite is `node app/tests/run-all.mjs`. The Playwright smoke test confirms changing grid size really changes the field data, verifies the preset controls, and rejects unsupported 3-D habits.

**Still NOT demonstrated:** robust sixfold growth, primary dendrite tips, secondary arms, independence from box size, independence from dx, convergence of vapor/heat fields, measured growth rates, calibrated SDAK, or habit transitions. This pass cannot guarantee a snowflake just by selecting "dendrite." A proper next step is a convergence experiment at fixed physical domain with dx refinement, (m_4/m_6) analysis, anisotropic attachment calibration and explicit lateral-tip branching success criteria; true needles require 3-D diffusion, mechanics and facet kinetics.
