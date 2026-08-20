# Milestone 2.5 — Seed to Coupled Ice Growth

Milestone 2.5 joins the previously separate nucleation/embryo work to the first continuously growing ice object.

## 1. Sub-grid postcritical bridge

A CNT-sized postcritical embryo is nanometre-scale and should not be dropped directly onto a micron-scale field grid. `app/solver/postseed.js` therefore grows the compact seed spherically until it reaches a chosen resolved radius.

The bridge solves a self-consistent surface state using:

- vapor diffusion,
- attachment kinetics,
- latent-heat warming,
- Ice-Ih density from the IAPWS EOS.

The coupled relations are

`rho_i v = D_v (c_inf - c_s) / r`

`v = alpha v_kin sigma_s`

`T_s = T_inf + rho_i L_s v r / k_air`.

The current vapor enthalpy contribution to sublimation enthalpy is still a compact atmospheric approximation; the condensed ice state is IAPWS Ice Ih.

## 2. Resolved basal-plane vapor field

Once `r >= resolvedRadius`, the same seed is handed to `app/solver/growth25d.js`.

The basal plane is represented by a Cartesian finite-volume grid with an evolving ice volume fraction. The external vapor field is quasi-steady:

`nabla^2 sigma = 0`.

The outer boundary is fixed at `sigma_inf`.

The warmed ice boundary uses the first-order Clausius–Clapeyron shift

`sigma_eq(T_s) ~= L_s DeltaT / (R_v T_inf^2)`.

The interface deposition drive is therefore the gas-side supersaturation relative to the locally warmed ice equilibrium, not raw ambient supersaturation.

## 3. Interface mass transfer

Each interface face uses diffusion and attachment kinetics as series resistances:

`j = Delta sigma / [dx/(D_v c_sat,i) + 1/(rho_i alpha v_kin)]`.

Mass flux advances the local ice volume fraction. A cell becomes resolved ice after its ice fraction approaches unity.

This gives an explicit local budget for:

- `Delta sigma`,
- attachment coefficient,
- mass flux,
- deposited mass,
- interface cell count,
- latent heat.

## 4. Heat field

Deposition creates a local volumetric heat source

`q''' = m_dot L_s / V_cell`.

The quasi-steady external heat field is solved from

`nabla^2 DeltaT = -q''' / k_air`

with `DeltaT = 0` on the outer boundary.

The resulting temperature field feeds back into the local ice-equilibrium vapor condition before the next growth step.

## 5. Attachment kinetics and sixfold response

The current facet law is

`alpha = A exp(-sigma0 / sigma_surface)`.

Separate `sigma0` controls are carried for prism and basal attachment.

The resolved basal plane includes an explicit sixfold modulation of prism attachment. A switchable SDAK research mode lowers local prism `sigma0` near under-coordinated edge cells. That SDAK option is deliberately labeled experimental and is not yet a calibrated morphology law.

The dense radial morphology diagnostic computes the actual `m=6` amplitude from 120 angular samples, avoiding the six-ray `m=6 <-> m=0` aliasing problem identified earlier.

## 6. Why this is 2.5-D

The basal plane is spatially resolved. The out-of-plane direction is not yet a 3-D mesh.

Instead the model carries an effective basal half-thickness `H` and advances it with basal attachment kinetics using the mean resolved interface supersaturation.

Therefore the live morphology state contains:

- resolved basal radius / contour,
- effective basal half-thickness,
- aspect ratio `R/H`,
- dense `m=6` amplitude.

This is enough to begin plate/column competition without misrepresenting the model as a full 3-D PDE.

## 7. Interactive laboratory

Open `app/milestone25.html`.

It displays in one synchronized state:

- vapor supersaturation field,
- latent-heat field,
- evolving basal-plane contour,
- sixfold diagnostic orientation,
- 2.5-D side profile,
- prism and basal attachment coefficients,
- deposition rate,
- interface supersaturation,
- maximum local temperature rise,
- accumulated latent energy,
- automatic sub-grid -> resolved handoff event.

Controls expose temperature, far-field ice supersaturation, postcritical seed radius, resolved handoff radius, prism/basal `sigma0`, sixfold anisotropy, field-relaxation depth, and an experimental SDAK switch.

## 8. Scientific boundary

Milestone 2.5 now contains a continuous chain

`supercooled droplet -> stochastic nucleation -> CNT postcritical seed -> sub-grid growth -> resolved vapor/heat growth`.

It still does **not** claim:

- a full 3-D vapor/temperature PDE,
- resolved basal facets in a volumetric mesh,
- convection/advection,
- calibrated SDAK,
- dislocation or stacking-disorder kinetics,
- material-specific immersion freezing,
- executable PCF,
- electric-field or flexoelectric feedback.

Those become the next layers rather than hidden knobs in 2.5.

## 9. PC testing

Serve the repository over a local HTTP server and open:

`app/labs.html`

For the full numerical regression bundle run:

`node app/tests/run-all.mjs`

The Milestone-2.5-specific regressions are:

- `app/tests/postseed.test.mjs`
- `app/tests/growth25d.test.mjs`

## 10. Main source references

- IAPWS R10-06(2009), Equation of State for H2O Ice Ih.
- IAPWS G12-15, Thermodynamic Properties of Supercooled Water.
- Murphy & Koop (2005), vapor pressure of ice and supercooled water.
- Koop et al. (2000), water activity and homogeneous ice nucleation.
- Libbrecht snow-crystal growth literature for attachment-kinetics form `v = alpha v_kin sigma` and nucleation-limited `alpha` parameterizations.
