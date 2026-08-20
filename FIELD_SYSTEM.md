# Canonical spatial field system — Milestone 2.5

The flake is governed by fields, not by decorative radial rules. The current resolved solver is a 2.5-D basal-plane finite-volume model with effective out-of-plane basal thickness.

## 1. Stored vapor scalar

The grid stores

`u(x,y) = c_v(x,y) / c_sat,i(T_inf) - 1`

rather than local ice supersaturation. For constant vapor diffusivity, quasi-steady diffusion is

`laplacian(u) = 0`

in gas cells.

The physical local ice supersaturation is reconstructed with the local temperature and curvature-adjusted equilibrium concentration.

## 2. Temperature scalar

The heat field stores

`Delta T(x,y) = T(x,y) - T_inf`.

With deposited latent heat represented as a volumetric source,

`laplacian(Delta T) = -q''' / k_air`.

Far-field temperature perturbation is zero.

## 3. Vapor mass-flux vector

`J_v = -D_v grad(c_v)`

and, because `c_v = c_ref (1+u)`, the implemented constant-property form is

`J_v = -D_v c_ref grad(u)`.

Units: `kg m^-2 s^-1`.

Streamlines show where water vapor is actually supplied to the crystal.

## 4. Conductive heat-flux vector

`q_h = -k_air grad(T)`.

Units: `W m^-2`.

These streamlines point away from warm latent-heat-producing regions toward the thermal reservoir.

## 5. Coupled growth-potential scalar

The main combined diagnostic is

`chi = ln[c_v / c_eq,i(T)]`.

Without capillarity this is equivalent to `ln(1 + sigma_i,local)`.

`chi > 0` means local vapor supply exceeds the local thermally shifted ice equilibrium concentration.

The diagnostic vector field is

`G_chi = -grad(chi)`.

It is not an additional physical force. It is a compact map of the net diffusion/thermal growth landscape.

## 6. Interface normal

The prism-interface orientation is computed from the local phase-fraction gradient using Sobel derivatives.

The outward normal is

`n = -grad(phi) / |grad(phi)|`.

No radial direction from the domain center is used for facet orientation.

This matters once branches or asymmetric perturbations appear.

## 7. Curvature / Gibbs–Thomson shift

In-plane curvature is reconstructed from derivatives of the phase fraction.

The local ice-equilibrium vapor density is corrected as

`c_eq = c_sat,i(T_surface) exp(delta_iv * kappa)`

with

`delta_iv = gamma_iv M_w / (rho_i R T)`.

The current `gamma_iv` is an explicit source-tagged reference parameter, not a hidden tuning constant.

## 8. Attachment anisotropy

Prism attachment uses a nucleation-limited form

`alpha_prism = A exp(-sigma0_prism / sigma_drive)`

multiplied by the explicit sixfold response

`1 + epsilon6 cos(6 theta_n)`.

Basal attachment is carried separately to evolve the out-of-plane half-thickness.

The optional SDAK switch modifies the local prism nucleation barrier near low-coordination edge cells and remains experimental/un-calibrated.

## 9. Interface deposition

For a local gas cell adjacent to ice, deposition uses combined diffusion and attachment resistances. The driving concentration is

`c_gas - c_eq(T_surface,kappa)`.

The deposited mass updates the phase fraction, and its latent heat becomes the source for the next heat solve.

Thus one growth iteration is

`heat relax -> vapor relax -> local equilibrium -> attachment/deposition -> latent source -> phase update`.

## 10. Numerical isotropy

Scalar vapor and heat relaxation use a 9-point isotropic Laplacian rather than the basic five-point stencil.

This is intended to suppress square-grid selection. The morphology spectrum therefore reports at least

`m = 1,2,3,4,6,8,12`.

`m4` is a mesh-contamination diagnostic; `m6` is the intended crystallographic mode. Six visible rays are never used to estimate `m6`.

## 11. Field-line integration

Streamlines use RK2 midpoint integration over bilinearly interpolated scalar gradients.

Available line families:

- vapor supply `u`,
- thermal field `Delta T`,
- combined growth potential `chi`.

Streamlines never advance growth and never alter the fields.

## 12. Residual diagnostics

`fieldlines.js` reports numerical checks including

- RMS divergence of vapor mass flux in gas cells,
- RMS divergence of conductive heat flux,
- RMS curl of the gradient-derived `chi` vector field,
- maximum vapor flux,
- maximum heat flux,
- maximum growth-potential gradient,
- interface inward vapor and outward heat flux summaries.

Residuals are not expected to be machine zero on a moving, source-driven finite grid, but they should converge downward with increased relaxation and grid refinement.

## 13. What still does not exist

The current neutral field system does not yet contain

- resolved 3-D vapor/temperature PDEs,
- airflow/advection,
- electrostatic field `E = -grad(phi_E)`,
- polarization/flexoelectric bound charge,
- mechanical strain field,
- anisotropic thermal/vapor transport tensors,
- surface premelting or explicit quasi-liquid layer,
- fully calibrated SDAK.

Those should be added as fields with their own conservation/equilibrium equations, not as visual multipliers on attachment.

## Relevant files

- `app/solver/growth25d.js` — evolving phase, vapor and heat solver
- `app/solver/fieldlines.js` — flux vectors, streamlines and residuals
- `app/spatialfield.html` — canonical field visualization
- `app/tests/growth25d.test.mjs`
- `app/tests/fieldlines.test.mjs`
