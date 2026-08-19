# Physics Instrument v2 — Activation Laboratory

Open `app/index.html` through a local/static HTTP server. The legacy single-file flexoelectric snowflake remains at the repository root while the physics-first rebuild advances independently.

## Milestone 1 is now an activation laboratory

The app no longer stops at drawing a κ-Köhler curve. One synchronized particle state now drives:

- the full Petters–Kreidenweis κ-Köhler equilibrium curve,
- numerical critical diameter and critical supersaturation,
- stable and unstable wet-equilibrium branches for the current ambient saturation,
- a state-plane trajectory `(wet diameter, ambient supersaturation)`,
- a molar transfer chemical-potential diagnostic `Δμ(l−v)=RT ln(S_eq/S_inf)`,
- a reversible-work landscape obtained by integrating `Δμ dN_w` along the wet-diameter coordinate,
- Maxwell–Mason spherical condensational growth with both vapor-diffusion and latent-heat resistances,
- controlled ambient histories (hold, supersaturation ramp, and cooling at fixed vapor partial pressure),
- a scrub-able particle trajectory/history,
- spherical vapor-field radial sampling in the sixfold HexaLens,
- a dense 72-angle Fourier mode spectrum that correctly distinguishes `m=6` from `m=0`,
- a dry-core / water-shell particle visualization in the HexaLens center,
- an inspectable local activation / transport budget.

## Scientific boundary

### Thermodynamic / reviewed-fit

- Murphy & Koop saturation vapor pressure over liquid water and ice.
- IAPWS surface tension of liquid water.
- Kelvin curvature term.

### Semi-empirical

- Petters–Kreidenweis single-parameter κ representation of water activity / CCN activity.

### Transport approximation

- Maxwell–Mason spherical condensational response.
- Temperature/pressure-scaled vapor diffusivity and compact air thermal-conductivity / latent-heat approximations.
- The current external vapor field is quasi-steady and spherical.

The ambient-history drivers are controlled-reservoir experiments. In particular, `cool_fixed_e` cools while holding vapor partial pressure fixed. It is intentionally **not** called a complete adiabatic cloud-parcel model because it does not yet solve parcel water conservation, vertical dynamics, aerosol population feedback, or supersaturation depletion by a droplet ensemble.

## Stability and reversible work

For a fixed ambient saturation ratio `S_inf`, roots of

`S_eq(D) = S_inf`

are found numerically. Linearizing `dD/dt ∝ S_inf - S_eq(D)` gives:

- positive `dS_eq/dD`: stable wet equilibrium,
- negative `dS_eq/dD`: unstable activation threshold.

The critical point is the maximum of the full κ-Köhler curve. When `S_inf >= S_crit`, the finite-size activation barrier disappears.

The reversible-work display uses

`dW = R T ln(S_eq/S_inf) dN_w`

with spherical liquid-water content

`dN_w/dD = rho_w π D² / (2 M_w)`.

This is a one-particle thermodynamic coordinate for the controlled reservoir, not a molecular ice-nucleation barrier and not a stochastic activation probability.

## Important sixfold symmetry note

The six arm axes are ideal for the visible HexaLens radial comparison, but they are **not sufficient to estimate an m=6 Fourier mode**. At `theta_i=i*pi/3`, `cos(6 theta_i)=1` on every arm, so `m=6` aliases the isotropic `m=0` component. The central mode spectrum therefore samples 72 angular directions while the visible instrument remains sixfold.

Milestone 1 is intentionally spherical. All six HexaLens profiles should coincide and all nonzero angular residual modes should be numerical zero. Later non-spherical vapor/heat/electric fields earn the right to break that symmetry.

## Current boundary

Still reserved and not faked:

- supercooled-droplet freezing / immersion / pore-condensation-freezing pathways,
- ice nucleation and mass-conserving ice-seed generation,
- non-spherical coupled vapor + heat field solver,
- basal/prism attachment kinetics and SDAK,
- growing 3-D crystal sourced from the interface solver,
- click-any-facet physical growth budget,
- electric potential and flexoelectric polarization feedback.

## Tests

Run:

```bash
node app/tests/thermo.test.mjs
node app/tests/activation.test.mjs
node app/tests/trajectory.test.mjs
```

The new tests cover:

- Murphy–Koop saturation-pressure checkpoints,
- IAPWS surface tension,
- numerical κ-Köhler critical point,
- stable/unstable equilibrium branches below critical saturation,
- disappearance of the finite-size barrier above `S_crit`,
- positive reversible-work barrier below critical,
- positive vapor and latent-heat transport resistances,
- controlled ambient-driver behavior,
- an end-to-end supersaturation-ramp trajectory that crosses both `S_crit` and `D_crit` and grows into the micron regime.

## Scientific references

- Petters, M. D. & Kreidenweis, S. M. (2007), *A single parameter representation of hygroscopic growth and cloud condensation nucleus activity*, Atmospheric Chemistry and Physics 7, 1961–1971. DOI: 10.5194/acp-7-1961-2007.
- Murphy, D. M. & Koop, T. (2005), *Review of the vapour pressures of ice and supercooled water for atmospheric applications*, QJRMS 131, 1539–1565. DOI: 10.1256/qj.04.94.
- IAPWS (2014), *Revised Release on Surface Tension of Ordinary Water Substance*, R1-76(2014).
