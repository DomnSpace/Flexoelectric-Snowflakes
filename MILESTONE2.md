# Milestone 2 — Nucleation and Ice Seed

Milestone 2 starts from the same supercooled, composition-resolved droplet state built in Milestone 1. It does **not** jump directly to a six-armed crystal.

## Executable homogeneous pathway

For homogeneous freezing, the current implementation uses the Koop water-activity fit in the form

`log10 J_hom[cm^-3 s^-1] = -906.7 + 8502 Δaw - 26924 Δaw^2 + 29180 Δaw^3`

with `Δaw = aw - aw_i(T)` and fit window `0.26 < Δaw < 0.34`.

The solver refuses to silently use the fit outside this range unless extrapolation is explicitly requested.

Given volumetric rate `J`, the one-droplet hazard rate is

`lambda = J V_liquid`

and the cumulative hazard is

`H(t)=integral lambda dt`.

A seeded stochastic clock samples one exponential threshold

`H* = -ln(U)`.

Nucleation occurs when `H >= H*`. Rendering never advances this clock; only physical time integration does.

## Immersion pathway

The code exposes an activity-based immersion-freezing form

`log10 J_het[cm^-2 s^-1] = m Δaw + c`

but `m,c` are material-specific and therefore no universal default is supplied. An immersed surface area and a documented material calibration are required before this pathway becomes executable.

## PCF pathway

Pore-condensation freezing remains a real pathway interface but is not yet executable. It requires at minimum:

- pore radius / geometry distribution,
- contact angle / wetting assumptions,
- capillary condensation condition,
- pore-water freezing condition,
- ice escape / growth condition.

No generic PCF probability is invented.

## Seed birth

A stochastic nucleation event creates a compact ice-seed state with:

- total H2O mass,
- ice mass,
- remaining liquid-water mass,
- ice mass fraction,
- equivalent compact-ice volume/diameter,
- remaining liquid composition,
- latent-heat pulse,
- mass residual diagnostic,
- nucleation pathway and timestamp.

Mass conversion is exact by construction. The current compact seed still uses temporary compact approximations for ice density and latent heat. The full IAPWS Ice Ih EOS and temperature-dependent phase energetics are the next refinement.

## Sixfold orientation frame

The seed can expose six prism-normal directions plus two basal normals as an **orientation frame only**. No arm growth, dendrites, SDAK, or field-generated branching occurs in Milestone 2 seed birth.

## Still deferred before post-seed crystal growth

- full IAPWS Ice Ih thermodynamics,
- ice–liquid interfacial free energy `gamma_il`,
- critical embryo radius / barrier from classical or non-classical nucleation theory,
- heterogeneous active-site calibration,
- executable PCF pore model,
- latent-heat redistribution into the coupled temperature field,
- non-spherical vapor + heat solver,
- basal/prism attachment kinetics,
- SDAK,
- 3-D faceted crystal interface sourced from the solver.

## Files

- `app/solver/nucleation.js` — rates and stochastic hazard clock
- `app/solver/ice-seed.js` — mass-conserving seed conversion
- `app/milestone2.html` — interactive phase-change lab
- `app/tests/nucleation.test.mjs` — rate/hazard/mass regression tests

## Primary references

- Koop, T., Luo, B., Tsias, A. & Peter, T. (2000), *Water activity as the determinant for homogeneous ice nucleation in aqueous solutions*, Nature 406, 611–614.
- Knopf, D. A. & Alpert, P. A. (2013), activity-based immersion freezing framework for atmospheric particles.
- IAPWS R10-06(2009), Revised Release on the Equation of State 2006 for H2O Ice Ih.
