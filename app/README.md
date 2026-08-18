# Physics Instrument v2

Open `app/index.html` through a local/static HTTP server. The legacy single-file flexoelectric snowflake prototype remains at the repository root while the physics-first rebuild is validated.

## What is real in Milestone 1

- Murphy & Koop saturation vapor pressure over liquid water and ice.
- IAPWS surface tension of liquid water.
- Full Petters–Kreidenweis κ-Köhler equilibrium relation.
- Numerical search of the full κ-Köhler critical point; the common analytic critical-point formula is retained only as a regression diagnostic.
- Quasi-steady spherical vapor-diffusion field around the wet particle.
- Vapor-diffusion-only condensational radius tendency. Latent-heat resistance is deliberately deferred to the coupled vapor/heat field solver so it is not hidden in a magic coefficient.
- Shared particle state and event history.
- Six crystallographic radial profiles plus a dense 72-angle Fourier mode decomposition.
- Inspectable activation budget in log-saturation / chemical-potential coordinates.

## Important symmetry note

The six arm axes are ideal for the HexaLens radial comparison, but they are **not sufficient to estimate an m=6 Fourier mode**: at θᵢ=iπ/3, cos(6θᵢ)=1 for every arm, so m=6 aliases the isotropic m=0 component. Therefore the central mode spectrum samples a dense angular ring while the visible instrument remains sixfold.

## Current boundary

The current field is intentionally spherical and neutral. Therefore all six HexaLens profiles coincide and residual angular modes are approximately zero. This is a physics test, not a missing visual effect.

The following are reserved for later milestones and are not faked in Milestone 1:

- freezing / pore-condensation-freezing / ice nucleation,
- vapor + heat PDE around a non-spherical ice interface,
- basal/prism attachment kinetics and SDAK,
- growing 3-D crystal sourced from the interface solver,
- interface-click physical budget,
- electric potential and flexoelectric polarization feedback.

## Tests

Run:

```bash
node app/tests/thermo.test.mjs
```

The regression suite checks 0 °C liquid/ice saturation pressures, 25 °C surface tension, numerical κ-Köhler critical-point behavior, agreement with the large-particle approximation within a tolerance, and zero net activation budget at the numerical critical point.

## Scientific references

- Petters, M. D. & Kreidenweis, S. M. (2007), *A single parameter representation of hygroscopic growth and cloud condensation nucleus activity*, Atmospheric Chemistry and Physics 7, 1961–1971. DOI: 10.5194/acp-7-1961-2007.
- Murphy, D. M. & Koop, T. (2005), *Review of the vapour pressures of ice and supercooled water for atmospheric applications*, QJRMS 131, 1539–1565. DOI: 10.1256/qj.04.94.
- IAPWS (2014), *Revised Release on Surface Tension of Ordinary Water Substance*, R1-76(2014).
