# Physics Instrument v2 — Activation + Supercool Boundary Laboratory

Open `app/index.html` for the activation laboratory and `app/supercool.html` for the supercooled-liquid handoff. The legacy single-file flexoelectric snowflake remains at the repository root while the physics-first rebuild advances independently.

## Milestone 1 activation laboratory

One synchronized particle state drives:

- full Petters–Kreidenweis κ-Köhler equilibrium curve,
- numerical critical diameter and critical supersaturation,
- stable and unstable wet-equilibrium branches,
- state-plane trajectory `(wet diameter, ambient supersaturation)`,
- molar transfer chemical-potential diagnostic,
- reversible-work landscape,
- Maxwell–Mason spherical condensational growth with vapor-diffusion and latent-heat resistances,
- controlled ambient histories,
- scrub-able particle trajectory/history,
- spherical sixfold HexaLens + dense 72-angle mode spectrum,
- dry-core / water-shell particle visualization,
- inspectable activation / transport budget.

## Supercooled-liquid boundary laboratory

`app/supercool.html` continues the *same activated droplet* below 0 °C rather than replacing it with an unrelated ice model.

It now exposes separately:

- saturation ratio with respect to liquid water `S_w`,
- saturation ratio with respect to ice `S_i`,
- κ-derived solution water activity `a_w`,
- ice-equilibrium water activity `a_w^i(T) = p_i,sat(T)/p_w,sat(T)`,
- `Δa_w = a_w - a_w^i(T)`,
- liquid-to-ice chemical-potential drive `Δμ = R T ln[a_w/a_w^i(T)]`,
- the liquid-stable / supercooled-metastable regime distinction,
- a T–a_w phase-state surface,
- fixed-vapor cooling trajectory,
- separate history events for sub-zero entry, ambient ice supersaturation, liquid metastability, and entry into the classic Koop water-activity coordinate window,
- exact conversion from an externally supplied volumetric nucleation rate `J` to one-particle freezing probability `P = 1-exp(-J V Δt)`.

### Important nucleation boundary

The classic Koop et al. water-activity result is used as a **coordinate system**, not as an undocumented rate law. The current code highlights the traditional `Δa_w` interval 0.26–0.34 and conservatively marks the original deep-supercooling temperature domain, but does not invent a universal `J(Δa_w)` polynomial.

An empirical homogeneous-nucleation-rate parameterization will be added only as a named, documented Milestone-2 plug-in. Heterogeneous ice-nucleating particles, immersion freezing, pore condensation/freezing, and embryo geometry also remain deferred.

## Scientific categories

### Thermodynamic / reviewed-fit

- Murphy & Koop saturation vapor pressure over liquid water and ice.
- IAPWS surface tension of liquid water.
- Kelvin curvature term.
- liquid/ice equilibrium water activity from the ratio of ice and liquid saturation vapor pressures.

### Semi-empirical

- Petters–Kreidenweis single-parameter κ representation of water activity / CCN activity.
- Koop water-activity organization of homogeneous ice nucleation is used as a diagnostic coordinate, not yet as a rate parameterization.

### Transport approximation

- Maxwell–Mason spherical condensational response.
- Temperature/pressure-scaled vapor diffusivity and compact air thermal-conductivity / latent-heat approximations.
- Current external vapor field remains quasi-steady and spherical.

The ambient-history drivers are controlled-reservoir experiments. `cool_fixed_e` cools while holding vapor partial pressure fixed. It is intentionally **not** a complete adiabatic cloud-parcel model because parcel water conservation, vertical dynamics, aerosol population feedback, and supersaturation depletion by a droplet ensemble are not yet solved.

## Stability and reversible work

For fixed ambient saturation ratio `S_inf`, roots of

`S_eq(D) = S_inf`

are found numerically. Linearizing `dD/dt ∝ S_inf - S_eq(D)` gives a stable root for positive `dS_eq/dD` and an unstable activation threshold for negative slope. The critical point is the maximum of the full κ-Köhler curve; when `S_inf >= S_crit`, the finite-size activation barrier disappears.

The reversible-work display uses

`dW = R T ln(S_eq/S_inf) dN_w`

with spherical liquid-water content. This is a one-particle controlled-reservoir coordinate, not a molecular ice-nucleation barrier.

## Important sixfold symmetry note

The six arm axes are ideal for the visible HexaLens radial comparison, but they are not sufficient to estimate an `m=6` Fourier mode. At `theta_i=i*pi/3`, `cos(6 theta_i)=1` on every arm, so `m=6` aliases `m=0`. The central mode spectrum therefore samples 72 angular directions.

Milestone 1 remains intentionally spherical. Non-spherical vapor/heat/electric fields must earn any later symmetry breaking.

## Tests

Run:

```bash
node app/tests/thermo.test.mjs
node app/tests/activation.test.mjs
node app/tests/trajectory.test.mjs
node app/tests/supercool.test.mjs
```

The suites cover saturation-pressure checkpoints, IAPWS surface tension, numerical κ-Köhler critical behavior, stable/unstable equilibria, reversible-work barrier, Maxwell–Mason resistances, ambient drivers, activation trajectories, liquid/ice water-activity equilibrium, supercooled metastability, ice supersaturation, phase-surface construction, and exact Poisson nucleation-hazard algebra.

## Scientific references

- Petters, M. D. & Kreidenweis, S. M. (2007), *A single parameter representation of hygroscopic growth and cloud condensation nucleus activity*, Atmospheric Chemistry and Physics 7, 1961–1971. DOI: 10.5194/acp-7-1961-2007.
- Murphy, D. M. & Koop, T. (2005), *Review of the vapour pressures of ice and supercooled water for atmospheric applications*, QJRMS 131, 1539–1565. DOI: 10.1256/qj.04.94.
- Koop, T., Luo, B., Tsias, A. & Peter, T. (2000), *Water activity as the determinant for homogeneous ice nucleation in aqueous solutions*, Nature 406, 611–614. DOI: 10.1038/35020537.
- IAPWS (2014), *Revised Release on Surface Tension of Ordinary Water Substance*, R1-76(2014).
- IAPWS (2015), *Guideline on Thermodynamic Properties of Supercooled Water*, G12-15.
