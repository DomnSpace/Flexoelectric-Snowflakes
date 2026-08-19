# Physics Instrument v2 — Activation → Supercooled Feedback Laboratory

Open `app/index.html` through a local/static HTTP server for the activation laboratory, then `app/supercool.html` for the metastable-liquid handoff. The legacy single-file flexoelectric snowflake remains at repository root while the physics-first rebuild advances independently.

## Milestone 1 activation laboratory

One synchronized particle state drives:

- full Petters–Kreidenweis κ-Köhler equilibrium curve,
- numerical critical diameter and supersaturation,
- stable and unstable wet-equilibrium branches,
- `(wet diameter, ambient supersaturation)` state trajectory,
- molar transfer chemical potential `Δμ(l−v)=RT ln(S_eq/S_inf)`,
- reversible-work landscape `∫ Δμ dN_w`,
- Maxwell–Mason condensational growth with vapor-diffusion and latent-heat resistances,
- controlled ambient histories,
- scrub-able particle history,
- spherical sixfold HexaLens and dense angular mode spectrum.

## Supercooled feedback laboratory

The activated droplet is continued below 0 °C as the same physical object rather than replaced by an unrelated ice mode.

### Conserved state

The evolved extensive variable is now **water mass**, with a conserved dry-solute inventory. Wet diameter is derived from:

1. dry-particle volume,
2. water mass,
3. pressure- and temperature-dependent liquid density,
4. Laplace pressure `Δp = 4 σ / D`.

This separates density-driven size changes from actual condensation/evaporation.

### Composition-derived water activity

`app/solver/composition.js` carries an explicit solute inventory:

- dry density,
- molar mass,
- effective van't Hoff factor,
- osmotic coefficient,
- selected water-activity model.

The default model is an explicit ideal effective-particle mole-fraction baseline:

`a_w = n_w / (n_w + ν φ n_s)`.

An exponential osmotic form is also available. `κ_effective` is now a diagnostic inferred from the current `(D,D_d,a_w)` state rather than the source of truth in this laboratory.

The default effective-solute constants are placeholders for experimentation, not a claim about aerosol chemical identity.

### Supercooled-water thermodynamic property layer

`app/solver/supercooled-eos.js` implements the IAPWS G12-15 Gibbs-energy formulation for cold/supercooled liquid water.

The module provides:

- density,
- isobaric and isochoric heat capacity,
- thermal expansivity,
- isothermal compressibility,
- sound speed,
- the IAPWS homogeneous-nucleation validity limit.

Density uses the analytic pressure derivative of the IAPWS Gibbs formulation. The remaining response properties are centered numerical derivatives of the same fundamental Gibbs/volume functions, avoiding disconnected empirical property fits.

The regression suite reproduces the published IAPWS verification states at 273.15 K / 0.101325 MPa, 235.15 K / 0.101325 MPa, 250 K / 200 MPa, 200 K / 400 MPa, and 250 K / 400 MPa.

### Laplace pressure + liquid/ice state

Liquid pressure is

`p_l = p_air + 4 σ_lv / D`.

The liquid IAPWS EOS is evaluated at this pressure. The low-pressure liquid/ice equilibrium activity is

`a_w^i(T) = p_i,sat(T) / p_w,sat(T)`.

A first **bulk equal-pressure** correction is included:

`ln a_w^i(T,p) = ln a_w^i(T,p_ref) + (V_i - V_l)(p-p_ref)/(R T)`.

`V_l` comes from IAPWS G12-15. `V_i` is currently an explicitly flagged constant-density Ice-Ih approximation (`ρ_i = 917 kg m⁻³`). The full IAPWS Ice Ih equation of state is the next refinement.

This is not yet the ice-embryo pressure correction. A finite embryo requires its own ice–liquid interfacial energy, geometry, curvature, and internal pressure; those belong to Milestone 2 nucleation physics.

### Surface tension boundary

The IAPWS R1-76(2014) surface-tension equation is used. Its stated supercooled extrapolation support extends only to about −25 °C, so the UI marks deeper-temperature use as an extrapolation rather than silently treating it as reference-quality.

### Water / ice saturation separation

The boundary laboratory keeps distinct:

- ambient `S_w`,
- ambient `S_i`,
- droplet surface equilibrium `S_eq`,
- composition water activity `a_w`,
- pressure-adjusted ice-equilibrium activity `a_w^i`,
- `Δa_w = a_w - a_w^i`,
- liquid-minus-ice chemical-potential drive.

These are different physical statements and are logged as separate events.

### Homogeneous-nucleation boundary

The Koop water-activity coordinate is exposed, together with exact Poisson conversion once a volumetric nucleation rate `J` is supplied:

`P_freeze = 1 - exp(-J V Δt)`.

No undocumented universal `J(T,a_w)` fit is inserted. Empirical homogeneous freezing, heterogeneous INP pathways, pore condensation/freezing, and ice-embryo geometry remain explicit Milestone-2 modules.

## Current model hierarchy

### Reference / reviewed formulation

- Murphy & Koop saturation vapor pressure over supercooled liquid water and ice.
- IAPWS R1-76 surface tension within its stated range.
- IAPWS G12-15 supercooled-water Gibbs equation of state.

### Explicit solution model

- conserved dry-solute and water inventories,
- effective van't Hoff factor and osmotic coefficient,
- ideal mole-fraction or exponential osmotic water-activity model.

### Transport approximation

- quasi-steady Maxwell–Mason spherical condensational growth,
- quasi-steady spherical external vapor field.

### Flagged approximation / deferred

- Ice-Ih molar volume currently from constant `ρ_i = 917 kg m⁻³`,
- surface tension below its stated −25 °C extrapolation range,
- no ice–liquid embryo curvature/pressure yet,
- no empirical homogeneous `J`, INP, or pore-freezing rate yet,
- no non-spherical vapor/heat field yet.

## Tests

Run:

```bash
node app/tests/thermo.test.mjs
node app/tests/activation.test.mjs
node app/tests/trajectory.test.mjs
node app/tests/supercool.test.mjs
node app/tests/supercooled-eos.test.mjs
node app/tests/composition-feedback.test.mjs
```

The new tests cover IAPWS published verification points, mass↔diameter inversion, density-driven size feedback at fixed water mass, condensation during cooling at fixed vapor pressure, composition-derived water activity, Laplace pressure, pressure-shifted liquid/ice equilibrium, and preservation of the exact nucleation-hazard algebra.

## Scientific references

- Petters, M. D. & Kreidenweis, S. M. (2007), *A single parameter representation of hygroscopic growth and cloud condensation nucleus activity*, Atmospheric Chemistry and Physics 7, 1961–1971. DOI: 10.5194/acp-7-1961-2007.
- Murphy, D. M. & Koop, T. (2005), *Review of the vapour pressures of ice and supercooled water for atmospheric applications*, QJRMS 131, 1539–1565. DOI: 10.1256/qj.04.94.
- Koop, T. et al. (2000), *Water activity as the determinant for homogeneous ice nucleation in aqueous solutions*, Nature 406, 611–614.
- IAPWS (2014), *Revised Release on Surface Tension of Ordinary Water Substance*, R1-76(2014).
- IAPWS (2015), *Guideline on Thermodynamic Properties of Supercooled Water*, G12-15.
- IAPWS (2009), *Revised Release on the Equation of State 2006 for H2O Ice Ih*, R10-06(2009) — identified as the next solid-phase property module.
