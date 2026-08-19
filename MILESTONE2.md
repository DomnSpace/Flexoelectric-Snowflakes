# Milestone 2 — Nucleation, Critical Embryo, and Ice Seed

Milestone 2 starts from the same supercooled, composition-resolved droplet state built in Milestone 1. It does **not** jump directly to a six-armed crystal.

## 1. Executable homogeneous pathway

The stochastic event model uses the Koop water-activity fit

`log10 J_hom[cm^-3 s^-1] = -906.7 + 8502 Δaw - 26924 Δaw^2 + 29180 Δaw^3`

with `Δaw = aw - aw_i(T)` and fit window `0.26 < Δaw < 0.34`. The solver refuses silent extrapolation.

`lambda = J V_liquid`

`H(t) = integral lambda dt`

A seeded exponential threshold `H* = -ln(U)` is sampled once. Nucleation occurs when `H >= H*`. Rendering never advances hazard.

## 2. Full Ice-Ih thermodynamics

`app/solver/ice-ih-eos.js` implements IAPWS R10-06(2009), the Gibbs-energy EOS for ordinary hexagonal Ice Ih.

It provides `g, h, s, rho, cp, alpha, beta, kappa_T, kappa_S` from the same thermodynamic potential used for liquid-solid equilibrium.

The supercooled-liquid EOS now also exposes Gibbs-derived `h` and `s`, so liquid and solid phase energetics share compatible IAPWS reference conventions.

## 3. Direct liquid–ice driving force

For solution water,

`g_l,solution = g_l,pure + R_v T ln(aw)`

and

`Delta mu_mass = g_l,solution - g_ice`.

Positive `Delta mu_mass` means Ice Ih is the lower-Gibbs state. A CNT bulk free-energy-density reference is

`Delta g_v = rho_ice Delta mu_mass`.

This replaces the old constant-density pressure correction as the primary embryo-driving diagnostic.

## 4. CNT embryo reference — deliberately separate from Koop J

Classical spherical CNT is used only as an interpretable geometry/barrier reference:

`r* = 2 gamma_il / Delta g_v`

`Delta G* = 16 pi gamma_il^3 / (3 Delta g_v^2)`.

The default `gamma_il = 0.029 J m^-2` is a source-tagged reference parameter motivated by TIP4P-family seeding/CNT literature and is user-adjustable. It is not claimed universal.

The empirical Koop rate and the CNT barrier are displayed side by side. The code does not tune one to force agreement with the other.

## 5. Post-critical seed birth

The arbitrary initial frozen-mass fraction has been removed.

After the Koop hazard crosses its stochastic threshold, the first seed is sized from a slightly post-critical CNT radius (`1.05 r*` by default), capped only by available water mass.

The resulting seed contains:

- total H2O mass,
- ice and liquid masses,
- exact mass residual,
- critical radius and critical molecule count,
- actual post-critical radius,
- Ice-Ih EOS density,
- remaining liquid composition,
- nucleation pathway and timestamp.

No arms or dendrites are created at seed birth.

## 6. EOS-derived latent heat

The previous fixed 333.55 kJ/kg constant has been removed from seed birth.

At the event state,

`L_f(T,p) = h_liquid(T,p) - h_ice(T,p)`

is evaluated from the two EOSs and

`Q_latent = m_ice L_f`.

The adiabatic temperature jump remains an upper-bound diagnostic until the released heat is coupled into the spatial temperature solver.

## 7. Other nucleation pathways

### Immersion

The code exposes

`log10 J_het[cm^-2 s^-1] = m Δaw + c`

but `m,c` are material-specific. No universal default is supplied.

### PCF

Pore-condensation freezing remains a real but non-executable pathway pending pore radius/geometry distributions, wetting/contact angle, condensation, pore-water freezing, and ice-escape conditions.

## 8. Visualization

- `app/milestone2.html`: stochastic hazard, empirical Koop rate, conserved phase inventory, seed birth.
- `app/embryo.html`: Ice-Ih/liquid Gibbs comparison, CNT `Delta G(r)`, `r*`, `Delta G*`, critical molecule count, post-critical seed mass, EOS latent heat.

The sixfold frame remains crystallographic orientation only.

## 9. Still deferred before faceted crystal growth

- spatial redistribution of latent heat,
- ice-liquid interfacial-energy anisotropy rather than scalar `gamma_il`,
- non-spherical embryo geometry / stacking disorder,
- material-calibrated immersion freezing,
- executable PCF pore model,
- coupled non-spherical vapor + heat solver,
- basal/prism attachment kinetics and SDAK,
- 3-D faceted interface growth.

## 10. Tests

- `app/tests/nucleation.test.mjs`
- `app/tests/ice-ih.test.mjs`
- `app/tests/embryo.test.mjs`
- plus the Milestone-1 thermodynamics/composition/supercooling tests.

## Primary references

- Koop, T., Luo, B., Tsias, A. & Peter, T. (2000), *Water activity as the determinant for homogeneous ice nucleation in aqueous solutions*, Nature 406, 611–614.
- IAPWS R10-06(2009), *Revised Release on the Equation of State 2006 for H2O Ice Ih*.
- IAPWS G12-15, *Guideline on Thermodynamic Properties of Supercooled Water*.
- Espinosa et al. (2013), *Homogeneous Ice Nucleation at Moderate Supercooling from Molecular Simulation*, JACS 135, 6108–6111.
- Knopf & Alpert (2013), activity-based immersion freezing framework.
