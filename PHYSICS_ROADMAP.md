# Flexoelectric Snowflakes — Physics-first v2 roadmap

This project is no longer organized around a decorative growth law. The target is a continuous microphysical trajectory from aerosol hydration through droplet activation, freezing / ice nucleation, and then anisotropic diffusion-limited crystal growth. Electric and flexoelectric effects are added only as explicitly marked extensions on top of the thermodynamic and transport baseline.

## 0. State variables

A particle/crystal state should carry enough information to survive the phase transition without replacing the model by an unrelated one:

- ambient: `T`, `p`, water-vapor partial pressure `e`, `S_w=e/e_sw(T)`, `S_i=e/e_si(T)`
- aerosol: dry diameter `D_d`, hygroscopicity `kappa_hyg`, soluble volume / optional explicit water activity model
- condensed water: wet diameter `D_w`, liquid water activity `a_w`, surface tension `sigma_lv`
- nucleation: phase fractions / state label, pore geometry where enabled, nucleation hazard / probability
- ice: interface geometry, principal curvatures, facet orientation, local surface supersaturation `sigma_surf`, surface attachment coefficient `alpha`
- transport fields: vapor concentration / supersaturation and temperature
- optional electro-mechanical extension: electric potential/field, polarization, strain and strain gradient

Use `kappa_hyg` for Petters–Kreidenweis hygroscopicity. Never reuse `kappa` for geometric curvature; use `H` (mean curvature) or `K_geom` instead.

## 1. Liquid activation: κ-Köhler

For wet diameter `D` and dry diameter `D_d`, use the Petters–Kreidenweis form

`S_w,eq(D) = [(D^3-D_d^3)/(D^3-D_d^3(1-kappa_hyg))] * exp(A/D)`

with

`A = 4 sigma_lv M_w / (R T rho_w)`.

The activation point is the local maximum of `S_w,eq(D)`. The simulation should solve it numerically, not hard-code the usual large-particle approximation. The UI should show:

- Raoult / hygroscopic contribution
- Kelvin curvature contribution
- full Köhler curve
- ambient `S_w` horizontal line
- critical wet diameter and critical supersaturation
- the actual particle state moving along the curve

This stage should make deliquescence/activation visible as a stability problem, not a binary threshold.

## 2. Beyond spherical Köhler: generalized interfacial chemical potential

The spherical Kelvin term is only a special case of curvature-controlled interfacial chemical potential. For a general interface use principal curvatures `1/r1 + 1/r2 = 2H` and write the capillary contribution schematically as

`Delta_mu_cap = gamma * v_m * (1/r1 + 1/r2)`.

The renderer should therefore treat the familiar Köhler curve as a 1-D slice through a more general free-energy / chemical-potential landscape. This is the bridge we need later when an ice surface is faceted, curved, branched, or confined in a pore.

## 3. Ice nucleation bridge

Do not jump directly from an activated droplet to a six-armed crystal. Provide selectable nucleation pathways:

1. homogeneous freezing of an aqueous droplet (water-activity based parameterization / hazard model),
2. immersion freezing on an ice-nucleating particle,
3. pore condensation and freezing (PCF) for porous solid aerosol,
4. manual `seed ice` mode for isolating post-nucleation growth physics.

For PCF, capillary condensation is computed from Kelvin / Young–Laplace geometry, then freezing occurs inside the pore, and macroscopic ice growth is only allowed when the ice phase can escape the pore under the ambient ice supersaturation.

The visualization should show the phase event explicitly: aerosol core -> hydrated droplet / porewater -> critical ice embryo -> faceted seed.

## 4. Vapor + heat transport around ice

Once ice exists, evolve fields rather than assigning growth independently to each arm.

Baseline continuum problem:

`partial_t c = D_v nabla^2 c - u . grad(c)`

`partial_t T = alpha_th nabla^2 T`

with latent-heat release coupled at the moving ice interface.

At the ice surface, mass flux and interface velocity are coupled to local supersaturation and attachment kinetics. The first performant browser implementation can use a 2-D/2.5-D hexagonal field or sparse voxel field; the 3-D mesh is then a visualization of the interface state rather than the source of truth.

Important: diffusion focusing at corners and tips must emerge from the field solution. Branching should not be triggered by an arbitrary curvature threshold once the field model exists.

## 5. Facet attachment kinetics

Use distinct basal and prism attachment coefficients. A useful measured/empirical form is

`alpha(sigma_surf,T) = A(T) * exp[-sigma_0(T)/sigma_surf]`.

Add structure-dependent attachment kinetics (SDAK) as an optional, clearly identified semi-empirical layer in which the nucleation barrier depends on local facet width / mesostructure. This is the route to plates, columns, hollow columns and dendrites without a hand-authored Nakaya lookup table.

## 6. Surface curvature, anisotropy and morphology stability

The local equilibrium vapor pressure over ice should include anisotropic capillarity / Gibbs–Thomson effects. Numerically this can be represented using an orientation-dependent surface energy `gamma(n)` and its stiffness, or with a phase-field / front-tracking approximation.

The system should expose three distinct causes of morphology:

- transport instability: diffusion / heat focusing,
- capillary regularization: curvature + anisotropic surface energy,
- attachment kinetics: basal/prism/SDAK response.

The visualization must let the user toggle these contributions independently.

## 7. Electric and flexoelectric extension — hypothesis layer

Only after the neutral baseline reproduces reasonable morphologies should electro-mechanical terms be enabled.

Track:

- electric potential `phi`, `E=-grad(phi)` with dielectric boundaries,
- polarization `P`, initially with a simple effective constitutive model,
- strain / curvature / strain-gradient observables from the ice interface.

Candidate extension:

`P_i = mu_ijkl * partial_l epsilon_jk`

and a free-energy contribution coupling polarization to the electric field, schematically

`f_E = -P . E + 1/2 epsilon |E|^2`.

Any resulting change to interfacial chemical potential / attachment kinetics must be reported separately from the baseline terms. Do not bury an arbitrary multiplicative `1 + alpha_E E` inside the growth law.

## 8. Visualization: physics becomes geometry

The app should have synchronized views:

### A. Microphysical trajectory
A phase-axis showing `dry aerosol -> hydrated aerosol -> activated droplet -> nucleation event -> ice seed -> growing crystal`.

### B. Köhler / free-energy view
Interactive curve(s) with the particle state moving through radius–supersaturation space. Once frozen, this view transforms into local interfacial chemical-potential / supersaturation diagnostics rather than disappearing.

### C. Field view
Slice through vapor supersaturation, temperature, electric potential (optional), and flux vectors around the crystal.

### D. Crystal view
3-D faceted ice geometry. Surface coloring must map to a selected physical scalar: `sigma_surf`, mass flux, `alpha`, mean curvature, temperature excess, polarization, or growth velocity.

### E. Budget inspector
At the selected surface element show contributions to the local driving force / velocity: ambient supersaturation, capillary shift, diffusion limitation, latent-heat limitation, attachment barrier, SDAK correction, electric/flexoelectric correction.

This inspector is the core anti-handwaving feature: every visible growth event should be decomposable into the modeled physics.

## 9. Implementation order

### Milestone A — Köhler laboratory
- robust saturation-vapor-pressure functions over liquid water and ice
- κ-Köhler curve and numerical critical point
- time integration of wet radius under prescribed ambient history
- interactive curve + particle marker

### Milestone B — phase bridge
- homogeneous / immersion / PCF pathway interfaces
- explicit phase transition state machine
- ice seed generation with mass conservation

### Milestone C — neutral ice growth kernel
- vapor diffusion field
- heat diffusion + latent heat
- anisotropic basal/prism attachment kinetics
- capillarity
- 2-D/2.5-D moving interface

### Milestone D — morphology
- SDAK / facet-width dependence
- field-generated branching and tip splitting
- map neutral kernel to the 3-D visualization

### Milestone E — electro/flexo hypothesis laboratory
- electric potential field
- polarization / strain-gradient observables
- explicit free-energy/chemical-potential coupling experiments
- side-by-side neutral vs electro/flexo runs

## 10. Scientific labeling

Every parameter and term in the UI gets one of:

- `thermodynamic` — established equilibrium relation
- `transport` — established continuum relation / material property
- `measured-fit` — fitted attachment/nucleation parameterization
- `semi-empirical` — morphology model such as SDAK
- `hypothesis` — proposed electric/flexoelectric coupling
- `visual-only` — rendering amplification with no physical feedback

No quantity may silently migrate between these categories.

## 11. Immediate refactor of the current prototype

The current prototype overloads `kappa` as curvature and uses a direct multiplicative electric growth bias. In v2:

- rename geometric curvature to `H` / `K_geom`,
- reserve `kappa_hyg` for aerosol hygroscopicity,
- remove the hand-authored electric multiplier from the neutral baseline,
- treat the current arm wiggle / flexoelectric system as an optional experimental module,
- make the continuum thermodynamic/transport kernel the source of growth velocity,
- make Three.js consume field/interface outputs rather than invent growth itself.

The old prototype remains useful as a rendering and interaction shell, but not as the final physics engine.
