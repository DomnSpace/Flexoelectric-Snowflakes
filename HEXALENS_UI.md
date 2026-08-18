# HexaLens — sixfold physics visualization

## Core idea

Exploit the crystal's sixfold symmetry as a measurement coordinate system, not merely as decoration.

The primary graph is a **sixfold radial comparison instrument** centered on the simulated crystal. Each of the six sectors corresponds to one crystallographic arm direction `theta_i = i*pi/3`.

For any selected scalar field `q`, sample the same radial profile along every arm direction:

`q_i(r) = q(r cos(theta_i), r sin(theta_i))`.

Render all six profiles on a common scale in six 60-degree sectors around the crystal. In an exactly symmetric state the six sectors are visually congruent. Any symmetry breaking therefore becomes immediately visible as a deformation of the sixfold graph.

This is more informative than a conventional radar chart: radius keeps its physical meaning (distance from the crystal / interface), and angle keeps its crystallographic meaning.

## Default HexaLens channels

The user can swap the scalar encoded by the radial traces / filled ribbons:

1. water-vapor supersaturation `sigma(r,theta)`
2. vapor flux magnitude `|J_v|`
3. temperature excess `Delta T`
4. local interface velocity `v_n`
5. attachment coefficient `alpha`
6. interfacial chemical-potential shift `Delta mu`
7. mean curvature `H`
8. electric potential `phi` / field projection `E.n` (experimental layer)
9. polarization / flexoelectric contribution (experimental layer)

For interface-only quantities (`v_n`, `alpha`, `H`, local budget terms), the radial trace collapses to a highlighted interface sample plus a short history trail or angular ribbon rather than inventing values away from the surface.

## The six arm scopes

Each sector should carry a tiny aligned profile plot rather than a disconnected mini-chart. Suggested construction:

- inner radius: crystal/interface
- outer radius: far-field sample radius
- radial axis: physical distance
- transverse displacement or ribbon thickness: selected scalar magnitude
- common zero/reference circle or hexagon
- shared normalization across all six sectors

This makes one glance answer:

- are all six arms seeing the same vapor field?
- which arm is diffusion-starved?
- where is latent heating strongest?
- which facet has the largest attachment barrier?
- did the electric/flexoelectric layer actually break symmetry?

## Hexagonal contour field

Behind the arm scopes, use concentric hexagonal isolines as a physically meaningful display grid. The isolines can be either:

- coordinate guides only, or
- sampled contours of the selected scalar field.

When a continuum solver is active, do not force the physics itself into hexagonal symmetry. The display grid may be hexagonal, while the solved field remains free to break symmetry under flow, external fields, asymmetric branching, substrate effects, etc.

## Center crystal

The center remains the actual crystal / seed, rendered either:

- as the 3-D Three.js object projected into the HexaLens plane, or
- as a 2-D interface slice from the physics solver.

Clicking an arm or surface element selects the corresponding sector and opens the local budget inspector.

## Sixfold budget petals

Around the interface, a second optional six-petal encoding can show the **local growth budget** for each arm tip. Each petal is a signed stacked decomposition, not a generic radar score:

`driving = ambient supersaturation`

minus / modified by:

- diffusion depletion
- latent-heating penalty
- capillary / Gibbs-Thomson shift
- attachment barrier
- SDAK correction
- electric correction (experimental)
- flexoelectric correction (experimental)

The final petal radius is proportional to the resulting `v_n`.

Thus the physical snowflake and the explanatory snowflake become geometrically homologous.

## Preserve the other views

HexaLens is the **overview**, not a replacement for detailed plots.

### 1. Köhler / activation view

Full 2-D graph of equilibrium saturation vs wet diameter, with:

- hygroscopic / Raoult contribution
- Kelvin contribution
- combined curve
- ambient saturation line
- critical point
- animated particle state

Before freezing, one HexaLens sector can act as the compact summary of the current Köhler state, but the full curve remains available.

### 2. Phase / nucleation trajectory

Timeline / state-space view:

`dry -> hydrated -> activated -> supercooled -> nucleating -> ice seed -> faceted growth`

with mass, water activity, phase fraction and nucleation hazard.

### 3. Field slice

Full 2-D scalar/vector view for vapor, heat, flow, electric potential, etc. This is where the continuum solution is inspected without sixfold projection.

### 4. Crystal view

Full 3-D geometry with selectable physical surface coloring.

### 5. Local budget inspector

Click any surface point and show the numerical decomposition of its growth velocity / chemical-potential drive.

### 6. History / bifurcation view

Time series of symmetry-breaking metrics and morphology descriptors.

## Symmetry metrics

Make symmetry itself quantitative.

For a sampled arm quantity `q_i`, compute:

`q_bar = (1/6) sum_i q_i`

`A_6 = sqrt[(1/6) sum_i (q_i-q_bar)^2] / (|q_bar|+eps)`

and optionally angular Fourier modes

`Q_m = (1/6) sum_i q_i exp(-i m theta_i)`.

Useful interpretation:

- `m=0`: common six-arm mode
- `m=1`: dipolar asymmetry, e.g. imposed flow / field direction
- `m=2`: quadrupolar distortion
- `m=3`: alternating-arm mode

With only six principal arm samples the modes are diagnostic rather than a complete spatial spectrum; dense angular sampling can extend this later.

The UI can show these modes as tiny glyphs beside the HexaLens rather than another conventional chart.

## Interaction model

- tap/click a sector -> select arm `i`
- drag radial cursor -> inspect distance `r`
- scroll / pinch -> change sampled far-field radius
- toggle `lock scales` to compare all six sectors honestly
- toggle `difference mode` to plot `q_i(r)-mean_i q_i(r)`
- toggle `rotate with crystal` vs `laboratory frame`
- toggle `symmetrized baseline` to compare the current state with its sixfold average
- tap center -> return to global crystal selection
- long-press / secondary click -> pin a sector for A/B comparison

## Coordinate frames

Always show which frame is active:

- **crystal frame**: arm 0 stays at the reference direction; external flow / E rotates around it
- **laboratory frame**: external forcing remains fixed; crystal orientation rotates

This matters for separating intrinsic anisotropy from imposed anisotropy.

## Visual hierarchy

The visual order should be:

1. crystal morphology
2. selected physical field / radial sixfold traces
3. interface / tip budgets
4. forcing vectors and reference geometry
5. labels and numerical values

Do not let decorative snowflake geometry obscure physical quantities.

## First prototype

Before coupling to the real solver, build a data-driven SVG HexaLens with a mock state object shaped exactly like the future solver output:

```js
{
  r: [...],
  arms: [
    {theta:0, sigma:[...], flux:[...], dT:[...], vTip:..., alphaTip:..., Htip:...},
    ... six arms ...
  ],
  forcing: {flowAngle:..., EAngle:...},
  symmetry: {A6:..., modes:{m1:...,m2:...,m3:...}}
}
```

The mock values must be explicitly labeled synthetic. Once the Köhler and continuum kernels exist, the SVG should consume their actual arrays without changing its plotting API.
