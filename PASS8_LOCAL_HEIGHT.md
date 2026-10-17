# Pass 8 — Local ice height, prism/basal growth and explicit open mass ledgers

## Scope and truth boundary

Pass 8 adds an **optional experimental reduced-order local height** in the already tested triangular 2.5-D ice-growth kernel. Pass 7, the existing classic solid and the five optical sheets, and both original numerical solvers remain intact when this experiment is off.

**What is actually solved:** Every occupied hexagonal lattice cell now has independent `topHeight[i]` and `bottomHeight[i]` in SI meters. Newly occupied prism cells inherit the weighted physical height of their adjacent ice facets. Basal attachment grows upper/lower surfaces by mass admitted through an explicitly **unresolved external vertical vapor reservoir**. The local ice volume and mass are computed from those actual heights and ice phase fractions; both contributions source the heat calculation.

**What is NOT solved:** The vapor concentration and temperature are still quasi-steady fields in **two dimensions**. There is no z-direction diffusion PDE, 3-D vapor flux, real vertically resolved temperature, ventilation or closed ice + vapor finite-volume storage balance. Claiming "fully mass-conserving snowflake physics" would therefore be incorrect. Ice-side mass closure is strictly against *separately recorded external source terms*.

The 3-D viewer still uses a global mean-depth extrusion or five **mean-depth** transparent sheets; **neither shows the newly computed local top/bottom shape yet**. Pass 9 will extract actual local z-contours from these fields. For Pass 8, select **Spatial Fields → Local ice thickness h [µm]** to inspect genuine computed height variation.

## Fields and geometry

At each physical cell of area `A_hex = (sqrt(3)/2) dx²`:

```
ice volume_i = A_hex * phi_i * max(0, topHeight_i - bottomHeight_i)
ice mass     = rhoIce * sum_i(ice volume_i)
```

The seed initializes `topHeight = +h0`, `bottomHeight = −h0` only for actually occupied ice cells, with `h0=basalHalfThickness`; gas cells have zero height. `phi` is the original per-cell ice fraction. No mass is synthesized by the new renderer, and no initial seed mass is counted as post-start deposition.

### Prism (lateral) mass

The existing 6-neighbor interface geometry, anisotropic attachment and Gibbs–Thomson vapor equilibrium remain the basis. Local prism face area is the summed occupied-neighbor **physical depth × hex-face length**, rather than a single imposed global slab depth. The admitted external prism mass is

```
DeltaM_prism = min(positive_flux * local_prism_face_area * dt,
                   rhoIce * A_hex * inherited_depth * (1 - phi_old))
```

The ice fraction increases by the admitted mass divided by the local available volume capacity, not by an unbounded growth multiplier. The local upper/lower elevations are inherited from the neighboring actual solid. A cell full of ice cannot accept more prism mass through the same face.

### Basal (vertical) mass

Upper and lower basal surfaces are fed by the **unresolved external vertical reservoir**, with prescribed `c_vertical = c_sat,ice(Tinf) * (1 + sigmaInf)`; local equilibrium `c_eq = c_sat,ice(Tinf + DeltaT[i])` and temperature-dependent basal attachment coefficient `alpha_basal(sigma_surface)`. A reduced-order diffusion/kinetic resistance computes the per-face vapor mass flux:

```
j_z = max(0,c_vertical - c_eq) /
      [ max(dx,local_thickness/2)/D_v + c_eq/(rhoIce * alpha_basal * v_kin) ]
```

The zero-attachment limit is identically zero; timestep depth increments are bounded by 0.1 dx per face as a numerical guard. For each occupied footprint,

```
delta h_top  = + j_z dt / rhoIce
delta h_bot  = - j_z dt / rhoIce
DeltaM_basal = rhoIce * A_hex * phi_i * (delta h_top - delta h_bot)
```

This closure does **not** deplete stored vapor above/below the plate. Its cumulative input appears as `basalVerticalInKg` in an explicitly **open** ledger. In a future physically closed coupled vapor/ice model, the corresponding gas control volume would lose exactly this mass.

### Heat and bookkeeping

New admitted prism and basal masses source the next local thermal relaxation as latent heating, each **once**. Cumulative `totalDepositedMass` and `totalLatentHeat` reflect both growth channels in the local-height mode (unlike the old half-height visual-only growth). The exact ice-side check is:

```
rhoIce * A_hex * sum_i(phi_i*(top_i-bottom_i))
    = M_initial + M_from_open_planar_reservoir + M_from_open_z_reservoir
      + epsilon_numerical
```

The solver records **residual mass in kilograms and relative residual** every step and rejects significant drift. This is a geometrical ice-side consistency condition, **not** a claim that the open gas boundary system satisfies closed-system conservation.

## Compatibility & domain expansion

- **Legacy mode (default):** Existing classic triangular and Cartesian kernels remain unchanged; snapshot schema v1, original `basalHalfThickness` scalar, existing electrical Cartesian capabilities.
- **Pass 8 enabled:** Hexagonal mode only; snapshot schema v2 contains copied `hTop[]`, `hBottom[]`, model flags, reservoir ledger and measured geometry. Legacy v1 replay is untouched.
- **Automatic 73→...→401 grid expansion:** Copies the entire top/bottom fields, phase, gas, temperature, heat source, accumulated ice mass and the source ledger to the enlarged hexagon; domain enlargement itself must not manufacture ice.
- **Render:** Solid uses a mean-height proxy. The five transparent planes use a mean-height guide, **not real per-cell z-level contours**. The thickness heatmap displays the actual local field in micrometers.

## Acceptance tests

```sh
node app/tests/pass8-thickness-budget.test.mjs
node app/tests/run-all.mjs
```

The Pass 8 suite must verify: initialized seed mass; prism mass capacity and actual occupied volume; equal top/bottom basal deposition; latent heat equals admitted mass times sublimation enthalpy; identically zero mass growth at zero kinetic accommodation; real time-evolving local height and positive prism/basal rates; no increase of ice mass without accounted external inputs; strict ice-mass ledger drift; 73→97 continuation without discontinuity; schema-v2 copied typed arrays; and correct honesty labels for unchanged mean-height optical five sheets.

The Chromium test additionally verifies opt-in v2 behavior, new local thickness spatial field, available CSV channels, legacy v1 return and the original electric solver as an independent control.

## Still missing (Pass 9 and beyond)

- Extract actual multiple z-contours (including holes) from spatial `topHeight`/`bottomHeight` + `phi`, using marching triangles, and illuminate real top/bottom gradients in five translucent physical slices. No photo-real height structures before that.
- Vertically resolved vapor supply and energy transport / reduced-order calibration against measured basal growth.
- Calibration of the microscopic attachment coefficient and real terrace-width SDAK.
- Fully closed finite-volume vapor + ice accounting (rather than current external-reservoir ice-side conservation), surface diffusion and stochastic fluctuation consistent with mass transport.
