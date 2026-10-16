# Model B — Five independently evolving ice slabs (Step 3/15)

**State:** experimental physics branch; measured runtime and mass conservation are regression-tested. The main 300 µm legacy Observatory is unchanged except for a link to an **independent small-grid laboratory** at `app/layered-lab.html`.

## Which physical mechanism is now implemented?

Five separate ice-eligible z slabs are embedded between two upper and two lower gas-only slabs. The computational geometry is a coarse extrusion of the 60° triangular lattice, using hex Voronoi areas. Real water-vapor mass concentration `c[k]` [kg/m³] is transported in a **single shared connected gas domain** with top/bottom/lateral open or sealed boundaries.

`layered-vapor-implicit.js` solves backward-Euler diffusion and attachment using a symmetric positive-definite matrix and Jacobi-preconditioned conjugate gradients. The matrix-vector multiplication is expressed in flux differences to reduce floating-point cancellation at large implicit step ratios. The linear solve is checked against its **true** residual and projected onto the exactly conserved common mass mode; this correction does not waive the full residual requirement. A failed solve throws *before* changing water masses, ice, time or deposition ledgers. With cell spacing 0.75 µm, full-scale, high-resolution convergence is **not yet demonstrated**.

For an active gas cell i with free volume fraction fᵢ, and an exposed ice surface, the equation has the reduced finite-volume form:

    f_i V (c_i^{n+1}-c_i^n) / dt
      = Σ_gas_neighbors K_ij(c_j^{n+1}-c_i^{n+1})
        + Σ_open_faces K_boundary(c_boundary-c_i^{n+1})
        - Σ_exposed_ice_faces k_attachment A_face (c_i^{n+1}-c_eq).

The condensation-only model requires `c >= c_eq`. It does not handle sublimation or adsorption hysteresis, and does not yet solve heat. The vapor/thermal references remain at the fixed ambient temperature. Each gas-to-ice face event records its physical mass and the ice slab it deposits into. Mass is subtracted exactly once from the common gas inventory.

## Five independent shape fields

`layered-ice-front.js` stores `iceFraction[k]` and `iceMassAtCellKg[k]` separately in all five ice-eligible slabs. These represent *coarse voxel volume occupation*, not individual molecular terraces or optically prescribed contours.

Every such slab can be initialized with the same small axial-radius seed. No layer has a preferred kinetic multiplier: upper/lower planes start symmetric. Growth proceeds from exposed gas/ice interfaces. Once a gas voxel acquires ice, its available vapor volume shrinks accordingly; its partial-ice content provides attachment sinks to exposed neighbors in subsequent steps. Where a voxel becomes completely full, it converts into solid; previously contained vapor is displaced to adjacent accessible gas volumes. If a voxel has no gas neighbors, the remnant is recorded in a separate `trappedVaporKg` inventory, **not destroyed**. Deposits on a fully occupied upper/lower basal voxel (facing a gas-only halo) are booked as added unresolved local ice thickness, not new ice permitted inside gas-only halo.

The method is a **first conservative front-motion closure**, not yet a calibrated phase-field interface. The volume fraction changes continuously, but individual cell faces remain coarse. Any excess mass beyond a completely filled voxel stays in its physical mass ledger rather than spawning fake adjacent cells.

The mass closure includes all gas inventories, accumulated newly deposited ice and any trapped vapor:

    M_water(t) = Σ_gas c_i V_i(free fraction)
                 + Σ_ice deposited mass + M_trappedVapor.

In an open domain, the change equals recorded reservoir inflows minus outflows; in sealed conditions it remains constant within tested tolerances. Independently, `Σ physicalIceMass - initialSeedIceMass = accumulatedNewIceMass`.

**Important:** the thermal solver is NOT implemented. The latent heat ledger contains deposited mass × Ls, but no temperature diffusion or feedback from this deposition occurs. So the five current patterns cannot yet be presented as calibrated snow-crystal ridges or genuine light-refracting thickness profiles.

## Separate browser experiment

The new [Model B laboratory](app/layered-lab.html) shows five **actually distinct numerical ice-fraction maps**, model time, total deposited mass, true CG residual, water-inventory residual and separately measured deposited mass by z slab.

Use the Step / 10 steps controls. The supplied contrast selector applies a *display-only* gain to very small subcell ice fractions (1×, 100×, 1000×); it is not added ice or a model parameter. Starting seed is the same in all five slabs; the outer two have directly exposed basal surfaces, while the central slab is more shielded. Depending on rates and gas exposure, the outer layers can grow faster. This refutes any assumption that a middle spine naturally has to dominate.

The legacy viewer's `ICE VIEW: 5 translucent height slices` remains a preview of one 2.5-D slab. It is **not** silently fed synthetic model B data. The experimental lab is separate until we have a defensible thickness-to-geometry mapping and the mass/heat coupling mature.

## Tests and limitations

Run:

    node app/tests/layered-implicit.test.mjs
    node app/tests/five-layer-ice-front.test.mjs
    node app/tests/run-all.mjs

Validation covers nonuniform 3D gas diffusion, positive attachment, shared resource competition, symmetric upper/lower growth, buried non-exposed ice, sealed and open boundaries, strict water ledgers, atomically refused nonconverged steps, physically separate ice fraction maps, non-negative gas inventories and conservative promotion of a nearly full ice voxel. Chromium smoke checks actual five-panel rendering, one numerical step, different exposed/central layer masses and that display gain does not change mass.

**Remaining before extending to the full 300 µm domain:**

- Full 3D temperature diffusion + latent-heat backreaction, including independent top and bottom ice surfaces.
- Consistently fitted basal/prism attachment kinetics and curvature/SDAK on moving 3D volume-fraction faces.
- A converged implicit solver at N=401 with 9 slabs, residual verification and memory benchmarks (PCG is currently only tested on small grids).
- Better interfacial geometric reconstruction and partial-cell surface areas than the present coarse occupancy closure.
- Physically meaningful vertical grid spacing/halo convergence against microscopic snow-crystal observations; no preference for a central layer is imposed.
- Integration of the independent five physical fields into the legacy 3D five-slice optical renderer and snapshot schema v2, after conservation and temperature budgets have been validated.

The existing `hex25d.js`, `growth25d.js`, `ice-snapshot-v1.js`, `uniform-five-layers.js` and original 300 µm viewer physics remain unchanged. Frozen baseline numerical fingerprints continue passing as a separate regression.
