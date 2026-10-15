# Step 2/15 — conservative gas finite-volume pilot (no five-layer crystal growth yet)

**Isolation:** `app/solver/layered-gas-fv.js` is an independent opt-in calculation with separate tests. It never imports, alters or replaces `hex25d.js`, the Cartesian growth solver, the Pass 7 snapshot schema or the five-layer *optical preview*. The existing Snowflake Observatory's Pages/3D viewer remains unchanged.

## Actual implementation

An equilateral triangular lattice on the basal plane, axial domain of odd width N=7..97, is extruded into `iceSlabs=5` ice-eligible z layers surrounded by `gasHaloSlabs=2` gas-only layers on **each side**. The default is 9 total layers. These are **candidate computational parameters**, not measured atmospheric z resolution.

Actual positions use (x=a(q+r/2), y=\sqrt{3}ar/2); each basal Voronoi area is (A_h=\sqrt{3}a^2/2). For gas slab thickness `dz`:

- Gas CV volume `V=A_h dz`, gas mass `m=cV`; `c` is absolute vapor mass density in kg/m³, not supersaturation.
- Lateral gas-gas face area `a dz / √3`, center spacing `a`.
- Vertical gas-gas face area `A_h`, center spacing `dz`.
- Signed gas-gas transfer: `dm = D * A_face / distance * (cA-cB) * dt`, identical debit and credit, calculated once per undirected face.
- Ice is a fixed `solid[k]` cell mask **only in five eligible slabs**. Exposed solid-gas faces receive flux requests `dm = k_attach * A_face * max(cGas - cEq,0)*dt`. Lateral faces use `kPrism`, vertical faces use `kBasal`. No deposition takes place on solid-solid interfaces, including fully buried ice. **Rates are placeholders** and not calibrated facet physics.
- The solid mask does **not** advance or change volume with deposited mass. Deposits accumulate in a separate ice mass ledger, so conservation is tested before moving-interface physics.
- Open axial/radial reservoir faces, if requested, use Dirichlet values at a half-cell distance. Boundary influx and outflux are explicit mass ledgers. Default is completely sealed.
- All possible gas-gas, ice-deposition and boundary outflows are evaluated against the **same old gas concentration**, then summed in one outgoing demand per gas cell. A cell cannot spend the same vapor twice. Default steps reject an excessive `dt` or overdraw **without modifying any state**. An explicit `allowFluxLimiter:true` experimental option distributes the available inventory proportionally across the competing demands, preserves positivity and records how much the requests were curtailed.
- Every kg deposited contributes `Ls` J to a **released latent heat ledger**. **No thermal PDE or heat removal is solved yet**, clearly marked `heatTransportSolved:false`.

The model contains the entire hexagonal 2-D transport geometry and coarse z slabs (a low-resolution **3-D numerical transport**). It does not simulate molecule-scale water, fully faceted moving interfaces or independent crystal shapes. Its slab count and dz need convergence tests before using large domains.

## Explicit time-step stiffness

At 0.75 µm lateral spacing, real vapor diffusivity implies a micro- to nanosecond-scale stability restriction for a naive explicit update. The prototype reports `layeredExplicitDtLimit(state)` and refuses larger `dt` by default. **A regular Observatory timestep of 0.5 s cannot be inserted here.** A future implicit diffusion solve or converged multigrid strategy is essential before browser-real-time coupled growth at N=401. The optional limiter is primarily a mass-accounting stress test and is not a substitute for temporal convergence.

The Step-2 cap `N<=97` is intentional to prevent accidental huge temporary candidate-flux arrays or misleading claims about runtime on the full 300 µm domain.

## Conservation contract

Define tracked inventory `M = Σ_gas c_i V_i + Σ_ice deposit_i`. The pre-existing (fixed) solid ice seed mass does not vary and cancels from the audit. For a closed boundary:

    M(t+dt) - M(t) = numerical roundoff.

For open boundaries:

    M(t+dt) - M(t) = boundaryInKg - boundaryOutKg + rounding.

Deposition also increments `latentReleasedJ += depositedKg * Ls`, without claiming the latent heat was thermally redistributed.

Test `node app/tests/layered-gas-fv.test.mjs` demonstrates:

1. Actual triangular-domain geometry and nine slabs with five ice-eligible levels.
2. Closed, nonuniform lateral/vertical diffusion and invariant total gas mass across many small steps.
3. A real, single fixed ice surface captures vapor from the same gas inventory, with opposite gas loss.
4. A **fully buried core** receives zero deposition, while exposed neighbors grow.
5. Top/bottom symmetric ice exposure produces symmetric vapor depletion.
6. Two distinct ice layers compete for the same gas cell, with no duplicate inventory spending and per-ice-slab flux ledgers.
7. Artificial extreme adhesion is rejected when `dt` violates stability; optional positive limiter records curtailed demands.
8. Open reservoir inflow/outflow accounts for all exchanged mass.
9. A change of `dz` changes gas volume and vertical stiffness.
10. The original frozen hex/Cartesian baselines remain numerically identical, covered separately by `five-layer-baseline.test.mjs`.

## Strict boundary of this milestone

- No moving interface, reconstructed five physical snowflake contours or middle-spine preferential rate yet.
- No local heat diffusion, ambient parcel dynamics or vertical latent-heat backreaction.
- No layer-specific SDAK, basal/prism temperature-dependent calibration, terrace width, stochastic deposition or true ice height field.
- The old viewer still shows five **identical transparent optical layers**. That is intentionally unchanged. They will become independently shaped only after this finite-volume kernel proves mass conservation and gains physically valid ice motion.

## Next: Step 3/15

Introduce five independent fixed-geometry ice-occupancy slabs + thermal proxy with **one shared FV gas solver** (not five separate 2-D vapor copies). Start with measured surface-specific deposition, then derive a conservative shape/volume update and a v2 snapshot. Stochastic, edge/facet effects remain excluded until the mass/heat exchanges are consistent and tested. Examine closed energy-budget controls before claiming thermodynamic conservation, which this Step 2 does not solve.
