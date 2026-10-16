# Step 4/15 — Model B thermal coupling and stellar 30° startup

**Scope:** default the *existing* successful stellar configuration in the main Observatory, and add an optional **conservative gas/ice thermal solver** to the independent Model B test laboratory. The original `hex25d.js`, Cartesian electric `growth25d.js`, five optical preview layers, and frozen numerical fingerprints are unchanged.

## 1. Why the triangular 30° stellar experiment is now the startup configuration

The user-supplied screenshot from 2026-10-16 showed six strong primary growth axes and interesting secondary branches, broadly closer to photographed stellar dendrites than the earlier split-tip cases. Previous controlled phase tests found an improved center-minus-flanks occupancy metric for a 30° kinetic phase under **one early-time domain and uncalibrated kinetic set**, but did not demonstrate mesh independence or universal correctness.

The start configuration now comes from one named `GROWTH_HABITS.dendrite` record, not two partly inconsistent hardcoded defaults:

| Parameter | Stellar default | Notes |
|---|---:|---|
| Topology | **Triangular (6-neighbor)** | Hex Voronoi transport controls, not hydrogen bonds |
| Kinetic phase | **30°** | Relative to local outward interface normal |
| Temperature | −15 °C | Stellar laboratory-inspired condition |
| Ice supersaturation | 14 % | Not fitted to observations |
| Sixfold kinetic anisotropy | 0.50 | Experimental parameter |
| Prism nucleation barrier | 0.0035 | Uncalibrated |
| Basal attachment prefactor | 0.004 | Uncalibrated |
| Seed | Round, 1.5 µm | No imposed six-arm skeleton |
| Surface-energy epsilon | 0 | Control remains regular; optional 0.01/0.02 |
| SDAK | On | Simplified existing local enhancement |
| Capillarity | On | Existing Gibbs–Thomson |
| Grid | 129 (96 µm) | Same a=0.75 µm |
| Hex domain expansion | On, up to 401 / 300 µm | Finite boundary validity checks remain active |
| Electric coupling | **Off** | Electric Poisson solve is Cartesian-only |

Loading `Stellar candidate` must also restore 30°, **never silently revert to 0°**. Plate/sector scenarios retain their own 0° control by default. All controls can be changed explicitly. The original square/electric branch remains selectable but is not misrepresented as triangular electrostatics.

## 2. Gas/ice heat equation

Model B keeps five separately occupied ice-eligible z-slabs and a shared nine-slab (5 ice eligible + 2 gas halo each side) vapor volume. It adds a separate `app/solver/layered-thermal.js` implicit finite-volume thermal update:

    C_i (T_i^{n+1}-T_i^n) =
       Δt Σ_{j exposed cell faces} K_ij (T_j^{n+1}-T_i^{n+1})
       + Q_latent,i
       + Q_ambientBoundary,i.

Here `T_i` is a local perturbation from fixed far-field `Tinf`, `Q_latent = ΔM_deposited × Ls` is in joules and **is released**, not consumed by diffusion. `K_ij = A_face/d * harmonic(k_i,k_j)` [W/K] and cell heat capacity is an explicit mixed-voxel approximation:

    C_i = V [ (1-φ_i) ρ_air cp_air + φ_i ρ_ice cp_ice ].

`ρ_ice, cp_ice` use the existing IAPWS Ice Ih module at the far-field state. Air density uses an ideal-dry-air approximation and cp_air = 1005 J/kg/K. Ice thermal conductivity is set to **2.2 W/(m·K), an illustrative fixed coefficient**, not facet-calibrated temperature anisotropy. The scalar conductivity `k_i` is the volume-weighted ice/gas blend; the shared face uses the harmonic mean. No convection, radiative transport, complex crystal-interface resistance, or detailed z microstructure is solved.

The primary stored thermal variable is *cell excess enthalpy in joules*. If the crystal occupies more of a control volume, the heat capacity changes but stored thermal energy is preserved before the next conductive update, avoiding spurious heating/cooling from a changing denominator.

External thermal boundaries match the gas boundary's `sealed` (zero heat flux) or `reservoir` (fixed ambient temperature) mode. A boundary flux is evaluated over the half-cell distance. Each internal heat flow is applied antisymmetrically to the two cells; the implicit matrix is symmetric positive definite with a Jacobi-preconditioned CG solve and verified true residual.

The closed-boundary excess-energy budget is:

    E_stored(t+dt) - E_stored(t) = ΔM_ice Ls.

With heat exchange at open boundaries:

    E_stored(t+dt)-E_stored(t) + Q_exported = ΔM_ice Ls.

`thermalLatentInjectedJ`, `thermalBoundaryHeatOutJ` and `thermalEnthalpyJ[i]` make the balance testable. The solver rejects non-finite/unconverged or out-of-Ice-Ih-range temperature fields before updating any heat state.

## 3. How the heat changes ice growth

The coupling is sequential within each growth timestep (operator splitting):

1. Solve the common gas field using the **previous** local ice-surface saturation concentration.
2. Count each gas→ice mass transfer once.
3. Update the five separate physical ice mass/fraction fields and displace gas inventory conservatively.
4. Inject the deposited latent energy into the accepting ice voxels.
5. Implicitly conduct heat through ice and gas; account ambient heat flux.
6. Recompute `c_eq,i=saturationMassDensityIce(Tinf+ΔT_ice,i)`.
7. Use the updated `c_eq,i` for the next vapor solve's exposed attachment faces.

The gas solve has a condensation-only active set: a surface whose local temperature raises its equilibrium above available vapor stops receiving ice deposition; **negative deposition is not interpreted as modeled sublimation**. Within the monotone active set it recomputes the vapor field rather than truncating a negative flux after solving an inconsistent equation. Phase-change coefficient, temperature-dependent k and D and within-step coupled nonlinear iterations are still missing. For further maturity, the next implementation must test timestep splitting convergence and actual joint thermal/vapor strong coupling.

Thermal coupling is **on by default in the separate Model-B laboratory** (`app/layered-lab.html`) with an explicit `No heat · control` setting. The legacy observatory continues using its original 2.5-D heat model; this new thermal closure is not backported to it.

The separate lab displays five independent ice occupations, actual deposition by layer, maximum ΔT, energy residual and real temperature-field pixels from the computed heat arrays. The visual gain selector only changes ice opacity, not deposition or temperature.

## 4. Measured reproducibility requirements

Numerical suite `app/tests/layered-thermal.test.mjs` asserts:

- Closed finite-volume heat diffusion conserves injected deposition latent energy.
- Open heat reservoirs export energy, with explicit balance to deposited latent input.
- Lateral and vertical thermal symmetry on an isotropic small hex grid.
- Warming raises the local ice/vapor equilibrium density.
- Four sequential five-slab vapor+ice+thermal updates retain the water ledger and latent heat ledger.
- With the same source and initial ice/vapor, thermal feedback should not cause artificial extra ice deposition.
- Failure to converge does not mutate heat state.
- Historic hex/Cartesian numerical fingerprints remain unchanged.

The exact **physical thickness of 0.75 µm per z slab is only a numerical discretization candidate**. Nine-slab confinement, boundary placement and heat conductivity approximations require z-convergence, proper ice kinetics and larger-domain validation before simulation can claim microphotographic internal ridge or stellar branch accuracy. The new feature is not molecularly resolved, even when high-resolution pixels show branches.

## 5. Following the user's 15-point implementation checklist

**Completed:** 1 frozen baseline; 2 shared finite-volume vapor ledger; 3 five independent ice fractions and mass-conserving moving coarse cells; **4 heat transport and temperature feedback (experimental, small grid)**.

**Next:** extend the tested kinetics on moving interfaces with consistently exposed basal/prism face areas, improve subcell ice surface reconstruction, then verify convergence as vertical halo/Δz and 300-µm xy resolution change. Later integrate true five-layer thickness contours into the existing main optical 3D viewer, retaining the legacy solver as an A/B control.
