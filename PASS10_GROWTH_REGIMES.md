# Pass 10 — Observable regimes, explicit experiments and cloud histories

## Seven deliverables

1. **Temperature-dependent cloud conditions**. The atmospheric panel computes the saturation line over liquid water expressed as ice supersaturation from the existing Murphy–Koop saturation-pressure functions. It also shows RH over liquid water. This is a reference for common mixed-phase clouds, NOT a hard law prohibiting transient/laboratory supersaturation above liquid water and NOT a morphology boundary.
2. **Candidate regime library**. Compact plate (−15 °C, 3% over ice), sectored plate (−15 °C, 9%), stellar candidate (−15 °C, 14%), warm plate (−2 °C, 1.1%), four-stage growth history. Observed needle/column (−5 °C, 4%) and capped column are listed **unsupported** because the live triangular 2.5-D solver cannot predict true column tips or caps. The previous 18% at −5 °C was 112.4% RH over liquid water, not a standard steady-cloud prescription.
3. **Growth-path experiment**. A single evolving ice crystal at −15 °C sees 3% initially, 13% from 15 s, 5% from 45 s and 14% from 65 s, with a **3 s linear ramp** into each new setpoint. A timeline is real physics input: it never resets the crystal or draws target arms. Its selected time intervals are demonstrator settings, not measured typical cloud residence times; the finite 300 µm domain may stop the simulation before the last stage.
4. **Bounds and provenance**. Atmospheric T −35..−1 °C, ice supersaturation 0..35%, pressure 60..110 kPa, instrument min/max plus a water-saturation line that depends on temperature. The kinetic sixfold amplitude, angle, nucleation barriers, basal prefactor and simple SDAK have **uncalibrated demonstrator bounds** and are not labeled as observational validity intervals.
5. **Corrected hidden lever**. \`basalSigma0\` was previously present in preset definitions but not actually consumed when loading; it now has an explicit UI slider and is passed to the basal kinetics. Every selected preset stages a full parameter object (electric flags reset off, seed, domain, local height, capillarity, σ0, SDAK, pressure, etc.) rather than silently inheriting the previous run's values.
6. **Better startup**. An experiment-card quartet, a versioned summary, pressure/saturation readout, condition lever groups and an explicit two-step workflow: \`Stage preset\` prepares controls **without restarting** the existing ice, and \`Run experiment / restart\` creates a new simulation. Unsupported columns are informative but cannot be staged as false 2.5-D morphologies.
7. **Versioning and reproducibility**. \`app/solver/experiment-regimes.js\` owns a declarative schema v1, immutable candidate presets, a cloneable per-run history with ramped stage IDs, environmental conditions captured in each solver snapshot and extra CSV columns. Custom changes discard the pending preset timeline rather than leaving invisible future forcing. An **Experiment JSON** export additionally records the *applied* complete physical/numerical configuration, all future stages (even if an early domain stop prevents reaching them), open-boundary status, and experimental provenance independently of any pending slider changes.

## Thermodynamic checks

The existing Murphy–Koop formulas determine

    sigma_i = e / e_sat,ice(T) − 1
    sigma_i,water-saturation(T) = e_sat,water(T) / e_sat,ice(T) − 1
    RH_water = 100% × (1 + sigma_i) / (1 + sigma_i,water-saturation)

These numbers are **thermodynamic**, not calibrated morphology boundaries. At −15 °C, saturation over liquid water corresponds to approximately 15.74% ice supersaturation; 14% ice supersaturation is around 98.5% RH over liquid water. At −5 °C, the water line corresponds to ~4.98% supersaturation over ice. The UI warns above the liquid-water saturation line but allows laboratory/transient conditions.

## Environmental evolution without rebasing artefacts

The original hex transport scalar is **dimensionless relative to the far-field ice saturation reference**:

    u_i(t) = c_v,i / c_sat,ice(T_inf(t)) − 1

Naively modifying T_inf while preserving u_i changes every existing gas-cell absolute vapor concentration instantaneously, inventing an apparent gas source. \`applyGrowthEnvironment\` instead conserves the *existing absolute local vapor density at the instant of the reference switch*:

    c_v,i,old = c_sat,ice(T_old) * (1+u_i,old)
    u_i,new = c_v,i,old / c_sat,ice(T_new) − 1

It then resets the **Dirichlet boundary** to the new externally prescribed far-field gas concentration. That boundary is an OPEN environment source/sink: preserving internal c_v during the rebase **does not** imply zero net gas exchange over subsequent relaxation steps.

The existing absolute thermal anomaly \`DeltaT\` remains unchanged; the heat solver subsequently relaxes it with the new ambient T. If the optional height model is active and T changes, ice top/bottom heights are density-rescaled \`rho_old/rho_new\` to maintain **the existing recorded ice mass** during a change in thermodynamic density reference. This is a reduced-order representation of thermal volume change, **not solved thermal expansion mechanics**, crystallographic deformation, or resolved 3D vapor flow.

The solver still uses the original post-seed open-reservoir mass/latent heat model; no new energy/ice mass is produced purely by rendering a growth-history stage.

## Experimental regime truth

The observational Nakaya habit diagram has broad temperature/humidity tendencies: plates near −2 °C, columns near −5 °C, plates near −15 °C. A specific crystal form also depends on total ice supersaturation, pressure, kinetic surface structure, turbulence/ventilation, impurities, latent heating and atmospheric history. Neither the model's stage schedule nor the facet coefficients is a quantitative fit to those observations. True needles, capped columns and 3-D ridge macrosteps remain outside the current live solver.

The history profile is a scientific *intervention test*: compact core → arm candidate → lower-supersaturation facets → renewed branching, not a command to produce those shapes. If the first stage's slow/fast growth differs from expectations, the renderer displays what the solver actually generates.

## Verification

    node app/tests/experiment-regimes.test.mjs
    node app/tests/run-all.mjs

The new tests verify: correct cloud thermodynamics; unsupported 3D habits; zero hidden preset kinetics; ramp timing and stage IDs; exact internal absolute-vapor rebase, separate open external boundary, and unchanged state phi/time/latent-energy at environmental transitions; consistency of the local-height ice-side mass under a temperature change; and deep-copied replay history.

The Chromium/WebGL smoke also verifies: default 14%/−15 °C, dynamic saturation comparison; staging without restart; eventual application to actual solver; versioned history & v2 snapshots for the history preset; editing clears the timeline; old Cartesian/electric and full-density height modes still work.

## Sources for physical interpretation

- Murphy, D. M. & Koop, T. (2005). Review of the vapour pressures of ice and supercooled water for atmospheric applications. *Q. J. R. Meteorol. Soc.* 131, 1539–1565, https://doi.org/10.1256/qj.04.94
- Libbrecht, K. G. (2005). The physics of snow crystals. *Reports on Progress in Physics* 68, 855–895, https://doi.org/10.1088/0034-4885/68/4/R03
- For ice attachment/temperature and SDAK limitations see the existing repository's \`SCENARIO_VALIDITY.md\`, \`ANISOTROPY_PHASE_STUDY.md\`, \`PASS8_LOCAL_HEIGHT.md\`.

## Further work

This pass does not implement measured attachment kinetics or a 3D cloud-parcel model. To correlate crystal habits with environmental trajectories quantitatively, add real measured T-specific basal/prism barriers, ice-surface diffusion and step/terrace kinetics, quantitative experimental morphology comparisons, and resolution/convergence studies at fixed physical domain size.
