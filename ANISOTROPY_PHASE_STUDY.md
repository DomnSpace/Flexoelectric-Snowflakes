# Pass 6: disentangling kinetic phase, capillary stiffness and initial crystal orientation

## The actual source of the visible split arms

The sixfold splitting in the **triangular** growth mode cannot be attributed to its dielectric Poisson field, because `app/solver/electric25d.js` is still strictly Cartesian and disabled in hex mode. The hex solver computes a physical vapor field, thermal feedback, interface normal from the *local phase gradient*, and attachment using the factor

    alpha = clamp[alpha0(sigma) * (1 + a*cos(6(theta - phase))), 0, 1]

where `theta` is the interface's **outward unit-normal angle**. This is NOT the polar location of a crystal tip.

Rotation of a sixfold kinetic preference by 30° gives

    cos(6(theta - pi/6)) = -cos(6 theta)

so the kinetic maxima and minima exchange. Whether this suppresses tip splitting is an **empirical question** for the present PDE and its boundary, not a mathematical consequence of the trigonometric identity. The hexagonal far-field Dirichlet boundary has sixfold symmetry of its own and can imprint an apparent `m6` signature even when all actual attachment anisotropy is removed.

Reference: Barrett, Garcke & Nürnberg (2012), "Numerical computations of faceted pattern formation in snow crystal growth", Physical Review E 86, 011604, DOI 10.1103/PhysRevE.86.011604; Demange et al. (2017), "A phase field model for snow crystal growth in three dimensions", npj Computational Materials 3, 15, https://doi.org/10.1038/s41524-017-0015-1.

## Three independently visible experimental variables

1. **Kinetic orientation phase (0°, 15°, 30°)** applies only to the prism attachment rule in hex mode. Zero is the existing baseline, 15° a non-aligned control and 30° swaps maxima/minima. It does not rotate the whole domain or force six trunks geometrically.
2. **Anisotropic capillary stiffness** uses the 2-D interface stiffness `gamma+gamma''`, rather than substituting an isotropically curved interface or incorrectly modifying the electric potential. For `gamma(theta)=gamma0[1+epsilon cos 6(theta-phase)]`, `gamma+gamma''=gamma0[1-35 epsilon cos 6(theta-phase)]`. Options are epsilon = 0, 0.01 and 0.02. The underlying function rejects |epsilon|>0.025, well below the unregularized positivity limit ~1/35 = .02857. This **weak smooth-anisotropy experiment does not create sharp crystalline cusps or true faceted Wulff shapes**; stronger faceting requires missing-orientation regularization.
3. **Seed occupancy** is either the original circular, unbiased physical seed, or a true hexagonal polygon initially occupied on the same physical numerical lattice. Its vertex azimuth is independently 0°, 15° or 30°. Small seeds remain discretized on the 0.75 µm triangular lattice, so rotation may be partially aliased at the initial size.

The current environmental presets retain the original defaults (circular seed, phase=0°, surface epsilon=0) for reproducibility. Optional experiments change real model input and are included in downloaded CSV. The older Cartesian Poisson solver has no new angular controls and remains unaffected.

## Falsification / measured consequences

The numerical test `app/tests/anisotropy-phase.test.mjs` verifies:

- identity of ice phase fields with kinetic anisotropy amplitude 0 and capillary epsilon 0 when switching 0° to 30° (phase alone is not a hidden forcing);
- exact extrema swapping under 30° for positive kinetic anisotropy;
- sixfold surface-stiffness bounds and rejection of unrealistically strong unregularized epsilon;
- physical seed occupation counts differ for circular vs hexagonal seeds of the same circumradius;
- 0°/15°/30° runs for otherwise identical temperature, supersaturation, seed, lattice and timestep;
- a **measured** center-ray minus ±10° flank phase-occupancy diagnostic both on fixed laboratory axes and on each run's independently measured *sixfold morphology orientation* (including the outer 55–95% of the growth arm). Comparing 30° results only on the fixed laboratory axes would mistake a rotated intact arm for a split arm. Negative aligned tip score is a possible splitter indicator, **not** proof of mechanism or natural crystal morphology;
- opt-in nonzero surface-energy anisotropy changes numerical deposition.

The comparisons do NOT assume 30° is better, and a negative spine/flank score is not treated as a failed test. The CI logs print outcomes for every phase. This is precisely why we **must not** merely flip the angular anisotropy sign without observing the resulting growth.

## What would then produce sharper, more realistic Eiskristalle?

- **Facet/cusp regularization** rather than unbounded surface anisotropy. Facet-limited attachment must be consistently calibrated with ice Ih prism orientations, not selected by the shape of the numerical Dirichlet boundary.
- **Coupled terrace-width-dependent SDAK** rather than the current simplified six-neighbor occupancy heuristic.
- **Conserved surface diffusion and spatially varying ice thickness**; separate basal and prism mass budgets before claiming internal ridges.
- **Environment history** (temperature, ice supersaturation and optional air velocity varying in time) using physically evaluated pressure/vapor pathways. Growth rings then arise from recorded history rather than painted bands.
- **Small mass-conserving, seeded stochastic forcing** tested against a zero-noise control. +/- 1% independent noise applied directly to each vapor cell each step is *not* by itself a correct diffusion PDE or a conserved thermal perturbation model.
- **Grid rotation, mesh refinement and bounded-domain convergence tests** so sixfold numerical symmetry cannot masquerade as crystalline sixfold physics.
- Optical 5-layer plate views remain visualizations until we solve an independent local thickness / mass-transport equation.

No 3-D rendering upgrade in this pass, to isolate the physics first.
