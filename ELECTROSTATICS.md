# Electrostatic experiment — numerical and physical contract

**Scope:** experimental electrostatic extension of the existing vapor / thermal / 2.5-D phase-field ice-growth solver. A **spatial potential field is solved**. The resulting electric field is not a cartoon arrow. Mechanical strain, however, is **prescribed** (not solved); the dielectric-to-growth closure is an **opt-in approximation**, not a validated snow-crystal growth law.

## 1. Equations and units

Electric potential V [volts] and field E [V/m]:

    E = -gradient(V)
    D = epsilon_0 epsilon_r E + P
    divergence(D) = 0

Therefore, in the absence of free charge,

    divergence(epsilon_r gradient(V)) = divergence(P) / epsilon_0

Here, epsilon_0 = 8.8541878128e-12 F/m. The ice phase fraction determines an effective spatial dielectric response,

    epsilon_r(x,y) = epsilon_air + phi_ice(x,y) (epsilon_ice-epsilon_air).

Normal fluxes use harmonic averages of face permittivities. The electric far-field Dirichlet condition is

    V_boundary(x,y) = -E_inf · (x,y).

The baseline epsilon_air = 1.0006 and editable illustrative ice values are 3.2 (fast dielectric response) and 90 (nominal low-frequency/static regime). The correct ice dielectric value depends strongly on temperature, frequency and material; no general-frequency constitutive fit is claimed.

There is **no** solved free-charge density, air ionization, conductivity, electric currents, electrohydrodynamic flow or bulk 3-D Maxwell system. The 2-D implementation treats out-of-plane variation as zero.

## 2. Prescribed flexoelectric experiment

Flexoelectric polarization is caused by **gradients of mechanical strain**, not by geometrical curvature alone.

    P_i = mu_ijkl partial_j(epsilon_kl)

For a deliberately simple 2-D applied-strain-gradient experiment we prescribe a constant vector G [m^-1] and use a reduced scalar effective coefficient mu:

    P = phi_ice * mu * G_vector
    |P| = mu [C/m] * G [m^-1] = [C/m²].

The default illustrative coefficient is 1.65 nC/m, broadly consistent with experimental bending measurements on ultrapure ice around 237 K, but not a calibrated scalar tensor law for every crystallographic direction or temperature. The prescribed uniform strain gradient is a laboratory forcing scenario, not a mechanical prediction from the snowflake geometry. Variation of phi produces interface bound-charge effects in the Poisson solution. When G = 0 or mu = 0 this source is exactly zero.

Experimental context:
- X. Wen et al. (2025), *Flexoelectricity and surface ferroelectricity of water ice*, Nature Physics 21, 1587–1593, https://doi.org/10.1038/s41567-025-02995-6
- X. Wen et al. (2025), *Streaming flexoelectricity in saline ice*, Nature Materials 24, 1533–1537, https://doi.org/10.1038/s41563-025-02332-5. Saline streaming enhancement should **not** be substituted for the pure-ice coefficient.

## 3. Opt-in fixed-field dielectric growth closure

An optional, deliberately limited equilibrium approximation treats dielectric stabilization of deposited ice in a locally fixed electric field as

    delta_g_e = -(1/2) epsilon_0 (epsilon_ice - epsilon_air) |E|^2 [J/m³]

and a dimensionless correction to interfacial equilibrium concentration is

    delta ln(c_eq) =
       delta_g_e * (m_H2O / rho_ice) / (k_B T_inf).

The updated interface concentration then enters the **existing** kinetic/diffusive resistance model:

    c_eq,electric = c_eq,neutral * exp(delta ln(c_eq))
    sigma_drive = max(0, c_gas/c_eq,electric - 1).

The model does **not** multiply the final growth velocity by an arbitrary field factor. With coupling disabled, dielectric calculations are diagnostic only; if field and gradient are zero, baseline growth is identical even when coupling is enabled.

**Important limitations:** the fixed-field energy-density model does not resolve the correct interfacial electric thermodynamic ensemble, Maxwell stress/traction, molecular adsorption orientation, contact charging, water vapor polarization dynamics, plasma ionization or surface premelting. It is not a measured ice-growth constitutive law. Exponents are hard-limited to -0.10..0.10 as a mathematical safeguard; values approaching that guardrail signal that the experimental closure is outside an appropriate parameter regime. Typical few-MV/m fields are expected to produce tiny changes in equilibrium compared with vapor supersaturation.

## 4. Solver cycle

At model start: assign ice/vapor dielectric contrast and imposed P; initialize potential from the far-field linear solution, and relax ∇·D = 0.

On each growth step:

1. Relax latent-heat and vapor transport as before.
2. Compute interface capillarity/temperature equilibrium and, only if enabled, add the stored electric equilibrium correction.
3. Deposit mass, update phase fractions, basal thickness and latent source.
4. Re-solve electric potential for the updated ice morphology so spatial potential/field maps are synchronized with the displayed crystal.
5. Capture the electric field, potential, correction, residual and polarization diagnostics along with the neutral physical fields.

The explicit solver reports RMS local voltage correction as a convergence diagnostic. This is **not** normalized Gauss-law error in physical C/m³.

## 5. Regression and falsifiability

Run from repo root:

    node app/tests/electric25d.test.mjs
    node app/tests/run-all.mjs

Tests must verify:
- a uniform dielectric and imposed boundary potential reproduce exact uniform E;
- zero E and zero G reproduce the neutral growth solver within floating-point tolerance;
- imposed flexoelectric polarization generates signed dipolar bound-charge potential;
- field calculations alone do not affect ice growth;
- enabling the optional dielectric equilibrium closure changes deposition with the expected sign at a stated, high laboratory field;
- no nonfinite potential or field and no unreported random forcing.

The browser smoke test also exercises electric controls and heatmap selection.

## 6. Next physics milestones

Before treating the result as a predictive electromechanical snowflake model, add a resolved elastic deformation/strain tensor satisfying force balance, a temperature/angle dependent ice flexoelectric tensor, validated dielectric dispersion, charge transport and Maxwell traction at a moving boundary. Then validate growth trends against published experimental electric-needle/snow-crystal measurements.
