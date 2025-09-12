You are a senior creative-coding engineer. Build a single-file HTML app (no bundlers) that simulates and visualizes a 3D snowflake whose arms wiggle and grow by flexoelectric ice accretion. Deliver one self-contained index.html with inline JS/CSS that runs in any modern browser.

Goals

3D visualization with Three.js.

A physics loop that updates arm bending, tip growth (length), and plating (arm width).

Flexoelectric coupling: strain-gradient → polarization → electric torque/accumulation → growth bias.

A control panel (dat.GUI or Tweakpane) exposing all key parameters (“levers”).

Smooth performance at 60 FPS on a laptop; numerical stability (clamps).


Model (keep simple but physically flavored)

Represent the snowflake as 6 primary arms, each a polyline of nodes:

Node state per arm segment s: position x, tangent t, curvature κ, width w.

Wiggling: damped Euler–Bernoulli-like single or few modes per arm:
Ä + 2ζω Ȧ + ω² A = f_gust(t); curvature at tip κ_tip ∝ A.

Flexoelectric polarization (effective scalar): P = μ_eff * κ_tip.

Electric field: uniform vector E.
Electric torque on arm axis using tip dipole p_tip = P * tip_area * lever, torque Tz = p × E → tiny rotation of whole flake (overdamped).

Growth laws:

Tip length growth: dL/dt = k0 * max(0, S - S*) * Φ, with S* = S0 + cκ|κ_tip| (Gibbs–Thomson-like) and Φ = 1 + α_E|E·n_tip| + α_P|P|.

Plating (width growth): dw/dt = k_plate * max(0, S - S_w) * ( β_E|E·n_seg| + β_P|P| + β_base ), applied along the outer third of each arm.

Branching: when |κ_tip| > κ_branch and distance since last branch > Δs_branch, spawn a short side branch at ±90° of local tangent; inherits reduced growth coefficients.


Ambient flow forcing: colored noise OU process for f_gust.


Controls (add a collapsible pane “Levers”)

Expose sliders/inputs with sane ranges & tooltips:

Field: E_mag (0–5e3 V/m), E_theta (0–180°), E_phi (−90–90°) for 3D direction.

Flexo: μ_eff (0–5e-6 C/m).

Supersaturation & kinetics: S (0–0.2), S0, k0 (0–1e-4 m/s).

Electro-enhancement: α_E (0–5e-4), α_P (0–5e-2).

Plating: k_plate, S_w, β_E, β_P, β_base, plating extent (fraction of arm length).

Branching: κ_branch, Δs_branch, side-branch scale factor.

Wiggle: ω (Hz), ζ (0–0.5), gust strength, gust correlation time.

Geometry: initial arm length, initial width, node count per arm.

Time: dt, sim speed multiplier.

View: toggle axes, wireframe, auto-rotate.

Run controls: Play/Pause/Step/Reset; “Randomize seeds”; “Export PNG”.


Rendering

Use Three.js (CDN from unpkg).

Geometry per arm: tubular mesh or thin ribbon whose radius uses current width w(s). Primary arms + any side branches.

Center plate: small hexagonal prism; optional faint disk.

Shader or material: basic/Phong is fine; add subtle transparency/alpha falloff with width.

Visual hints:

Color arm tips by |P| (e.g., brighter for higher polarization).

Optionally draw a faint world-space arrow for E.

Show tiny dots when a branch spawns.



Architecture

initScene(), initGUI(), resetState(), step(dt), rebuildMeshes(), updateArmGeometry(arm).

Keep state in JS objects:

const params = { E_mag: 1500, E_theta: 45, E_phi: 10, mu_eff: 1e-6, S: 0.08, S0: 0.02,
  k0: 2e-5, alpha_E: 1e-4, alpha_P: 1e-2, k_plate: 5e-6, S_w: 0.01,
  beta_E: 0.4, beta_P: 0.6, beta_base: 0.1, kappa_branch: 80, ds_branch: 1e-4,
  omega: 12, zeta: 0.12, gustSigma: 2e-8, gustTau: 0.25, dt: 0.016, speed: 1.0,
  armNodes: 24, L0: 6e-4, w0: 3e-5, platingExtent: 0.35, seed: 1
};
const snowflake = { arms: [ /* 6 arms, each with nodes, branches array, modal A,A_dot */ ],
                    angleYaw: 0, anglePitch: 0 };

Numerical safety: clamp A, A_dot, dL, dw per step; rebuild geometry only when dirty; GC-friendly pools.


UI/UX

Left: scene; right: control pane.

Buttons: Play, Pause, Step, Reset, Randomize.

Keyboard: Space=Play/Pause, R=Reset, S=Step, L=toggle auto-rotate.

A small status HUD: FPS, mean |κ_tip|, mean growth rates, branch count.


Acceptance tests (self-check before finalizing)

1. With S=0.00, no growth; with arms extend and get visibly wider near tips over 10–20 seconds.


2. Increasing E_mag increases anisotropy: arms aligned with E grow faster (both length and plating).


3. Increasing μ_eff at fixed E biases growth via flexo term (tips with larger curvature grow faster; see brighter tip color).


4. Setting a low κ_branch produces side branches; raising it suppresses them.


5. Wiggle: raising gustSigma or changing ω visibly changes tip motion amplitude.


6. Reset returns to initial hexapod; Randomize changes noise realization.



Delivery

Produce a single index.html with:

<script src="https://unpkg.com/three@0.158.0/build/three.min.js"></script>

<script src="https://unpkg.com/three@0.158.0/examples/js/controls/OrbitControls.js"></script>

<script src="https://unpkg.com/tweakpane/dist/tweakpane.min.js"></script> (or dat.GUI)

All logic inline; commented sections; no external assets.

Start the sim automatically; include a helpful README comment at top with controls.


Stretch (if time allows): simple CSV export of time series (avg length, avg width, branch count, mean |P|).
