export function createInitialState() {
  return {
    clock: { t: 0, running: false, speed: 1, dt: 0.02 },
    ambient: { T: 298.15, supersaturationPct: 0.05, pressurePa: 101325 },
    aerosol: {
      dryDiameter: 80e-9,
      kappaHyg: 0.30,
      wetDiameter: 82e-9,
      surfaceTensionMode: 'iapws',
    },
    experiment: {
      mode: 'supersat_ramp',
      startSupersaturationPct: 0.05,
      targetSupersaturationPct: 0.45,
      rampDuration: 12,
      coolingRateKPerMin: 4,
      targetTemperatureC: 5,
      driver: null,
      barrierCrossedAt: null,
    },
    supercool: {
      enabled: false,
      mode: 'cool_fixed_e',
      targetTemperatureC: -38,
      coolingRateKPerMin: 6,
      enteredBelowZeroAt: null,
      iceSupersaturatedAt: null,
      metastableAt: null,
      koopWindowAt: null,
      nucleationModel: 'none',
      selectedJ_m3_s: null,
      observationWindowS: 1,
    },
    phase: { name: 'dry_aerosol', activated: false, activatedAt: null },
    view: {
      field: 'supersaturation',
      hexaMode: 'absolute',
      selectedArm: 0,
      selectedRadiusFactor: 2.2,
      inspectIndex: -1,
      liveInspect: true,
      workMode: 'work',
      phaseSurfaceScalar: 'deltaMu',
      phaseSurfaceCursor: null,
    },
    diagnostics: {
      critical: null, saturation: null, rate: 0, budget: null,
      modes: [], radial: [], curve: null, landscape: null, equilibria: null,
      supercool: null, phaseSurface: null,
    },
    history: [
      { t: 0, phase: 'dry_aerosol', label: 'dry aerosol', detail: 'Initial dry particle before controlled humidification.' },
    ],
    trajectory: [],
    sampleAccumulator: 0,
  };
}

export function pushHistory(state, phase, label, detail) {
  if (state.history.some(e => e.phase === phase)) return;
  state.history.push({ t: state.clock.t, phase, label, detail });
}

export function pushTrajectory(state, sample, force = false) {
  const last = state.trajectory[state.trajectory.length - 1];
  if (!force && last && state.clock.t - last.t < 0.10) return;
  state.trajectory.push(sample);
  if (state.trajectory.length > 2400) state.trajectory.shift();
}
