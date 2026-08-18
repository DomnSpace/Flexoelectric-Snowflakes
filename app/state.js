export function createInitialState() {
  return {
    clock: { t: 0, running: false, speed: 1, dt: 0.02 },
    ambient: {
      T: 298.15,
      supersaturationPct: 0.45,
      pressurePa: 101325,
    },
    aerosol: {
      dryDiameter: 80e-9,
      kappaHyg: 0.30,
      wetDiameter: 82e-9,
      surfaceTensionMode: 'iapws',
    },
    phase: {
      name: 'dry_aerosol',
      activated: false,
      activatedAt: null,
    },
    view: {
      field: 'supersaturation',
      hexaMode: 'absolute',
      selectedArm: 0,
      selectedRadiusFactor: 2.2,
    },
    diagnostics: {
      critical: null,
      saturation: null,
      rate: 0,
      budget: null,
      modes: [],
      radial: [],
    },
    history: [
      { t: 0, phase: 'dry_aerosol', label: 'dry aerosol', detail: 'Initial hygroscopic dry particle.' },
    ],
  };
}

export function pushHistory(state, phase, label, detail) {
  if (state.history.some(e => e.phase === phase)) return;
  state.history.push({ t: state.clock.t, phase, label, detail });
}
