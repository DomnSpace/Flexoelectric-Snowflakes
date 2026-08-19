import { supercooledState } from './supercool.js';

/**
 * Milestone 2 nucleation engine.
 *
 * Homogeneous pathway:
 * Koop et al. (2000) water-activity parameterization, as commonly written
 *   log10 J_hom[cm^-3 s^-1] = -906.7 + 8502 Δaw - 26924 Δaw^2 + 29180 Δaw^3
 * with fitted range 0.26 < Δaw < 0.34.
 *
 * Immersion pathway:
 * generic ABIFM form from Knopf & Alpert (2013)
 *   log10 J_het[cm^-2 s^-1] = m Δaw + c
 * but m,c are material-specific and therefore MUST be supplied explicitly.
 *
 * Stochastic execution is exact for piecewise-constant rate: a single exponential
 * threshold H*=-ln(U) is sampled once, then cumulative hazard H=∫λdt advances.
 */

export const KOOP2000 = Object.freeze({
  deltaAwMin: 0.26,
  deltaAwMax: 0.34,
  coefficients: Object.freeze([-906.7, 8502, -26924, 29180]),
  sourceUnits: 'cm^-3 s^-1',
});

export function koopHomogeneousRate(deltaAw, { allowExtrapolation = false } = {}) {
  const valid = deltaAw > KOOP2000.deltaAwMin && deltaAw < KOOP2000.deltaAwMax;
  if (!valid && !allowExtrapolation) {
    return {
      valid: false,
      deltaAw,
      Jcm3s: null,
      Jm3s: null,
      log10Jcm3s: null,
      status: deltaAw <= KOOP2000.deltaAwMin ? 'below_fit_range' : 'above_fit_range',
    };
  }
  const [b0,b1,b2,b3] = KOOP2000.coefficients;
  const log10Jcm3s = b0 + b1*deltaAw + b2*deltaAw**2 + b3*deltaAw**3;
  const Jcm3s = 10 ** log10Jcm3s;
  return {
    valid,
    deltaAw,
    log10Jcm3s,
    Jcm3s,
    Jm3s: Jcm3s * 1e6,
    status: valid ? 'within_fit_range' : 'extrapolated',
  };
}

export function abifmImmersionRate(deltaAw, calibration) {
  if (!calibration || !Number.isFinite(calibration.m) || !Number.isFinite(calibration.c)) {
    return { valid:false, status:'material_calibration_required', Jcm2s:null, Jm2s:null, log10Jcm2s:null };
  }
  const log10Jcm2s = calibration.m * deltaAw + calibration.c;
  const Jcm2s = 10 ** log10Jcm2s;
  return {
    valid:true,
    status:'material_calibration_supplied',
    calibration,
    log10Jcm2s,
    Jcm2s,
    Jm2s:Jcm2s * 1e4,
  };
}

export function volumetricHazardRate(Jm3s, volumeM3) {
  return Math.max(0, Jm3s ?? 0) * Math.max(0, volumeM3 ?? 0);
}

export function surfaceHazardRate(Jm2s, areaM2) {
  return Math.max(0, Jm2s ?? 0) * Math.max(0, areaM2 ?? 0);
}

export function freezingProbabilityFromHazard(hazard) {
  return 1 - Math.exp(-Math.max(0, hazard));
}

function xorshift32(seed) {
  let x = (seed | 0) || 0x6d2b79f5;
  x ^= x << 13; x ^= x >>> 17; x ^= x << 5;
  return x >>> 0;
}

export function makeNucleationClock(seed = 1) {
  const nextSeed = xorshift32(seed);
  const u = (nextSeed + 0.5) / 4294967296;
  return {
    seed: nextSeed,
    uniform: u,
    threshold: -Math.log(Math.max(1e-16, 1-u)),
    cumulativeHazard: 0,
    crossed: false,
    crossedAt: null,
    pathway: null,
  };
}

export function advanceHazardClock(clock, hazardRatePerSecond, dt, t, pathway) {
  if (clock.crossed) return clock;
  const dH = Math.max(0, hazardRatePerSecond) * Math.max(0, dt);
  clock.cumulativeHazard += dH;
  if (clock.cumulativeHazard >= clock.threshold) {
    clock.crossed = true;
    clock.crossedAt = t;
    clock.pathway = pathway;
  }
  return clock;
}

export function homogeneousNucleationDiagnostic({ T, ambientSw, Dd, waterMass, pAir=101325, inventory }) {
  const sc = supercooledState({T, ambientSw, Dd, waterMass, pAir, inventory});
  const rate = koopHomogeneousRate(sc.deltaAw);
  const liquidVolume = sc.waterMass / sc.liquidProperties.rho;
  const hazardRate = rate.Jm3s == null ? 0 : volumetricHazardRate(rate.Jm3s, liquidVolume);
  return {
    state: sc,
    rate,
    liquidVolume,
    hazardRate,
    oneSecondProbability: freezingProbabilityFromHazard(hazardRate),
  };
}

export function immersionNucleationDiagnostic({ deltaAw, immersedAreaM2, calibration }) {
  const rate = abifmImmersionRate(deltaAw, calibration);
  const hazardRate = rate.Jm2s == null ? 0 : surfaceHazardRate(rate.Jm2s, immersedAreaM2);
  return {
    rate,
    immersedAreaM2,
    hazardRate,
    oneSecondProbability: freezingProbabilityFromHazard(hazardRate),
  };
}
