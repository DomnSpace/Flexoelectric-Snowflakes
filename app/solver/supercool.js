import {
  CONSTANTS,
  kohlerTerms,
  saturationState,
  saturationVaporPressureIceMK,
  saturationVaporPressureWaterMK,
  waterActivityKappa,
} from './thermo.js';

/**
 * Supercooled-liquid diagnostics connecting Milestone 1 activation to Milestone 2 ice nucleation.
 *
 * Thermodynamic identities:
 *   a_w^i(T) = p_i,sat(T) / p_w,sat(T)
 *   Delta a_w = a_w - a_w^i
 *   Delta mu_(liq->ice drive) = R T ln(a_w / a_w^i)
 *
 * The water-activity coordinate follows Koop et al. (2000). We deliberately do NOT
 * evaluate an empirical homogeneous nucleation rate J here; this module exposes the
 * thermodynamic coordinate and exact Poisson hazard algebra so a chosen, documented
 * J parameterization can be plugged in during Milestone 2 without being hidden.
 */

export const KOOP_ACTIVITY_WINDOW = Object.freeze({ min: 0.26, max: 0.34 });
export const KOOP_REFERENCE_T_MAX = 240; // K: conservative boundary for the classic solution-droplet parameterization

export function iceEquilibriumWaterActivity(T) {
  return saturationVaporPressureIceMK(T) / saturationVaporPressureWaterMK(T);
}

export function liquidIceChemicalPotentialDrive(T, aw) {
  const awIce = iceEquilibriumWaterActivity(T);
  return {
    awIce,
    deltaAw: aw - awIce,
    deltaMuLiquidMinusIce: CONSTANTS.R * T * Math.log(Math.max(aw, 1e-15) / Math.max(awIce, 1e-15)),
  };
}

export function supercooledState({ T, D, Dd, kappaHyg, ambientSw }) {
  const sat = saturationState(T, ambientSw);
  const terms = kohlerTerms(D, Dd, kappaHyg, T);
  const aw = waterActivityKappa(D, Dd, kappaHyg);
  const ice = liquidIceChemicalPotentialDrive(T, aw);
  const surfaceVaporPressure = terms.Seq * sat.esw;
  const surfaceSi = surfaceVaporPressure / sat.esi;
  const bulkSolutionSi = aw * sat.esw / sat.esi;
  const belowMelting = T < 273.15;
  const iceFavored = ice.deltaMuLiquidMinusIce > 0;
  const metastableLiquid = belowMelting && iceFavored;
  const inKoopActivityWindow = ice.deltaAw >= KOOP_ACTIVITY_WINDOW.min && ice.deltaAw <= KOOP_ACTIVITY_WINDOW.max;
  const inClassicKoopTemperatureRange = T <= KOOP_REFERENCE_T_MAX;

  let regime = 'liquid_stable';
  if (belowMelting && !iceFavored) regime = 'solution_liquid_stable_against_ice';
  if (metastableLiquid) regime = 'supercooled_metastable_liquid';
  if (metastableLiquid && inKoopActivityWindow && inClassicKoopTemperatureRange) regime = 'homogeneous_nucleation_coordinate_window';

  return {
    T,
    aw,
    awIce: ice.awIce,
    deltaAw: ice.deltaAw,
    deltaMuLiquidMinusIce: ice.deltaMuLiquidMinusIce,
    belowMelting,
    iceFavored,
    metastableLiquid,
    inKoopActivityWindow,
    inClassicKoopTemperatureRange,
    regime,
    ambientSw: sat.Sw,
    ambientSi: sat.Si,
    ambientWaterSupersaturationPct: (sat.Sw - 1) * 100,
    ambientIceSupersaturationPct: (sat.Si - 1) * 100,
    surfaceSeq: terms.Seq,
    surfaceSi,
    bulkSolutionSi,
    esw: sat.esw,
    esi: sat.esi,
  };
}

/** Exact stochastic conversion once a volumetric homogeneous nucleation rate J is supplied. */
export function freezingProbabilityFromRate(J_m3_s, diameter, dt) {
  const V = Math.PI * diameter ** 3 / 6;
  const hazard = Math.max(0, J_m3_s) * V * Math.max(0, dt);
  return 1 - Math.exp(-hazard);
}

/** J needed for a chosen one-particle freezing probability over dt. */
export function nucleationRateForProbability(probability, diameter, dt) {
  const p = Math.min(1 - 1e-15, Math.max(0, probability));
  const V = Math.PI * diameter ** 3 / 6;
  if (!(V > 0) || !(dt > 0)) return Infinity;
  return -Math.log(1 - p) / (V * dt);
}

export function makeSupercoolPhaseBoundary({ Tmin = 200, Tmax = 278, points = 260 } = {}) {
  const data = [];
  for (let i = 0; i < points; i++) {
    const T = Tmin + (Tmax - Tmin) * i / (points - 1);
    data.push({
      T,
      Tc: T - 273.15,
      awIce: iceEquilibriumWaterActivity(T),
      deltaAwLow: iceEquilibriumWaterActivity(T) + KOOP_ACTIVITY_WINDOW.min,
      deltaAwHigh: iceEquilibriumWaterActivity(T) + KOOP_ACTIVITY_WINDOW.max,
    });
  }
  return data;
}

/** Dense T-a_w grid for a 2-D phase-state surface / contour renderer. */
export function makePhaseStateSurface({ Tmin = 200, Tmax = 278, awMin = 0.55, awMax = 1.0, nT = 72, nAw = 64 } = {}) {
  const cells = [];
  for (let i = 0; i < nT; i++) {
    const T = Tmin + (Tmax - Tmin) * i / (nT - 1);
    const awIce = iceEquilibriumWaterActivity(T);
    for (let j = 0; j < nAw; j++) {
      const aw = awMin + (awMax - awMin) * j / (nAw - 1);
      const deltaAw = aw - awIce;
      const deltaMu = CONSTANTS.R * T * Math.log(aw / awIce);
      cells.push({
        T,
        Tc: T - 273.15,
        aw,
        awIce,
        deltaAw,
        deltaMu,
        metastable: T < 273.15 && deltaMu > 0,
        koopWindow: T <= KOOP_REFERENCE_T_MAX && deltaAw >= KOOP_ACTIVITY_WINDOW.min && deltaAw <= KOOP_ACTIVITY_WINDOW.max,
      });
    }
  }
  return { cells, nT, nAw, Tmin, Tmax, awMin, awMax };
}

export function supercoolSnapshot({ t, T, supersaturationPct, D, Dd, kappaHyg }) {
  const ambientSw = 1 + supersaturationPct / 100;
  return { t, D, ...supercooledState({ T, D, Dd, kappaHyg, ambientSw }) };
}
