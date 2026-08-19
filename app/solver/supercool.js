import {
  CONSTANTS,
  saturationState,
  saturationVaporPressureIceMK,
  saturationVaporPressureWaterMK,
} from './thermo.js';
import {
  DEFAULT_SOLUTE,
  solutionDropletState,
  soluteInventoryFromDryParticle,
} from './composition.js';

/**
 * Supercooled-liquid diagnostics connecting Milestone 1 activation to Milestone 2 ice nucleation.
 *
 * Baseline liquid/ice equilibrium at low pressure:
 *   a_w^i(T) = p_i,sat(T) / p_w,sat(T)
 *
 * Pressure-adjusted bulk liquid/ice equilibrium:
 *   ln a_w^i(T,p) = ln a_w^i(T,p_ref) + (V_i - V_l) (p - p_ref)/(R T)
 *
 * Here p is the common bulk pressure of candidate liquid and ice phases. The liquid
 * molar volume comes from IAPWS G12-15. Ice molar volume is currently an explicit
 * constant-density approximation (rho_ice = 917 kg m^-3) until the IAPWS Ice-Ih EOS
 * is added. Ice-embryo interfacial Laplace pressure is NOT included yet because it
 * requires embryo geometry and gamma_il; that is Milestone 2 nucleation physics.
 */

export const KOOP_ACTIVITY_WINDOW = Object.freeze({ min: 0.26, max: 0.34 });
export const KOOP_REFERENCE_T_MAX = 240;
export const ICE_DENSITY_APPROX = 917.0;

export function iceEquilibriumWaterActivity(T) {
  return saturationVaporPressureIceMK(T) / saturationVaporPressureWaterMK(T);
}

export function pressureAdjustedIceEquilibriumWaterActivity(T, pCommon, liquidDensity, pRef = 101325) {
  const aw0 = iceEquilibriumWaterActivity(T);
  const Vi = CONSTANTS.MW / ICE_DENSITY_APPROX;
  const Vl = CONSTANTS.MW / liquidDensity;
  const lnCorrection = (Vi - Vl) * (pCommon - pRef) / (CONSTANTS.R * T);
  return {
    awReference: aw0,
    awPressure: aw0 * Math.exp(lnCorrection),
    Vi,
    Vl,
    lnCorrection,
    pressureCorrectionApproximation: 'liquid density = IAPWS G12-15; ice density fixed at 917 kg m^-3; no ice-embryo interfacial pressure yet',
  };
}

export function liquidIceChemicalPotentialDrive(T, aw, options = {}) {
  const pCommon = options.pCommon ?? 101325;
  const liquidDensity = options.liquidDensity ?? 1000;
  const pressure = pressureAdjustedIceEquilibriumWaterActivity(T, pCommon, liquidDensity, options.pRef ?? 101325);
  return {
    awIce: pressure.awPressure,
    awIceReference: pressure.awReference,
    deltaAw: aw - pressure.awPressure,
    deltaMuLiquidMinusIce: CONSTANTS.R * T * Math.log(Math.max(aw, 1e-15) / Math.max(pressure.awPressure, 1e-15)),
    pressure,
  };
}

export function supercooledState({
  T,
  ambientSw,
  Dd,
  waterMass,
  pAir = 101325,
  inventory,
  solute = DEFAULT_SOLUTE,
}) {
  const inv = inventory ?? soluteInventoryFromDryParticle(Dd, solute);
  const drop = solutionDropletState({ waterMass, Dd, T, pAir, inventory: inv });
  const sat = saturationState(T, ambientSw);
  const aw = drop.composition.aw;
  const ice = liquidIceChemicalPotentialDrive(T, aw, {
    pCommon: drop.pLiquid,
    liquidDensity: drop.waterProperties.rho,
    pRef: pAir,
  });

  const surfaceVaporPressure = drop.Seq * sat.esw;
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
    drop,
    D: drop.D,
    radius: drop.radius,
    waterMass,
    aw,
    awIce: ice.awIce,
    awIceReference: ice.awIceReference,
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
    surfaceSeq: drop.Seq,
    surfaceSi,
    bulkSolutionSi,
    esw: sat.esw,
    esi: sat.esi,
    laplacePressure: drop.laplacePressure,
    pLiquid: drop.pLiquid,
    liquidProperties: drop.waterProperties,
    composition: drop.composition,
    effectiveKappa: drop.kappaEffective,
    kelvinLn: drop.kelvinLn,
    surfaceTension: drop.surface,
    liquidIcePressureCorrection: ice.pressure,
  };
}

export function freezingProbabilityFromRate(J_m3_s, diameter, dt) {
  const V = Math.PI * diameter ** 3 / 6;
  const hazard = Math.max(0, J_m3_s) * V * Math.max(0, dt);
  return 1 - Math.exp(-hazard);
}

export function nucleationRateForProbability(probability, diameter, dt) {
  const p = Math.min(1 - 1e-15, Math.max(0, probability));
  const V = Math.PI * diameter ** 3 / 6;
  if (!(V > 0) || !(dt > 0)) return Infinity;
  return -Math.log(1 - p) / (V * dt);
}

export function makeSupercoolPhaseBoundary({ Tmin = 200, Tmax = 278, points = 260, pCommon = 101325, liquidDensity = 1000 } = {}) {
  const data = [];
  for (let i = 0; i < points; i++) {
    const T = Tmin + (Tmax - Tmin) * i / (points - 1);
    const awIce = pressureAdjustedIceEquilibriumWaterActivity(T, pCommon, liquidDensity).awPressure;
    data.push({
      T,
      Tc: T - 273.15,
      awIce,
      deltaAwLow: awIce + KOOP_ACTIVITY_WINDOW.min,
      deltaAwHigh: awIce + KOOP_ACTIVITY_WINDOW.max,
    });
  }
  return data;
}

export function makePhaseStateSurface({ Tmin = 200, Tmax = 278, awMin = 0.55, awMax = 1.0, nT = 72, nAw = 64, pCommon = 101325, liquidDensity = 1000 } = {}) {
  const cells = [];
  for (let i = 0; i < nT; i++) {
    const T = Tmin + (Tmax - Tmin) * i / (nT - 1);
    const awIce = pressureAdjustedIceEquilibriumWaterActivity(T, pCommon, liquidDensity).awPressure;
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

export function supercoolSnapshot({ t, T, ambientSw, Dd, waterMass, pAir = 101325, inventory, solute }) {
  return { t, ...supercooledState({ T, ambientSw, Dd, waterMass, pAir, inventory, solute }) };
}
