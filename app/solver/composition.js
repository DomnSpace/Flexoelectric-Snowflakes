import { CONSTANTS, maxwellMasonRadiusRate, saturationVaporPressureWaterMK, surfaceTensionWaterIAPWS } from './thermo.js';
import { supercooledWaterProperties } from './supercooled-eos.js';

const SPHERE = Math.PI / 6;

export const DEFAULT_SOLUTE = Object.freeze({
  label: 'effective dissolved solute',
  molarMass: 0.13214,
  dryDensity: 1770,
  vantHoff: 3.0,
  osmoticCoefficient: 1.0,
  activityModel: 'ideal_mole_fraction',
});

export function soluteInventoryFromDryParticle(Dd, solute = DEFAULT_SOLUTE) {
  const dryVolume = SPHERE * Dd ** 3;
  const dryMass = solute.dryDensity * dryVolume;
  const soluteMoles = dryMass / solute.molarMass;
  return { ...solute, dryVolume, dryMass, soluteMoles };
}

export function compositionFromWaterMass(waterMass, inventory) {
  const nWater = Math.max(0, waterMass) / CONSTANTS.MW;
  const effectiveSoluteMoles = inventory.vantHoff * inventory.osmoticCoefficient * inventory.soluteMoles;
  let aw = 0;
  if (nWater > 0) {
    aw = inventory.activityModel === 'osmotic_exponential'
      ? Math.exp(-effectiveSoluteMoles / nWater)
      : nWater / (nWater + effectiveSoluteMoles);
  }
  const totalMass = waterMass + inventory.dryMass;
  return {
    waterMass,
    waterMoles: nWater,
    soluteMoles: inventory.soluteMoles,
    effectiveSoluteMoles,
    aw,
    waterMassFraction: totalMass > 0 ? waterMass / totalMass : 0,
    soluteMassFraction: totalMass > 0 ? inventory.dryMass / totalMass : 1,
    molality: waterMass > 0 ? inventory.soluteMoles / waterMass : Infinity,
  };
}

export function effectiveKappaFromComposition(D, Dd, aw) {
  const Vd = Dd ** 3;
  const Vw = Math.max(0, D ** 3 - Vd);
  if (!(Vd > 0) || !(aw > 0)) return Infinity;
  return Vw / Vd * (1 / aw - 1);
}

export function surfaceTensionSupercooled(T) {
  const sigma = surfaceTensionWaterIAPWS(T);
  return {
    sigma,
    referenceSupported: T >= 248.15,
    note: T >= 248.15
      ? 'IAPWS R1-76(2014) metastable extrapolation within stated support.'
      : 'IAPWS surface-tension equation extrapolated below its stated supercooled support (−25 °C).',
  };
}

export function diameterFromWaterMass({ waterMass, Dd, T, pAir = 101325, inventory, iterations = 16 }) {
  const inv = inventory ?? soluteInventoryFromDryParticle(Dd);
  const surf = surfaceTensionSupercooled(T);
  let rhoGuess = 1000;
  let D = Math.cbrt(Math.max(Dd ** 3, Dd ** 3 + 6 * Math.max(0, waterMass) / (Math.PI * rhoGuess)));
  let props = null;
  for (let i = 0; i < iterations; i++) {
    const laplacePressure = 4 * surf.sigma / Math.max(D, 1e-12);
    const pLiquid = pAir + laplacePressure;
    props = supercooledWaterProperties(T, pLiquid);
    const waterVolume = Math.max(0, waterMass) / props.rho;
    const newD = Math.cbrt(Math.max(Dd ** 3, 6 * (inv.dryVolume + waterVolume) / Math.PI));
    if (Math.abs(newD - D) / Math.max(D, 1e-12) < 1e-11) { D = newD; break; }
    D = 0.5 * (D + newD);
  }
  const laplacePressure = 4 * surf.sigma / Math.max(D, 1e-12);
  const pLiquid = pAir + laplacePressure;
  props = supercooledWaterProperties(T, pLiquid);
  return { D, radius: D / 2, laplacePressure, pLiquid, waterProperties: props, surface: surf, inventory: inv };
}

export function waterMassFromWetDiameter({ D, Dd, T, pAir = 101325, inventory }) {
  const inv = inventory ?? soluteInventoryFromDryParticle(Dd);
  const surf = surfaceTensionSupercooled(T);
  const laplacePressure = 4 * surf.sigma / Math.max(D, 1e-12);
  const pLiquid = pAir + laplacePressure;
  const props = supercooledWaterProperties(T, pLiquid);
  const waterVolume = Math.max(0, SPHERE * (D ** 3 - Dd ** 3));
  return { waterMass: waterVolume * props.rho, laplacePressure, pLiquid, waterProperties: props, inventory: inv, surface: surf };
}

export function solutionDropletState({ waterMass, Dd, T, pAir = 101325, inventory }) {
  const geom = diameterFromWaterMass({ waterMass, Dd, T, pAir, inventory });
  const comp = compositionFromWaterMass(waterMass, geom.inventory);
  const { D, waterProperties: wp } = geom;
  const kelvinLn = 4 * geom.surface.sigma * CONSTANTS.MW / (CONSTANTS.R * T * wp.rho * D);
  const Seq = comp.aw * Math.exp(kelvinLn);
  const kappaEffective = effectiveKappaFromComposition(D, Dd, comp.aw);
  return { ...geom, composition: comp, kelvinLn, kelvinFactor: Math.exp(kelvinLn), Seq, kappaEffective };
}

export function integrateWaterMassStep({ waterMass, Dd, T, pAir = 101325, vaporPressure, inventory, dt, maxFraction = 0.01 }) {
  let mass = Math.max(0, waterMass);
  let remaining = dt;
  let substeps = 0;
  let last = null;
  while (remaining > 1e-12 && substeps < 20000) {
    const drop = solutionDropletState({ waterMass: mass, Dd, T, pAir, inventory });
    const Sw = vaporPressure / saturationVaporPressureWaterMK(T);
    const transport = maxwellMasonRadiusRate(drop.radius, T, Sw, drop.Seq, pAir, { rhoW: drop.waterProperties.rho });
    const massRate = 4 * Math.PI * drop.radius ** 2 * drop.waterProperties.rho * transport.rate;
    const scale = Math.max(mass, drop.inventory.dryMass * 1e-3, 1e-24);
    const h = Math.min(remaining, Math.abs(massRate) > 0 ? maxFraction * scale / Math.abs(massRate) : remaining);
    mass = Math.max(0, mass + massRate * h);
    remaining -= h;
    substeps++;
    last = { drop, transport, massRate, Sw };
  }
  const drop = solutionDropletState({ waterMass: mass, Dd, T, pAir, inventory });
  return { waterMass: mass, drop, transport: last?.transport ?? null, massRate: last?.massRate ?? 0, Sw: vaporPressure / saturationVaporPressureWaterMK(T), substeps };
}
