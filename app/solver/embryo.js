import { CONSTANTS } from './thermo.js';
import { supercooledWaterProperties } from './supercooled-eos.js';
import { iceIhProperties } from './ice-ih-eos.js';

/**
 * Classical-nucleation-theory reference layer.
 *
 * IMPORTANT: this is a diagnostic/reference model, not a replacement for the
 * empirical Koop homogeneous nucleation rate. The two are intentionally shown
 * side by side.
 *
 * Default gamma_il is a source-tagged reference value near melting, 0.029 J m^-2,
 * motivated by Espinosa et al. (2013) / related TIP4P-family estimates. It is
 * explicitly configurable and should not be treated as universal.
 */

export const GAMMA_IL_REFERENCE = Object.freeze({
  value: 0.029,
  uncertainty: 0.003,
  units: 'J m^-2',
  status: 'reference-parameter',
  source: 'Espinosa et al. 2013, JACS 135, 6108–6111; extrapolated near coexistence',
});

export function solutionWaterChemicalPotentialSpecific({ T, pPa, aw }) {
  const liquid = supercooledWaterProperties(T, pPa);
  const gSolutionWater = liquid.g + CONSTANTS.RV * T * Math.log(Math.max(aw, 1e-15));
  return { liquid, gSolutionWater };
}

export function iceLiquidDrivingForce({ T, pPa, aw }) {
  const { liquid, gSolutionWater } = solutionWaterChemicalPotentialSpecific({ T, pPa, aw });
  const ice = iceIhProperties(T, Math.min(pPa, 210e6));
  const deltaMuMass = gSolutionWater - ice.g; // J kg^-1; positive means ice is lower Gibbs state
  const deltaGv = deltaMuMass * ice.rho;       // J m^-3, spherical CNT bulk driving magnitude
  return {
    T,
    pPa,
    aw,
    liquid,
    ice,
    gSolutionWater,
    gIce: ice.g,
    deltaMuMass,
    deltaGv,
    iceFavored: deltaMuMass > 0,
    referenceConsistency: 'liquid G12-15 and Ice Ih R10-06 use IAPWS triple-point/IAPWS-95 reference conventions',
  };
}

export function cntCriticalEmbryo({ T, pPa, aw, gammaIl = GAMMA_IL_REFERENCE.value }) {
  const drive = iceLiquidDrivingForce({ T, pPa, aw });
  if (!(drive.deltaGv > 0) || !(gammaIl > 0)) {
    return {
      ...drive,
      gammaIl,
      criticalRadius: Infinity,
      barrierJ: Infinity,
      barrierKBT: Infinity,
      criticalVolume: Infinity,
      criticalMass: Infinity,
      criticalMolecules: Infinity,
      status: 'ice_not_thermodynamically_favored',
    };
  }
  const criticalRadius = 2 * gammaIl / drive.deltaGv;
  const barrierJ = 16 * Math.PI * gammaIl ** 3 / (3 * drive.deltaGv ** 2);
  const criticalVolume = 4 * Math.PI * criticalRadius ** 3 / 3;
  const criticalMass = criticalVolume * drive.ice.rho;
  const criticalMolecules = criticalMass / CONSTANTS.MW * 6.02214076e23;
  const barrierKBT = barrierJ / (1.380649e-23 * T);
  return {
    ...drive,
    gammaIl,
    criticalRadius,
    barrierJ,
    barrierKBT,
    criticalVolume,
    criticalMass,
    criticalMolecules,
    status: 'cnt_reference_ready',
  };
}

export function postCriticalSeedFromCNT({
  T,
  pPa,
  aw,
  totalWaterMass,
  gammaIl = GAMMA_IL_REFERENCE.value,
  radiusFactor = 1.05,
}) {
  const cnt = cntCriticalEmbryo({ T, pPa, aw, gammaIl });
  if (!Number.isFinite(cnt.criticalRadius)) {
    return { cnt, iceMass: 0, radius: 0, volume: 0, cappedByAvailableWater: false };
  }
  const radius = cnt.criticalRadius * Math.max(1, radiusFactor);
  const volume = 4 * Math.PI * radius ** 3 / 3;
  const requestedIceMass = volume * cnt.ice.rho;
  const iceMass = Math.min(Math.max(0, totalWaterMass), requestedIceMass);
  const actualVolume = iceMass / cnt.ice.rho;
  const actualRadius = actualVolume > 0 ? Math.cbrt(3 * actualVolume / (4 * Math.PI)) : 0;
  return {
    cnt,
    iceMass,
    radius: actualRadius,
    volume: actualVolume,
    requestedIceMass,
    cappedByAvailableWater: iceMass < requestedIceMass,
    radiusFactor,
  };
}

export function phaseEnthalpyDifference({ T, pPa }) {
  const liquid = supercooledWaterProperties(T, pPa);
  const ice = iceIhProperties(T, Math.min(pPa, 210e6));
  const latentHeat = liquid.h - ice.h; // J kg^-1 released for liquid -> ice at same T,p
  return { liquid, ice, latentHeat };
}
