import { compositionFromWaterMass } from './composition.js';
import { postCriticalSeedFromCNT, phaseEnthalpyDifference } from './embryo.js';

/**
 * Post-critical ice seed born from a stochastic nucleation event.
 *
 * Seed mass is now set by a source-tagged CNT reference embryo radius rather than
 * an arbitrary fraction of droplet mass. The stochastic Koop hazard remains the
 * event generator; CNT is used only to size/diagnose the first post-critical seed.
 */
export function makeMassConservingIceSeed({
  waterMass,
  inventory,
  T,
  pPa = 101325,
  aw,
  pathway,
  t,
  gammaIl,
  radiusFactor = 1.05,
}) {
  const initialWaterMass = Math.max(0, waterMass);
  const embryo = postCriticalSeedFromCNT({
    T,
    pPa,
    aw,
    totalWaterMass: initialWaterMass,
    gammaIl,
    radiusFactor,
  });
  const iceMass = embryo.iceMass;
  const liquidWaterMass = initialWaterMass - iceMass;
  const thermo = phaseEnthalpyDifference({ T, pPa });
  const latentHeatFusion = Math.max(0, thermo.latentHeat);
  const latentHeatReleased = iceMass * latentHeatFusion;
  const liquidComposition = compositionFromWaterMass(liquidWaterMass, inventory);

  return {
    createdAt: t,
    pathway,
    initialWaterMass,
    totalWaterMass: initialWaterMass,
    iceMass,
    liquidWaterMass,
    iceMassFraction: initialWaterMass > 0 ? iceMass / initialWaterMass : 0,
    iceVolume: embryo.volume,
    equivalentDiameter: embryo.volume > 0 ? Math.cbrt(6 * embryo.volume / Math.PI) : 0,
    seedRadius: embryo.radius,
    rhoIce: embryo.cnt.ice.rho,
    latentHeatFusion,
    latentHeatReleased,
    liquidComposition,
    massResidual: initialWaterMass - (iceMass + liquidWaterMass),
    cnt: embryo.cnt,
    embryo,
    phaseThermo: thermo,
    status: 'postcritical_cnt_seed',
    scientificBoundary: 'Koop hazard triggers the event; CNT reference sizes the first post-critical seed. gamma_il remains an explicit source-tagged parameter.',
  };
}

/** Advance frozen mass while conserving total H2O mass. Post-seed growth law is still deferred. */
export function advanceIceFraction(seed, targetIceFraction, inventory) {
  const frac = Math.max(seed.iceMassFraction, Math.min(1, targetIceFraction));
  const total = seed.totalWaterMass;
  const iceMass = total * frac;
  const liquidWaterMass = total - iceMass;
  const deltaIceMass = iceMass - seed.iceMass;
  const iceVolume = iceMass / seed.rhoIce;
  return {
    ...seed,
    iceMass,
    liquidWaterMass,
    iceMassFraction: frac,
    iceVolume,
    equivalentDiameter: Math.cbrt(6 * iceVolume / Math.PI),
    latentHeatReleased: seed.latentHeatReleased + Math.max(0, deltaIceMass) * seed.latentHeatFusion,
    liquidComposition: compositionFromWaterMass(liquidWaterMass, inventory),
    massResidual: total - (iceMass + liquidWaterMass),
  };
}

/** Adiabatic upper-bound temperature jump if the released heat stayed in the remaining liquid. */
export function latentHeatTemperatureJump(seed, liquidCp, liquidMass = seed.liquidWaterMass) {
  const heatCapacity = Math.max(1e-30, liquidMass * liquidCp);
  return seed.latentHeatReleased / heatCapacity;
}

export function makeHexagonalSeedFrame(radius = 1) {
  const normals=[];
  for(let i=0;i<6;i++){const theta=i*Math.PI/3;normals.push({theta,x:Math.cos(theta),y:Math.sin(theta),z:0});}
  return {radius,symmetry:6,prismNormals:normals,basalNormals:[{x:0,y:0,z:1},{x:0,y:0,z:-1}],status:'orientation_frame_only'};
}
