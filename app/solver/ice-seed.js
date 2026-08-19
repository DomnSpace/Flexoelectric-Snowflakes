import { compositionFromWaterMass } from './composition.js';
import { postCriticalSeedFromCNT, phaseEnthalpyDifference } from './embryo.js';

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
  const initialComposition = compositionFromWaterMass(initialWaterMass, inventory);
  const awResolved = Number.isFinite(aw) ? aw : initialComposition.aw;
  const embryo = postCriticalSeedFromCNT({
    T,
    pPa,
    aw: awResolved,
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
    createdAt:t,pathway,initialWaterMass,totalWaterMass:initialWaterMass,
    initialComposition,awAtNucleation:awResolved,
    iceMass,liquidWaterMass,iceMassFraction:initialWaterMass>0?iceMass/initialWaterMass:0,
    iceVolume:embryo.volume,equivalentDiameter:embryo.volume>0?Math.cbrt(6*embryo.volume/Math.PI):0,
    seedRadius:embryo.radius,rhoIce:embryo.cnt.ice.rho,
    latentHeatFusion,latentHeatReleased,liquidComposition,
    massResidual:initialWaterMass-(iceMass+liquidWaterMass),cnt:embryo.cnt,embryo,phaseThermo:thermo,
    status:'postcritical_cnt_seed',
    scientificBoundary:'Koop hazard triggers the event; CNT reference sizes the first post-critical seed. gamma_il remains an explicit source-tagged parameter.',
  };
}

export function advanceIceFraction(seed,targetIceFraction,inventory){const frac=Math.max(seed.iceMassFraction,Math.min(1,targetIceFraction)),total=seed.totalWaterMass,iceMass=total*frac,liquidWaterMass=total-iceMass,deltaIceMass=iceMass-seed.iceMass,iceVolume=iceMass/seed.rhoIce;return{...seed,iceMass,liquidWaterMass,iceMassFraction:frac,iceVolume,equivalentDiameter:Math.cbrt(6*iceVolume/Math.PI),latentHeatReleased:seed.latentHeatReleased+Math.max(0,deltaIceMass)*seed.latentHeatFusion,liquidComposition:compositionFromWaterMass(liquidWaterMass,inventory),massResidual:total-(iceMass+liquidWaterMass)}}
export function latentHeatTemperatureJump(seed,liquidCp,liquidMass=seed.liquidWaterMass){const heatCapacity=Math.max(1e-30,liquidMass*liquidCp);return seed.latentHeatReleased/heatCapacity}
export function makeHexagonalSeedFrame(radius=1){const normals=[];for(let i=0;i<6;i++){const theta=i*Math.PI/3;normals.push({theta,x:Math.cos(theta),y:Math.sin(theta),z:0})}return{radius,symmetry:6,prismNormals:normals,basalNormals:[{x:0,y:0,z:1},{x:0,y:0,z:-1}],status:'orientation_frame_only'}}
