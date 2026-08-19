import { CONSTANTS } from './thermo.js';
import { compositionFromWaterMass } from './composition.js';

export const ICE_SEED_DEFAULTS = Object.freeze({
  iceDensity: 917.0,                  // kg m^-3; temporary until IAPWS Ice-Ih EOS module
  latentHeatFusionRef: 333550,        // J kg^-1 at ~0 C; explicit compact approximation
  minSeedMassFraction: 1e-6,
  maxInitialSeedMassFraction: 0.05,
});

/**
 * Compact initial ice seed. This is NOT yet classical embryo thermodynamics:
 * gamma_il and critical-radius energetics belong to the next refinement.
 * It is the first resolved post-nucleation mass state after a stochastic event.
 */
export function makeMassConservingIceSeed({
  waterMass,
  inventory,
  T,
  pathway,
  t,
  seedMassFraction = ICE_SEED_DEFAULTS.minSeedMassFraction,
  options = {},
}) {
  const rhoIce = options.iceDensity ?? ICE_SEED_DEFAULTS.iceDensity;
  const Lf = options.latentHeatFusion ?? ICE_SEED_DEFAULTS.latentHeatFusionRef;
  const frac = Math.max(
    ICE_SEED_DEFAULTS.minSeedMassFraction,
    Math.min(options.maxInitialSeedMassFraction ?? ICE_SEED_DEFAULTS.maxInitialSeedMassFraction, seedMassFraction)
  );
  const initialWaterMass = Math.max(0, waterMass);
  const iceMass = Math.min(initialWaterMass, initialWaterMass * frac);
  const liquidWaterMass = initialWaterMass - iceMass;
  const iceVolume = iceMass / rhoIce;
  const equivalentDiameter = iceVolume > 0 ? Math.cbrt(6 * iceVolume / Math.PI) : 0;
  const latentHeatReleased = iceMass * Lf;
  const liquidComposition = compositionFromWaterMass(liquidWaterMass, inventory);

  return {
    createdAt: t,
    pathway,
    initialWaterMass,
    totalWaterMass: initialWaterMass,
    iceMass,
    liquidWaterMass,
    iceMassFraction: initialWaterMass > 0 ? iceMass / initialWaterMass : 0,
    iceVolume,
    equivalentDiameter,
    rhoIce,
    latentHeatFusion: Lf,
    latentHeatReleased,
    liquidComposition,
    massResidual: initialWaterMass - (iceMass + liquidWaterMass),
    status: 'compact_ice_seed',
    scientificBoundary: 'Seed mass conversion is exact; ice density and latent heat are compact approximations; gamma_il/critical embryo geometry not yet loaded.',
  };
}

/** Advance frozen mass while conserving total H2O mass. */
export function advanceIceFraction(seed, targetIceFraction, inventory, options = {}) {
  const frac = Math.max(seed.iceMassFraction, Math.min(1, targetIceFraction));
  const total = seed.totalWaterMass;
  const iceMass = total * frac;
  const liquidWaterMass = total - iceMass;
  const deltaIceMass = iceMass - seed.iceMass;
  const Lf = options.latentHeatFusion ?? seed.latentHeatFusion;
  const rhoIce = options.iceDensity ?? seed.rhoIce;
  const iceVolume = iceMass / rhoIce;
  return {
    ...seed,
    iceMass,
    liquidWaterMass,
    iceMassFraction: frac,
    iceVolume,
    equivalentDiameter: Math.cbrt(6 * iceVolume / Math.PI),
    latentHeatReleased: seed.latentHeatReleased + Math.max(0, deltaIceMass) * Lf,
    liquidComposition: compositionFromWaterMass(liquidWaterMass, inventory),
    massResidual: total - (iceMass + liquidWaterMass),
  };
}

/** Sensible warming estimate if latent heat is initially trapped in the remaining droplet. */
export function latentHeatTemperatureJump(seed, liquidCp, liquidMass = seed.liquidWaterMass) {
  const heatCapacity = Math.max(1e-30, liquidMass * liquidCp);
  return seed.latentHeatReleased / heatCapacity;
}

/** Useful for an eventual interface seed: six equivalent facet normals, but no arm growth yet. */
export function makeHexagonalSeedFrame(radius = 1) {
  const normals = [];
  for (let i = 0; i < 6; i++) {
    const theta = i * Math.PI / 3;
    normals.push({ theta, x: Math.cos(theta), y: Math.sin(theta), z: 0 });
  }
  return {
    radius,
    symmetry: 6,
    prismNormals: normals,
    basalNormals: [{x:0,y:0,z:1},{x:0,y:0,z:-1}],
    status: 'orientation_frame_only',
  };
}
