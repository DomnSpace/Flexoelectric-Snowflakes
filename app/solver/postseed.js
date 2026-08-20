import {
  CONSTANTS,
  diffusivityWaterVaporAir,
  saturationVaporPressureIceMK,
  thermalConductivityAir,
} from './thermo.js';
import { iceIhProperties } from './ice-ih-eos.js';

/**
 * Milestone 2.5: compact post-critical seed growth before it is large enough to
 * resolve on the spatial hex field. This is a spherical bridge only; it must hand
 * off to the faceted field solver once radius >= resolvedRadius.
 *
 * Coupling:
 *   rho_i v = D_v (c_inf - c_s)/r
 *   v = alpha v_kin sigma_s
 *   T_s = T_inf + rho_i L_s v r / k_air
 *
 * where c_s = c_sat,i(T_s)(1+sigma_s). The three relations are iterated to a
 * self-consistent surface temperature/supersaturation/growth velocity.
 */

export function saturationMassDensityIce(T) {
  return saturationVaporPressureIceMK(T) * CONSTANTS.MW / (CONSTANTS.R * T);
}

export function iceKineticVelocity(T, rhoIce) {
  const cSat = saturationMassDensityIce(T);
  const molecularNormalSpeed = Math.sqrt(CONSTANTS.RV * T / (2 * Math.PI));
  return cSat / rhoIce * molecularNormalSpeed;
}

export function sublimationEnthalpy(T, pPa = 101325) {
  const ice = iceIhProperties(T, Math.min(210e6, Math.max(0, pPa)));
  // Vapor enthalpy relative to condensed phase is represented by the usual
  // atmospheric approximation Lv + Lf; the ice EOS supplies Lf through h_i.
  // The exact vapor EOS can replace this in Milestone 3 without changing APIs.
  const LvApprox = 2.501e6 - 2.37e3 * (T - 273.15);
  const LfApprox = Math.max(0, 333550 + 1200 * (273.15 - T) / 40);
  return { Ls: LvApprox + LfApprox, ice, status: 'vapor enthalpy compact approximation; ice state from IAPWS Ice Ih' };
}

export function attachmentCoefficient(sigmaSurface, model = {}) {
  const mode = model.mode ?? 'constant';
  if (mode === 'nucleation') {
    const A = Math.max(0, Math.min(1, model.A ?? 1));
    const sigma0 = Math.max(0, model.sigma0 ?? 0.003);
    if (!(sigmaSurface > 0)) return 0;
    return A * Math.exp(-sigma0 / sigmaSurface);
  }
  return Math.max(0, Math.min(1, model.alpha ?? 0.1));
}

export function solveSphericalIceSurface({
  radius,
  Tinf,
  vaporPressure,
  pPa = 101325,
  attachment = { mode: 'constant', alpha: 0.1 },
  iterations = 16,
}) {
  const ice = iceIhProperties(Tinf, Math.min(pPa, 210e6));
  const rhoIce = ice.rho;
  const Dv = diffusivityWaterVaporAir(Tinf, pPa);
  const kAir = thermalConductivityAir(Tinf);
  const { Ls, status: enthalpyStatus } = sublimationEnthalpy(Tinf, pPa);
  const cInf = vaporPressure * CONSTANTS.MW / (CONSTANTS.R * Tinf);
  let Ts = Tinf;
  let v = 0;
  let sigmaSurface = 0;
  let alpha = attachmentCoefficient(0.01, attachment);
  let cSatSurface = saturationMassDensityIce(Ts);

  for (let i = 0; i < iterations; i++) {
    cSatSurface = saturationMassDensityIce(Ts);
    const vkin = iceKineticVelocity(Ts, rhoIce);
    const sigmaGuess = Math.max(-0.99, (cInf / cSatSurface) - 1);
    alpha = attachmentCoefficient(Math.max(0, sigmaGuess), attachment);
    if (alpha <= 0 || vkin <= 0) {
      v = 0;
      sigmaSurface = sigmaGuess;
      break;
    }
    // Combined diffusion + surface-kinetic resistance at the current Ts.
    const denominator = rhoIce * radius + Dv * cSatSurface / (alpha * vkin);
    v = Dv * (cInf - cSatSurface) / Math.max(1e-30, denominator);
    sigmaSurface = v / (alpha * vkin);
    const newTs = Tinf + rhoIce * Ls * v * radius / Math.max(1e-12, kAir);
    if (Math.abs(newTs - Ts) < 1e-8) { Ts = newTs; break; }
    Ts = 0.5 * (Ts + newTs);
  }

  const cSurface = saturationMassDensityIce(Ts) * (1 + sigmaSurface);
  const massFlux = rhoIce * v;
  const heatFlux = massFlux * Ls;
  return {
    radius,Tinf,Ts,v,sigmaSurface,alpha,cInf,cSurface,Dv,kAir,rhoIce,Ls,
    massFlux,heatFlux,
    surfaceTemperatureRise: Ts-Tinf,
    enthalpyStatus,
    status: 'spherical_postcritical_bridge',
  };
}

export function advanceCompactSeed({
  seed,
  Tinf,
  vaporPressure,
  pPa = 101325,
  dt,
  attachment,
  resolvedRadius = 1e-6,
}) {
  const radius0 = Math.max(seed.seedRadius ?? seed.equivalentDiameter/2 ?? 1e-9, 1e-12);
  const surface = solveSphericalIceSurface({radius:radius0,Tinf,vaporPressure,pPa,attachment});
  const dr = surface.v * Math.max(0, dt);
  const radius = Math.max(radius0, radius0 + dr);
  const iceMass = 4*Math.PI*radius**3/3 * surface.rhoIce;
  const addedIceMass = Math.max(0, iceMass - (seed.iceMass ?? 0));
  return {
    ...seed,
    seedRadius: radius,
    equivalentDiameter: 2*radius,
    iceMass,
    addedIceMass,
    vaporDepositedMass: (seed.vaporDepositedMass ?? 0) + addedIceMass,
    latentHeatFromVapor: (seed.latentHeatFromVapor ?? 0) + addedIceMass*surface.Ls,
    surface,
    spatiallyResolved: radius >= resolvedRadius,
    resolvedRadius,
    status: radius >= resolvedRadius ? 'ready_for_hex_field_handoff' : 'subgrid_compact_seed_growth',
  };
}
