import assert from 'node:assert/strict';
import {
  KOOP_ACTIVITY_WINDOW,
  freezingProbabilityFromRate,
  iceEquilibriumWaterActivity,
  liquidIceChemicalPotentialDrive,
  makePhaseStateSurface,
  nucleationRateForProbability,
  pressureAdjustedIceEquilibriumWaterActivity,
  supercooledState,
} from '../solver/supercool.js';
import {
  soluteInventoryFromDryParticle,
  waterMassFromWetDiameter,
} from '../solver/composition.js';

const near=(a,b,tol,msg)=>assert.ok(Math.abs(a-b)<=tol,`${msg}: ${a} vs ${b}`);

near(iceEquilibriumWaterActivity(273.15),1,2e-3,'ice equilibrium water activity near 0 C');

const pure=liquidIceChemicalPotentialDrive(258.15,1.0);
assert.ok(pure.awIce<1,'a_w^i below melting must be below unity');
assert.ok(pure.deltaAw>0,'pure water must sit above ice-equilibrium water activity');
assert.ok(pure.deltaMuLiquidMinusIce>0,'liquid->ice thermodynamic drive should be positive');

const awEq=iceEquilibriumWaterActivity(245.15);
near(liquidIceChemicalPotentialDrive(245.15,awEq).deltaMuLiquidMinusIce,0,1e-10,'zero low-pressure liquid-ice drive on aw_i line');

// Equal common pressure penalizes ice because ice Ih has the larger molar volume.
const shifted=pressureAdjustedIceEquilibriumWaterActivity(258.15,5e6,1000);
assert.ok(shifted.awPressure>shifted.awReference,'pressure should raise equilibrium aw for ice in the current bulk approximation');

const Dd=80e-9,inventory=soluteInventoryFromDryParticle(Dd);
const init=waterMassFromWetDiameter({D:5e-6,Dd,T:258.15,pAir:101325,inventory});
const state=supercooledState({T:258.15,Dd,waterMass:init.waterMass,pAir:101325,inventory,ambientSw:1.0});
assert.equal(state.belowMelting,true);
assert.equal(state.metastableLiquid,true);
assert.ok(state.ambientSi>1,'water-saturated air below freezing should be supersaturated wrt ice');
assert.ok(state.laplacePressure>0,'Laplace pressure wired into state');
assert.ok(state.composition.aw>0&&state.composition.aw<=1,'composition-derived water activity');
assert.ok(Number.isFinite(state.effectiveKappa),'diagnostic effective kappa');

const J=1e14,D=10e-6,dt=1;
const p=freezingProbabilityFromRate(J,D,dt);
assert.ok(p>0&&p<1);
near(nucleationRateForProbability(p,D,dt),J,J*1e-10,'hazard inversion');

const surface=makePhaseStateSurface({nT:12,nAw:11});
assert.equal(surface.cells.length,132);
assert.ok(surface.cells.some(c=>c.metastable));
assert.ok(KOOP_ACTIVITY_WINDOW.min<KOOP_ACTIVITY_WINDOW.max);

console.log('supercool.test.mjs PASS');
