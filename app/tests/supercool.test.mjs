import assert from 'node:assert/strict';
import {
  KOOP_ACTIVITY_WINDOW,
  freezingProbabilityFromRate,
  iceEquilibriumWaterActivity,
  liquidIceChemicalPotentialDrive,
  makePhaseStateSurface,
  nucleationRateForProbability,
  supercooledState,
} from '../solver/supercool.js';

const near=(a,b,tol,msg)=>assert.ok(Math.abs(a-b)<=tol,`${msg}: ${a} vs ${b}`);

// At the triple/melting neighborhood p_i ~= p_w, hence a_w^i ~= 1.
near(iceEquilibriumWaterActivity(273.15),1,2e-3,'ice equilibrium water activity near 0 C');

// Pure liquid water below melting must be thermodynamically metastable to ice.
const pure=liquidIceChemicalPotentialDrive(258.15,1.0);
assert.ok(pure.awIce<1,'a_w^i below melting must be below unity');
assert.ok(pure.deltaAw>0,'pure water must sit above ice-equilibrium water activity');
assert.ok(pure.deltaMuLiquidMinusIce>0,'liquid->ice thermodynamic drive should be positive');

// A solution adjusted exactly to the ice-equilibrium water activity has zero bulk drive.
const awEq=iceEquilibriumWaterActivity(245.15);
near(liquidIceChemicalPotentialDrive(245.15,awEq).deltaMuLiquidMinusIce,0,1e-10,'zero liquid-ice drive on aw_i line');

const state=supercooledState({T:258.15,D:5e-6,Dd:80e-9,kappaHyg:0.3,ambientSw:1.0});
assert.equal(state.belowMelting,true);
assert.equal(state.metastableLiquid,true);
assert.ok(state.ambientSi>1,'water-saturated air below freezing should be supersaturated wrt ice');

// Exact Poisson conversion between volumetric rate and one-particle probability.
const J=1e14,D=10e-6,dt=1;
const p=freezingProbabilityFromRate(J,D,dt);
assert.ok(p>0&&p<1);
near(nucleationRateForProbability(p,D,dt),J,J*1e-10,'hazard inversion');

const surface=makePhaseStateSurface({nT:12,nAw:11});
assert.equal(surface.cells.length,132);
assert.ok(surface.cells.some(c=>c.metastable));
assert.ok(KOOP_ACTIVITY_WINDOW.min<KOOP_ACTIVITY_WINDOW.max);

console.log('supercool boundary tests: ok');
