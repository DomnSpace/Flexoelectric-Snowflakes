import assert from 'node:assert/strict';
import {
  saturationVaporPressureWaterMK,
  saturationVaporPressureIceMK,
  surfaceTensionWaterIAPWS,
  findKohlerCriticalPoint,
  kohlerCriticalApproximation,
  kohlerSaturationRatio,
  activationBudget,
} from '../solver/thermo.js';

const near = (a,b,rel=5e-4) => Math.abs(a-b) <= rel*Math.max(Math.abs(a),Math.abs(b),1);
assert.ok(near(saturationVaporPressureWaterMK(273.15), 611.21, 3e-3));
assert.ok(near(saturationVaporPressureIceMK(273.15), 611.15, 3e-3));
assert.ok(near(surfaceTensionWaterIAPWS(298.15), 0.07197, 5e-3));

const Dd=80e-9, kappa=0.3, T=298.15;
const crit=findKohlerCriticalPoint(Dd,kappa,T);
const approx=kohlerCriticalApproximation(Dd,kappa,T);
assert.ok(crit.Dcrit > Dd);
assert.ok(crit.Scrit > 1);
assert.ok(Math.abs(crit.Dcrit-approx.Dcrit)/crit.Dcrit < 0.08);
assert.ok(Math.abs(crit.sCritPct-approx.sCritPct)/crit.sCritPct < 0.08);

const below=kohlerSaturationRatio(Dd*1.01,Dd,kappa,T);
assert.ok(below < 1);
const budget=activationBudget(crit.Dcrit,Dd,kappa,T,crit.Scrit);
assert.ok(Math.abs(budget.net) < 5e-8);
console.log('thermo tests passed', {Dcrit_nm:crit.Dcrit*1e9, sCrit_pct:crit.sCritPct});
