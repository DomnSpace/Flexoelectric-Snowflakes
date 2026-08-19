import assert from 'node:assert/strict';
import {
  saturationVaporPressureWaterMK,
  saturationVaporPressureIceMK,
  surfaceTensionWaterIAPWS,
  findKohlerCriticalPoint,
  kohlerCriticalApproximation,
  maxwellMasonRadiusRate,
} from '../solver/thermo.js';
import {
  findKohlerEquilibria,
  makeActivationLandscape,
  createAmbientDriver,
  ambientFromDriver,
} from '../solver/activation.js';

const T = 298.15, Dd = 80e-9, kappaHyg = 0.30;
const c = findKohlerCriticalPoint(Dd, kappaHyg, T);
const ca = kohlerCriticalApproximation(Dd, kappaHyg, T);
assert.ok(Math.abs(c.Dcrit * 1e9 - 470.57) < 0.8, `Dcrit ${c.Dcrit*1e9}`);
assert.ok(Math.abs(c.sCritPct - 0.29836) < 0.003, `scrit ${c.sCritPct}`);
assert.ok(Math.abs(c.Dcrit - ca.Dcrit) / c.Dcrit < 0.025);
assert.ok(Math.abs(c.sCritPct - ca.sCritPct) / c.sCritPct < 0.025);

assert.ok(Math.abs(saturationVaporPressureWaterMK(273.15) - 611.21) < 0.5);
assert.ok(Math.abs(saturationVaporPressureIceMK(273.15) - 611.15) < 0.5);
assert.ok(Math.abs(surfaceTensionWaterIAPWS(298.15) - 0.07197) < 0.0003);

const below = findKohlerEquilibria(Dd, kappaHyg, T, 1.0015);
assert.equal(below.roots.length, 2);
assert.equal(below.roots.filter(r=>r.stable).length, 1);
assert.equal(below.roots.filter(r=>!r.stable).length, 1);
assert.equal(below.barrierGone, false);
const above = findKohlerEquilibria(Dd, kappaHyg, T, 1.0045);
assert.equal(above.roots.length, 0);
assert.equal(above.barrierGone, true);

const landscape = makeActivationLandscape({Dd,kappaHyg,T,Sinf:1.0015});
assert.ok(landscape.barrier && landscape.barrier.heightJ > 0);
assert.ok(landscape.barrier.heightKBT > 0);

const mm = maxwellMasonRadiusRate(0.5e-6,T,1.005,1.002,101325);
assert.ok(mm.rate > 0);
assert.ok(mm.Fd > 0 && mm.Fk > 0);

const driver = createAmbientDriver('supersat_ramp',{T,supersaturationPct:0.05},{targetSupersaturationPct:0.45,rampDuration:10});
assert.ok(Math.abs(ambientFromDriver(driver,5).supersaturationPct - 0.25) < 1e-12);
const cool = createAmbientDriver('cool_fixed_e',{T,supersaturationPct:0},{coolingRateKPerMin:6,targetTemperatureK:288.15});
const coolState = ambientFromDriver(cool,60);
assert.ok(coolState.T < T);
assert.ok(coolState.Sw > 1);

console.log('Activation laboratory tests passed');
console.log({Dcrit_nm:c.Dcrit*1e9,scrit_pct:c.sCritPct,barrier_kBT:landscape.barrier.heightKBT});
