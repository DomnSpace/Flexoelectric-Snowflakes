import assert from 'node:assert/strict';
import { createAmbientDriver, ambientFromDriver, initialEquilibriumDiameter, integrateWetDiameter } from '../solver/activation.js';
import { findKohlerCriticalPoint } from '../solver/thermo.js';

const Dd = 80e-9, kappaHyg = 0.30, pressurePa = 101325;
const ambient0 = { T: 298.15, supersaturationPct: 0.05 };
const driver = createAmbientDriver('supersat_ramp', ambient0, {
  startSupersaturationPct: 0.05,
  targetSupersaturationPct: 0.45,
  rampDuration: 12,
});
let D = initialEquilibriumDiameter(Dd, kappaHyg, ambient0.T, 1 + ambient0.supersaturationPct / 100);
let crossedS = false, crossedD = false;
let t = 0;
while (t < 18) {
  const a = ambientFromDriver(driver, t);
  const c = findKohlerCriticalPoint(Dd, kappaHyg, a.T);
  if (a.Sw >= c.Scrit) crossedS = true;
  if (D >= c.Dcrit) crossedD = true;
  const step = integrateWetDiameter({ D, Dd, kappaHyg, T: a.T, Sinf: a.Sw, pressurePa, dt: 0.02 });
  D = step.D;
  t += 0.02;
}
assert.ok(crossedS, 'ambient trajectory should cross critical supersaturation');
assert.ok(crossedD, 'particle should cross critical diameter');
assert.ok(D > 1e-6, `expected micron-scale activated droplet, got ${D}`);
console.log('Trajectory test passed', { finalDiameter_um: D * 1e6 });
