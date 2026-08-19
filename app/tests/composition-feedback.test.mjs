import assert from 'node:assert/strict';
import { saturationVaporPressureWaterMK } from '../solver/thermo.js';
import {
  integrateWaterMassStep,
  soluteInventoryFromDryParticle,
  solutionDropletState,
  waterMassFromWetDiameter,
} from '../solver/composition.js';

const Dd = 80e-9;
const T0 = 273.15;
const pAir = 101325;
const inventory = soluteInventoryFromDryParticle(Dd);
const init = waterMassFromWetDiameter({ D: 5e-6, Dd, T: T0, pAir, inventory });
const initial = solutionDropletState({ waterMass: init.waterMass, Dd, T: T0, pAir, inventory });
assert.ok(Math.abs(initial.D - 5e-6) / 5e-6 < 1e-7, 'mass↔diameter inversion');
assert.ok(initial.composition.aw > 0 && initial.composition.aw <= 1, 'composition water activity');
assert.ok(initial.laplacePressure > 0, 'positive Laplace pressure');
assert.ok(initial.waterProperties.rho > 900, 'IAPWS liquid density wired');

const coldMassConserved = solutionDropletState({ waterMass: init.waterMass, Dd, T: 253.15, pAir, inventory });
assert.ok(Math.abs(coldMassConserved.D - initial.D) > 1e-12, 'cooling density feedback changes size at fixed water mass');

let mass = init.waterMass;
let T = T0;
const e0 = saturationVaporPressureWaterMK(T0);
for (let i = 0; i < 1200; i++) {
  T = Math.max(248.15, T - (6 / 60) * 0.05);
  const stepped = integrateWaterMassStep({ waterMass: mass, Dd, T, pAir, vaporPressure: e0, inventory, dt: 0.05 });
  mass = stepped.waterMass;
}
const final = solutionDropletState({ waterMass: mass, Dd, T, pAir, inventory });
assert.ok(mass > init.waterMass, 'cooling at fixed vapor pressure condenses water');
assert.ok(final.D > initial.D, 'condensation grows the droplet');
assert.ok(final.composition.aw >= initial.composition.aw, 'dilution raises water activity for fixed solute inventory');

console.log('composition-feedback.test.mjs PASS');
