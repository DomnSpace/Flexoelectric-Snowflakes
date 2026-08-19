import assert from 'node:assert/strict';
import { supercooledWaterProperties } from '../solver/supercooled-eos.js';

const verification = [
  [273.15, 0.101325e6, 999.84229, -0.683042e-4, 5.088499e-10, 4218.3002, 1402.3886],
  [235.15, 0.101325e6, 968.09999, -29.633816e-4, 11.580785e-10, 5997.5632, 1134.5855],
  [250.00, 200e6, 1090.45677, 3.267768e-4, 3.361311e-10, 3708.3902, 1668.2020],
  [200.00, 400e6, 1185.02800, 6.716009e-4, 2.567237e-10, 3338.5250, 1899.3294],
  [250.00, 400e6, 1151.71517, 4.929927e-4, 2.277029e-10, 3757.2144, 2015.8782],
];

const rel = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(b));
for (const [T, p, rho, alpha, kappaT, cp, speed] of verification) {
  const q = supercooledWaterProperties(T, p);
  assert.ok(rel(q.rho, rho) < 2e-7, `rho @ ${T} K`);
  assert.ok(Math.abs(q.alphaP - alpha) < 2e-8, `alpha_p @ ${T} K`);
  assert.ok(rel(q.kappaT, kappaT) < 2e-6, `kappa_T @ ${T} K`);
  assert.ok(rel(q.cp, cp) < 2e-5, `cp @ ${T} K`);
  assert.ok(rel(q.speedSound, speed) < 2e-5, `sound speed @ ${T} K`);
}

console.log('supercooled-eos.test.mjs PASS');
