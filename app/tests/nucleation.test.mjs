import assert from 'node:assert/strict';
import {
  koopHomogeneousRate,
  freezingProbabilityFromHazard,
  makeNucleationClock,
  advanceHazardClock,
} from '../solver/nucleation.js';
import { makeMassConservingIceSeed } from '../solver/ice-seed.js';
import { soluteInventoryFromDryParticle } from '../solver/composition.js';

for (const d of [0.27,0.30,0.33]) {
  const r=koopHomogeneousRate(d);
  assert.equal(r.valid,true);
  assert.ok(Number.isFinite(r.log10Jcm3s));
  assert.ok(r.Jm3s>0);
}
assert.equal(koopHomogeneousRate(0.20).valid,false);
assert.equal(koopHomogeneousRate(0.40).valid,false);

const H=2.3;
assert.ok(Math.abs(freezingProbabilityFromHazard(H)-(1-Math.exp(-H)))<1e-15);

const clock=makeNucleationClock(42);
const threshold=clock.threshold;
advanceHazardClock(clock,threshold/2,1,1,'homogeneous');
assert.equal(clock.crossed,false);
advanceHazardClock(clock,threshold/2+1e-12,1,2,'homogeneous');
assert.equal(clock.crossed,true);
assert.equal(clock.pathway,'homogeneous');

const inv=soluteInventoryFromDryParticle(80e-9);
const seed=makeMassConservingIceSeed({waterMass:1e-13,inventory:inv,T:235,pathway:'homogeneous',t:12});
assert.ok(seed.iceMass>0);
assert.ok(seed.liquidWaterMass>0);
assert.ok(Math.abs(seed.massResidual)<1e-28);
assert.ok(seed.latentHeatReleased>0);
assert.equal(seed.pathway,'homogeneous');

console.log('Milestone 2 nucleation tests passed');
