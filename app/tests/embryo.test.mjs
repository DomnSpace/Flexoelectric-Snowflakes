import assert from 'node:assert/strict';
import { cntCriticalEmbryo, phaseEnthalpyDifference, postCriticalSeedFromCNT } from '../solver/embryo.js';

const T=235, p=101325, aw=0.98;
const cnt=cntCriticalEmbryo({T,pPa:p,aw});
assert.ok(cnt.iceFavored);
assert.ok(Number.isFinite(cnt.criticalRadius) && cnt.criticalRadius>0);
assert.ok(Number.isFinite(cnt.barrierJ) && cnt.barrierJ>0);
assert.ok(Number.isFinite(cnt.criticalMass) && cnt.criticalMass>0);
assert.ok(Number.isFinite(cnt.criticalMolecules) && cnt.criticalMolecules>0);

const post=postCriticalSeedFromCNT({T,pPa:p,aw,totalWaterMass:1e-13,radiusFactor:1.05});
assert.ok(post.iceMass>0);
assert.ok(post.radius>=cnt.criticalRadius || post.cappedByAvailableWater);
assert.ok(post.iceMass<=1e-13);

const latent=phaseEnthalpyDifference({T,pPa:p});
assert.ok(latent.latentHeat>0);
assert.ok(Number.isFinite(latent.liquid.h));
assert.ok(Number.isFinite(latent.ice.h));

console.log('CNT embryo tests passed');
