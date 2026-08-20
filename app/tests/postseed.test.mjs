import assert from 'node:assert/strict';
import { saturationVaporPressureIceMK } from '../solver/thermo.js';
import { solveSphericalIceSurface, advanceCompactSeed } from '../solver/postseed.js';

const T=258.15,p=101325,r=50e-9;
const e=saturationVaporPressureIceMK(T)*1.08;
const surf=solveSphericalIceSurface({radius:r,Tinf:T,vaporPressure:e,pPa:p,attachment:{mode:'constant',alpha:0.2}});
assert.ok(surf.v>0,'ice seed should grow in ice-supersaturated vapor');
assert.ok(surf.Ts>=T,'latent heat should not cool the interface');
assert.ok(surf.sigmaSurface>=0,'surface supersaturation should be nonnegative during deposition');
assert.ok(surf.massFlux>0&&surf.heatFlux>0);

let seed={seedRadius:r,equivalentDiameter:2*r,iceMass:4*Math.PI*r**3/3*surf.rhoIce,vaporDepositedMass:0,latentHeatFromVapor:0};
const next=advanceCompactSeed({seed,Tinf:T,vaporPressure:e,pPa:p,dt:0.1,attachment:{mode:'constant',alpha:0.2},resolvedRadius:1e-6});
assert.ok(next.seedRadius>seed.seedRadius);
assert.ok(next.iceMass>seed.iceMass);
assert.ok(next.vaporDepositedMass>0);
assert.equal(next.spatiallyResolved,false);

console.log('postseed.test.mjs PASS');
