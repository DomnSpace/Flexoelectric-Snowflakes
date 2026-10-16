import assert from 'node:assert/strict';
import {
  createLayeredGasFV,layeredGasMassKg,layeredTrackedWaterKg,
  assembleLayeredImplicit,solveLayeredImplicit,stepLayeredImplicit
} from '../solver/layered-vapor-implicit.js';
import {layeredGasIndex,layeredExplicitDtLimit} from '../solver/layered-gas-fv.js';
import {saturationMassDensityIce} from '../solver/postseed.js';

const b={N:11,iceSlabs:5,gasHaloSlabs:2,dx:.75e-6,dz:.75e-6,
  Tinf:258.15,pressurePa:101325,sigmaInf:.10};
function near(a,c,rel=3e-9,abs=1e-28){
  assert.ok(Number.isFinite(a)&&Math.abs(a-c)<Math.max(abs,rel*Math.abs(c)),
    'expected '+a+' ~= '+c);
}
const closed=createLayeredGasFV(b);
const cell=layeredGasIndex(closed,0,0,4),up=layeredGasIndex(closed,0,0,5);
closed.c[cell]*=1.20;
const initial=layeredTrackedWaterKg(closed);
const g=layeredExplicitDtLimit(closed);
const long=stepLayeredImplicit(closed,g*100000,{tolerance:1e-11});
assert.ok(long.iterations>0&&long.relativeResidual<2e-11);
assert.ok(closed.c[cell]<closed.c[up]*1.2);
near(layeredTrackedWaterKg(closed),initial);
assert.equal(closed.iceDepositedKg.reduce((a,v)=>a+v),0);
assert.equal(long.heatTransportSolved,false);

const ceq=saturationMassDensityIce(b.Tinf);
const ice=createLayeredGasFV({...b,
  iceCells:[{q:0,r:0,z:4}],
  basalAttachmentMPerS:3e-3,prismAttachmentMPerS:3e-3
});
const waterBefore=layeredTrackedWaterKg(ice);
const first=stepLayeredImplicit(ice,.05,{tolerance:1e-10});
assert.ok(first.depositedKg>0,'surface must receive physically depleted vapor');
assert.ok(first.relativeResidual<=1.25e-10);
assert.equal(first.depositionEvents.length,8);
near(layeredTrackedWaterKg(ice),waterBefore,2e-9,1e-28);
near(ice.totalDepositedKg,first.perIceSlabDepositedKg.reduce((a,b)=>a+b),1e-10);
near(ice.latentReleasedJ,ice.totalDepositedKg*ice.Ls,1e-10,1e-26);
assert.ok(first.dtSeconds/layeredExplicitDtLimit(ice)>100000,
  'backward-Euler must step far beyond explicit nanosecond CFL');

const compete=createLayeredGasFV({...b,
  iceCells:[{q:0,r:0,z:3},{q:0,r:0,z:5}],
  basalAttachmentMPerS:5e-3,prismAttachmentMPerS:0
});
const before=layeredTrackedWaterKg(compete);
const shared=stepLayeredImplicit(compete,.03,{tolerance:1e-10});
assert.ok(shared.perIceSlabDepositedKg[1]>0);
near(shared.perIceSlabDepositedKg[1],shared.perIceSlabDepositedKg[3],3e-8);
near(layeredTrackedWaterKg(compete),before,3e-9);
const center=layeredGasIndex(compete,0,0,4);
assert.ok(compete.c[center]<ceq*1.10);
assert.ok(compete.c[center]>=ceq);

const sheltered=createLayeredGasFV({...b,
  iceCells:[...Array.from({length:3},(_,l)=>l+3).flatMap(z=>
    [[0,0],[1,0],[0,1],[1,-1],[-1,0],[0,-1],[-1,1]]
      .map(([q,r])=>({q,r,z})))],
  prismAttachmentMPerS:5e-3,basalAttachmentMPerS:5e-3
});
stepLayeredImplicit(sheltered,.02,{tolerance:1e-10});
assert.equal(sheltered.iceDepositedKg[layeredGasIndex(sheltered,0,0,4)],0,
  'buried ice cell has no exposed surface and no uptake');

const open=createLayeredGasFV({...b,
  boundaries:{lateral:'sealed',bottom:'sealed',top:'reservoir',
    topConcentrationKgM3:ceq*1.16},
  initialVaporKgM3:ceq*1.08
});
const firstInventory=layeredTrackedWaterKg(open);
const boundary=stepLayeredImplicit(open,.02);
assert.ok(boundary.boundaryMassInKg>0);
near(layeredTrackedWaterKg(open)-firstInventory,
  boundary.boundaryMassInKg-boundary.boundaryMassOutKg,3e-9);

const bad=createLayeredGasFV({...b,iceCells:[{q:0,r:0,z:4}],
  basalAttachmentMPerS:1e-3});
const beforeBad=bad.c.slice(),tm=bad.timeSeconds;
assert.throws(()=>stepLayeredImplicit(bad,.5,{maxIterations:1,tolerance:1e-14}),
  /did not converge/);
assert.deepEqual(bad.c,beforeBad,'failed implicit convergence must not commit vapor');
assert.equal(bad.timeSeconds,tm);
assert.equal(bad.totalDepositedKg,0);

const max=createLayeredGasFV({...b,N:25,
  iceCells:[{q:0,r:0,z:4}],basalAttachmentMPerS:2e-3});
const wall=performance.now();
const run=stepLayeredImplicit(max,.05,{tolerance:1e-9,maxIterations:1600});
const runtimeMs=performance.now()-wall;
assert.ok(run.relativeResidual<2e-9);
assert.ok(runtimeMs<15000,'small 9-slab PCG should finish promptly');
console.log('layered-implicit.test.mjs PASS',{
  diffusionRatio:long.dtSeconds/g,closedMassKg:initial,
  sharedLayerMassKg:shared.perIceSlabDepositedKg,
  singleSurfaceMassKg:first.depositedKg,iterations:run.iterations,
  runtimeMs:Number(runtimeMs.toFixed(1)),massResidualKg:first.massResidualKg
});
