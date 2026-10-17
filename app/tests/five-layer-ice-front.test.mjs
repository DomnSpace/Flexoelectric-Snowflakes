import assert from 'node:assert/strict';
import {
  createFiveLayerCrystal,fiveLayerIceMassKg,fiveLayerOccupancy,
  fiveLayerDiagnostics,stepFiveLayerCrystal,settleFiveLayerFront
} from '../solver/layered-ice-front.js';
import {layeredGasIndex,axialHexRadius} from '../solver/layered-gas-fv.js';
import {layeredImplicitGasMassKg,layeredImplicitTrackedWaterKg} from '../solver/layered-vapor-implicit.js';
const close=(a,b,relative=2e-8,absolute=1e-25)=>
  assert.ok(Math.abs(a-b)<=Math.max(absolute,Math.max(Math.abs(a),Math.abs(b))*relative),
    'expected '+a+' ~= '+b);
const opts={N:11,dx:.75e-6,dz:.75e-6,seedRadiusCells:0,Tinf:258.15,
  sigmaInf:.14,prismAttachmentMPerS:5e-3,basalAttachmentMPerS:6e-3};
const state=createFiveLayerCrystal(opts);
assert.equal(state.nz,9);
assert.equal(state.iceSlabs,5);
assert.equal(state.model,'snowflake/model-b-five-independent-slabs-v1');
assert.equal(state.heatTransportSolved,false);
assert.equal(state.interfaceMotionSolved,true);
assert.equal(state.seedIceMassKg,5*state.cellIceCapacityKg);
assert.equal(state.seedSolidMask.reduce((x,y)=>x+y,0),5);
let layers=fiveLayerOccupancy(state);
assert.equal(layers.length,5);
assert.ok(layers.every(l=>l instanceof Float64Array));
assert.notEqual(layers[0],layers[1],'ice arrays must not alias');
for(let z=0;z<5;z++){
  assert.equal(layers[z][5*11+5],1);
  assert.equal(layers[z][5*11+6],0);
}
close(fiveLayerIceMassKg(state),state.seedIceMassKg,1e-12,1e-30);
const before=layeredImplicitTrackedWaterKg(state);
let sumDeposited=0;
let last;
for(let step=0;step<8;step++){
  last=stepFiveLayerCrystal(state,.015,{tolerance:1e-8,maxIterations:1600});
  assert.ok(last.relativeResidual<1.25e-8);
  assert.ok(last.depositedKg>=0);
  assert.equal(last.independentLayerMaps,true);
  assert.equal(last.heatTransportSolved,false);
  assert.equal(last.layerDiagnostics.massKg.length,5);
  sumDeposited+=last.depositedKg;
}
close(state.totalDepositedKg,sumDeposited,2e-9,1e-30);
close(fiveLayerIceMassKg(state)-state.seedIceMassKg,
  state.totalDepositedKg,1e-7,1e-25);
close(state.latentReleasedJ,state.totalDepositedKg*state.Ls,1e-10,1e-25);
close(layeredImplicitTrackedWaterKg(state)-before,
  state.externalMassInKg-state.externalMassOutKg,5e-8,1e-25);
assert.equal(state.seedSolidMask.reduce((x,y)=>x+y,0),5);
assert.ok(state.c.every(Number.isFinite));
assert.ok(state.iceFraction.every(p=>p>=0&&p<=1));
const out=fiveLayerDiagnostics(state);
assert.ok(out.newGrowthSites.some(x=>x>0),'at least one voxel gains physical ice occupancy');
layers=fiveLayerOccupancy(state);
const flankIndex=5*state.N+6;
assert.ok(layers.some(l=>l[flankIndex]>0),'ice advances from exposed seed face');
const different=layers.some(l=>Math.abs(l[flankIndex]-layers[2][flankIndex])>1e-14);
assert.ok(different,'five layers must develop distinct growth under shared gas transport');
close(out.massKg[0],out.massKg[4],1e-6,1e-27);
close(out.massKg[1],out.massKg[3],1e-6,1e-27);
assert.ok(out.massKg[0]>out.massKg[2],
  'exposed outer basal surfaces can capture more water than the shielded middle');
assert.ok(last.layerDiagnostics.occupiedSites.every(x=>x>=1));

// No attachment = no growth despite large supersaturation.
const neutral=createFiveLayerCrystal({...opts,prismAttachmentMPerS:0,basalAttachmentMPerS:0});
const np=stepFiveLayerCrystal(neutral,.02);
assert.equal(np.depositedKg,0);
assert.deepEqual(fiveLayerOccupancy(neutral),fiveLayerOccupancy(createFiveLayerCrystal({
  ...opts,prismAttachmentMPerS:0,basalAttachmentMPerS:0
})));

// Closed reservoir: no water appears. It may move between gas and ice.
const closed=createFiveLayerCrystal({...opts,
  boundaries:{lateral:'sealed',top:'sealed',bottom:'sealed'}});
const initialClosed=layeredImplicitTrackedWaterKg(closed);
for(let k=0;k<4;k++)stepFiveLayerCrystal(closed,.01);
close(layeredImplicitTrackedWaterKg(closed),initialClosed,1e-8,1e-27);
assert.ok(closed.totalDepositedKg>0);

// Stress a near-complete independent ice voxel with externally known mass
// (avoids waiting years of physical growth for one coarse 0.75µm voxel).
const completion=createFiveLayerCrystal({...opts,seedLayerIndices:[2]});
const site=layeredGasIndex(completion,1,0,completion.gasHaloSlabs+2);
assert.equal(completion.iceFraction[site],0);
completion.iceFraction[site]=.999;
completion.iceMassAtCellKg[site]=completion.cellIceCapacityKg*.999;
const gasOld=layeredImplicitGasMassKg(completion);
const mass=.002*completion.cellIceCapacityKg;
const moved=settleFiveLayerFront(completion,[{
  fromGasIndex:site,toIceIndex:site,iceLayer:2,massKg:mass
}]);
assert.equal(moved.fullCellsAdded,1);
assert.equal(completion.solid[site],1);
assert.equal(completion.c[site],0);
close(layeredImplicitGasMassKg(completion)+completion.trappedVaporKg,
  gasOld,1e-9,1e-27);
assert.ok(completion.iceMassAtCellKg[site]>=completion.cellIceCapacityKg,
  'excess deposition remains in physical ice-mass ledger');

console.log('five-layer-ice-front.test.mjs PASS',{
  massDepositedKg:state.totalDepositedKg,slabMassKg:out.massKg,
  flankPhase:layers.map(l=>l[flankIndex]),
  layersDiffer:different,outerBasalDominatesInThisControlledCase:true,
  simulationSeconds:state.timeSeconds,gasMassKg:out.vaporInventoryKg
});
