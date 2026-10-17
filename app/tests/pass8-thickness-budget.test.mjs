import assert from 'node:assert/strict';
import {makeHexGrowthState,resolvedHexGrowthStep,enlargeHexDomain,
  auditHexDomain} from '../solver/hex25d.js';
import {iceIhProperties} from '../solver/ice-ih-eos.js';
import {localIceMass,localHeightStats,localCellArea,seedLocalHeight,
  acceptPrismMass,acceptBasalMass,finishLocalHeightStep} from '../solver/ice-thickness.js';
import {versionLocalIceSnapshot,verifyLocalIceSnapshot,
  localThicknessAt} from '../solver/ice-snapshot-v2.js';
import {makeUniformLayerPlan} from '../render/uniform-five-layers.js';

const temperature=258.15,pressure=101325,dx=.75e-6;
const rho=iceIhProperties(temperature,pressure).rho;
const nearly=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<=
  Math.max(1e-24,Math.abs(a),Math.abs(b))*tol,`expected ${a} ~= ${b}`);
const opts={N:73,dx,Tinf:temperature,pressurePa:pressure,
  sigmaInf:.09,seedRadius:1.5e-6,localHeight:true,
  kinetics:{prism:{mode:'nucleation',A:1,sigma0:.0035},
    basal:{mode:'nucleation',A:.01,sigma0:.01},
    anisotropy6:.35,kineticPhaseDeg:30,surfaceEnergyEpsilon:.01,
    surfacePhaseDeg:30,sdak:false,capillarity:true,gammaIv:.109}};
const local=makeHexGrowthState(opts);
assert.equal(local.localHeight,true);
assert.equal(local.topHeight.length,73*73);
assert.equal(local.heightLedger.initialMassKg,localIceMass(local,rho));
assert.ok(local.heightLedger.initialMassKg>0);
assert.equal(local.heightLedger.basalVerticalInKg,0);
assert.equal(local.totalDepositedMass,0);
assert.ok(localHeightStats(local,rho).maxThicknessM>0);

// Demonstrate geometric, rather than cosmetic, mass admission in a controlled
// cell: prism growth is capped by available φ and height is inherited.
const test=makeHexGrowthState({...opts,N:41,seedRadius:1.5e-6});
const m=(test.N-1)/2,free=m*test.N+m+4;
assert.equal(test.phi[free],0);
const minDepth=2*test.basalHalfThickness;
const capacity=rho*localCellArea(test)*minDepth;
const p=acceptPrismMass(test,free,capacity*1.5,rho,
  test.basalHalfThickness,-test.basalHalfThickness);
nearly(p,capacity);
assert.equal(test.phi[free],1);
assert.equal(test.ice[free],1);
const basal=acceptBasalMass(test,free,rho*1e-9,.5,rho);
assert.ok(basal>0);
nearly(basal,rho*localCellArea(test)*1e-9);
const before=test.heightLedger.initialMassKg;
const result=finishLocalHeightStep(test,{
  rhoIce:rho,Ls:2.8e6,dt:.5,prismKg:p,basalKg:basal
});
assert.ok(result.relativeLedgerResidual<2e-10);
nearly(result.physicalMassKg-before,p+basal);
nearly(test.totalDepositedMass,p+basal);
nearly(test.totalLatentHeat,(p+basal)*2.8e6);
assert.equal(test.heightLedger.prismInKg,p);
assert.equal(test.heightLedger.basalVerticalInKg,basal);

// Zero kinetic attachment -> strictly no growth, no thermal mass creation.
// The legacy slab must not accidentally acquire the new source accounting.
const zero=makeHexGrowthState({...opts,kinetics:{
  ...opts.kinetics,
  prism:{mode:'nucleation',A:0,sigma0:.003},
  basal:{mode:'nucleation',A:0,sigma0:.01}
}});
const baseline=makeHexGrowthState({...opts,localHeight:false});
for(let step=0;step<6;step++){
  resolvedHexGrowthStep(zero,.5,{vaporIterations:25,heatIterations:25});
  resolvedHexGrowthStep(baseline,.5,{vaporIterations:25,heatIterations:25});
}
nearly(zero.totalDepositedMass,0);
assert.equal(zero.heightLedger.prismInKg,0);
assert.equal(zero.heightLedger.basalVerticalInKg,0);
assert.equal(baseline.localHeight,undefined);
assert.equal(baseline.topHeight,undefined);
assert.ok(baseline.totalDepositedMass>0,'legacy physical solver still grows');

// Real live height-coupled growth: each admitted mass appears in actual ice
// geometry; incoming prism and basal sources are separate and heat is counted.
let oldMass=localIceMass(local,rho);
let previousTotal=0,oldLatent=0;
for(let k=0;k<20;k++){
  resolvedHexGrowthStep(local,.5,{vaporIterations:28,heatIterations:32});
  const nextMass=localIceMass(local,rho);
  const ledger=local.heightLedger,expected=ledger.prismInKg+
    ledger.basalVerticalInKg;
  assert.ok(local.totalDepositedMass>previousTotal,'real ice deposits in both directions');
  assert.ok(local.totalLatentHeat>oldLatent,'both sources create latent heat');
  assert.ok(nextMass>=oldMass,'without sublimation the ice cannot lose mass');
  nearly(nextMass-ledger.initialMassKg,expected,2e-10);
  nearly(local.totalDepositedMass,expected);
  nearly(local.totalLatentHeat,local.totalDepositedMass*local.lastBudget.Ls);
  assert.ok(local.lastBudget.basalMassRate>0);
  assert.ok(local.lastBudget.prismMassRate>0);
  assert.ok(local.lastBudget.relativeMassResidual<2e-10);
  assert.ok(local.lastBudget.measuredThicknessMaxM>0);
  assert.ok(Number.isFinite(local.lastBudget.maxDeltaT));
  oldMass=nextMass;previousTotal=local.totalDepositedMass;
  oldLatent=local.totalLatentHeat;
}
const stats=localHeightStats(local,rho);
assert.ok(stats.occupiedCells>0);
assert.ok(stats.maxThicknessM>=stats.meanThicknessM);
assert.ok(stats.minThicknessM<=stats.meanThicknessM);
assert.ok(local.lastBudget.basalMassRate>0);
assert.ok(local.lastBudget.prismMassRate>0);
assert.ok(auditHexDomain(local).haloCells>6);

// Domain growth must NOT reset the local top/bottom profiles or the ledger.
const expanded=enlargeHexDomain(local,97);
assert.equal(expanded.localHeight,true);
assert.equal(expanded.t,local.t);
assert.equal(expanded.totalDepositedMass,local.totalDepositedMass);
assert.equal(expanded.heightLedger.initialMassKg,local.heightLedger.initialMassKg);
assert.equal(expanded.heightLedger.prismInKg,local.heightLedger.prismInKg);
assert.equal(expanded.heightLedger.basalVerticalInKg,local.heightLedger.basalVerticalInKg);
nearly(localIceMass(expanded,rho),localIceMass(local,rho),5e-11);
resolvedHexGrowthStep(expanded,.5,{vaporIterations:20,heatIterations:20});
assert.ok(expanded.lastBudget.relativeMassResidual<2e-10);

// Reproducible snapshot v2: copied height fields not aliases, and five-layer
// renderer remains an explicitly identified MEAN-height guide until Pass 9.
const morph={radius:Array.from({length:180},()=>2e-6),
  halfThickness:expanded.basalHalfThickness};
const frame=versionLocalIceSnapshot({
  n:expanded.N,dx:expanded.dx,t:expanded.t,topology:'triangular',
  phi:new Float32Array(expanded.phi),ice:new Uint8Array(expanded.ice),
  u:new Float32Array(expanded.sigma),dT:new Float32Array(expanded.dT),
  mask:new Uint8Array(expanded.mask),morph
},expanded,rho);
assert.ok(verifyLocalIceSnapshot(frame));
assert.equal(frame.snapshotSchema.version,2);
assert.equal(frame.iceStructure.localHeightSolved,true);
assert.equal(frame.iceStructure.verticalVaporTransportSolved,false);
assert.equal(frame.iceStructure.basalMassCoupled,true);
nearly(frame.iceStructure.iceMassKg,localIceMass(expanded,rho));
const k=expanded.ice.findIndex(x=>x>0);
nearly(localThicknessAt(frame,k),
  expanded.topHeight[k]-expanded.bottomHeight[k],1e-7);
frame.hTop[k]=0;
assert.ok(expanded.topHeight[k]>0,'snapshot arrays are copied');
const plan=makeUniformLayerPlan(frame);
assert.equal(plan.sliceCount,5);
assert.equal(plan.mode,'mean-height-guide-not-local-slices');
assert.equal(plan.sliceGeometryFromLocalHeight,false);
assert.equal(plan.physicalLocalHeightSolved,true);
assert.throws(()=>verifyLocalIceSnapshot({...frame,snapshotSchema:
  {...frame.snapshotSchema,version:3}}),/Invalid/);
console.log('pass8-thickness-budget.test.mjs PASS',{
  steps:20,netDepositPg:local.totalDepositedMass*1e12,
  prismPg:local.heightLedger.prismInKg*1e12,
  externalBasalPg:local.heightLedger.basalVerticalInKg*1e12,
  meanThicknessUm:stats.meanThicknessM*1e6,
  maxThicknessUm:stats.maxThicknessM*1e6,
  openReservoirResidual:local.heightLedger.relativeResidual,
  expandedN:expanded.N,previewIsMeanHeightOnly:true
});
