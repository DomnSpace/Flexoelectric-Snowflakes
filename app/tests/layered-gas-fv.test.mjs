import assert from 'node:assert/strict';
import {
  FIVE_ICE_SLABS,DEFAULT_HALO_SLABS,phaseCellArea,axialHexRadius,
  createLayeredGasFV,layeredGasIndex,layeredGasMassKg,
  layeredTrackedWaterKg,layeredExplicitDtLimit,stepLayeredGasFV
} from '../solver/layered-gas-fv.js';
import {saturationMassDensityIce} from '../solver/postseed.js';

const close=(got,want,atol=1e-27,rtol=1e-10)=>{
  assert.ok(Number.isFinite(got)&&Math.abs(got-want)<=atol+rtol*Math.max(Math.abs(got),Math.abs(want)),
    `expected ${got} approximately ${want}`);
};
const base={N:9,dx:.75e-6,dz:.6e-6,iceSlabs:5,gasHaloSlabs:2,
  Tinf:258.15,pressurePa:101325,sigmaInf:.10};
const s=createLayeredGasFV(base);
assert.equal(s.nz,9);
assert.equal(s.iceSlabs,FIVE_ICE_SLABS);
assert.equal(s.gasHaloSlabs,DEFAULT_HALO_SLABS);
assert.equal(s.basalSites,1+3*4*5);
assert.equal(s.hexMask.reduce((a,b)=>a+b),61);
assert.equal(layeredGasIndex(s,0,0,4),4*s.nxy+4*s.N+4);
assert.equal(layeredGasIndex(s,5,0,4),-1);
assert.equal(axialHexRadius(2,-3),3);
close(s.area,Math.sqrt(3)/2*s.dx*s.dx,1e-27);
assert.ok(layeredExplicitDtLimit(s)>0);
assert.equal(s.heatTransportSolved,false);
assert.equal(s.interfaceMotionSolved,false);
const separateBoundary=createLayeredGasFV({...base,surfaceEquilibriumKgM3:0});
assert.ok(separateBoundary.c[layeredGasIndex(separateBoundary,0,0,4)]>0,
  'changing ice-surface equilibrium must not erase the independent atmospheric vapor inventory');
assert.equal(separateBoundary.ceq,0);

assert.throws(()=>createLayeredGasFV({...base,N:8}),/odd basal/);
assert.throws(()=>createLayeredGasFV({...base,iceCells:[{q:0,r:0,z:0}]}),/eligible/);

// Closed, *nonuniform* 3D triangular vapor diffusion. Every gas-to-gas face
// contributes the same signed mass transfer to both volumes.
const k0=layeredGasIndex(s,0,0,4);
const kTop=layeredGasIndex(s,0,0,5),kRight=layeredGasIndex(s,1,0,4);
s.c[k0]*=1.35;
const initialGas=layeredGasMassKg(s);
const cStart=s.c[k0],cTop=s.c[kTop],cRight=s.c[kRight];
for(let i=0;i<40;i++){
  const rep=stepLayeredGasFV(s,layeredExplicitDtLimit(s)*.75);
  assert.equal(rep.limitedCells,0);
  assert.equal(rep.depositedKg,0);
  assert.equal(rep.boundaryMassInKg,0);
  assert.equal(rep.boundaryMassOutKg,0);
}
close(layeredGasMassKg(s),initialGas,1e-27,2e-12);
assert.ok(s.c[k0]<cStart,'center must transfer vapor to neighbors');
assert.ok(s.c[kTop]>cTop,'vertical vapor neighbor must receive mass');
assert.ok(s.c[kRight]>cRight,'triangular lateral neighbor must receive mass');
assert.ok(s.c.every(x=>x>=0&&Number.isFinite(x)));

// A single ice surface is a genuine sink of the very same gas inventory.
// ICE GEOMETRY remains fixed, while deposited mass is ledgered independently.
const ice=createLayeredGasFV({...base,
  iceCells:[{q:0,r:0,z:4}],
  prismAttachmentMPerS:2e-3,basalAttachmentMPerS:2e-3
});
const totalBefore=layeredTrackedWaterKg(ice);
const gasBefore=layeredGasMassKg(ice);
const rep=stepLayeredGasFV(ice,layeredExplicitDtLimit(ice)*.7);
assert.ok(rep.depositedKg>0);
assert.equal(rep.limitedCells,0);
assert.equal(rep.boundaryMassInKg,0);
assert.equal(rep.boundaryMassOutKg,0);
close(gasBefore-layeredGasMassKg(ice),rep.depositedKg,1e-26,4e-12);
close(layeredTrackedWaterKg(ice),totalBefore,1e-26,4e-12);
close(ice.iceDepositionPerSlabKg[4],rep.depositedKg,1e-30);
close(rep.latentReleasedJ,rep.depositedKg*ice.Ls,1e-22);
close(ice.latentReleasedJ,rep.latentReleasedJ,1e-22);
assert.ok(ice.solid[layeredGasIndex(ice,0,0,4)]===1,'solid support is not fabricated');
assert.equal(ice.c[layeredGasIndex(ice,0,0,4)],0,'solid cells have no vapor inventory');
assert.equal(rep.heatTransportSolved,false);

// Buried ice cells may have deposited water on their external PRISM faces,
// but their completely closed vertical faces must not absorb hidden gas.
// Test a 7-cell central hex block at 3 consecutive z layers: the core
// z=4 is surrounded laterally and vertically by other solid cells.
const block=[];for(const z of [3,4,5])for(const [q,r] of
  [[0,0],[1,0],[0,1],[1,-1],[-1,0],[0,-1],[-1,1]])block.push({q,r,z});
const shielded=createLayeredGasFV({...base,iceCells:block,
  prismAttachmentMPerS:1e-3,basalAttachmentMPerS:1e-3});
const b=layeredTrackedWaterKg(shielded);
stepLayeredGasFV(shielded,layeredExplicitDtLimit(shielded)*.8);
assert.equal(shielded.iceDepositedKg[layeredGasIndex(shielded,0,0,4)],0,
  'buried core receives no unaccounted deposition');
assert.ok(shielded.totalDepositedKg>0);
close(layeredTrackedWaterKg(shielded),b,1e-26,5e-12);

// Upper/lower symmetry check: symmetric vapor and one planar ice element,
// with attachment ONLY on basal-facing interfaces.
const symmetric=createLayeredGasFV({...base,iceCells:[{q:0,r:0,z:4}],
  prismAttachmentMPerS:0,basalAttachmentMPerS:4e-3});
stepLayeredGasFV(symmetric,layeredExplicitDtLimit(symmetric)*.8);
close(symmetric.c[layeredGasIndex(symmetric,0,0,3)],
  symmetric.c[layeredGasIndex(symmetric,0,0,5)],1e-17);
assert.ok(symmetric.totalDepositedKg>0);
assert.equal(symmetric.iceDepositionPerSlabKg[3],0);
assert.equal(symmetric.iceDepositionPerSlabKg[5],0);

// Two distinct solid layers request water from the SAME middle gas volume.
// Their accesses are assembled simultaneously and can never spend water twice.
const compete=createLayeredGasFV({...base,
  iceCells:[{q:0,r:0,z:3},{q:0,r:0,z:5}],
  prismAttachmentMPerS:0,basalAttachmentMPerS:2e-3});
const preShared=compete.c[layeredGasIndex(compete,0,0,4)];
const preCompete=layeredTrackedWaterKg(compete);
const shared=stepLayeredGasFV(compete,layeredExplicitDtLimit(compete)*.75);
assert.ok(shared.perIceSlabDepositedKg[1]>0);
assert.ok(shared.perIceSlabDepositedKg[3]>0);
close(shared.perIceSlabDepositedKg.reduce((a,v)=>a+v,0),shared.depositedKg,1e-28);
assert.ok(compete.c[layeredGasIndex(compete,0,0,4)]<preShared);
close(layeredTrackedWaterKg(compete),preCompete,1e-26,4e-12);

// Explicitly test oversubscription by extreme artificial sticking.
const limitCase=createLayeredGasFV({...base,
  iceCells:[{q:0,r:0,z:3},{q:0,r:0,z:5}],
  basalAttachmentMPerS:1e9,prismAttachmentMPerS:0,
  surfaceEquilibriumKgM3:0,
  initialVaporKgM3:saturationMassDensityIce(258.15)*1.1});
const cCopy=limitCase.c.slice(),beforeLimited=layeredTrackedWaterKg(limitCase);
assert.throws(()=>stepLayeredGasFV(limitCase,1e-7),/explicitly stable/);
assert.deepEqual(limitCase.c,cCopy,'refused step must not alter any physical state');
const limited=stepLayeredGasFV(limitCase,1e-7,{allowFluxLimiter:true});
assert.ok(limited.limitedCells>0);
assert.ok(limited.limiterMinimum<1);
assert.ok(limitCase.c.every(x=>x>=0&&Number.isFinite(x)));
close(layeredTrackedWaterKg(limitCase),beforeLimited,1e-25,5e-12);

// Open top/bottom reservoirs: track the exact signed incoming/outgoing flux.
const refCon=saturationMassDensityIce(258.15);
for(const sign of [-1,1]){
  const open=createLayeredGasFV({...base,
    initialVaporKgM3:refCon,
    boundaries:{lateral:'sealed',bottom:'sealed',top:'reservoir',
      topConcentrationKgM3:refCon*(sign>0?1.3:.7)}});
  const prev=layeredTrackedWaterKg(open);
  const response=stepLayeredGasFV(open,layeredExplicitDtLimit(open)*.4);
  if(sign>0){assert.ok(response.boundaryMassInKg>0);assert.equal(response.boundaryMassOutKg,0);}
  else{assert.ok(response.boundaryMassOutKg>0);assert.equal(response.boundaryMassInKg,0);}
  close(layeredTrackedWaterKg(open)-prev,
    response.boundaryMassInKg-response.boundaryMassOutKg,1e-27,3e-12);
}

// Changing Δz changes both gas volume and inter-slab diffusive conductance.
const different=createLayeredGasFV({...base,dz:1.2e-6});
close(different.volume/s.volume,2,1e-15);
assert.ok(layeredExplicitDtLimit(different)>layeredExplicitDtLimit(s));
console.log('layered-gas-fv.test.mjs PASS',{
  triangleSites:s.basalSites,
  cellSlabs:s.nz,dtLimitSeconds:layeredExplicitDtLimit(s),
  closedMassKg:initialGas,singleSurfaceMassKg:rep.depositedKg,
  sharedLayerDepositsKg:shared.perIceSlabDepositedKg,
  extremeLimiter:limited.limiterMinimum
});
