import assert from 'node:assert/strict';
import {makeHexGrowthState,resolvedHexGrowthStep,hexGrowthMorphology,
  sixfoldAttachmentMultiplier,sixfoldSurfaceStiffness,pointInsideHexSeed,
  sampleHexLinear,auditHexDomain} from '../solver/hex25d.js';

const near=(x,y,tol=1e-12)=>assert.ok(Math.abs(x-y)<=tol,`expected ${x} ~= ${y}`);

// θ is the LOCAL OUTWARD interface NORMAL, not the global polar arm angle.
// Shifting it by π/6 reverses the 6-fold extrema, without any change to vapor PDE.
const a=.5,theta=.13;
near(sixfoldAttachmentMultiplier(theta,a,30),
  2-sixfoldAttachmentMultiplier(theta,a,0));
near(sixfoldAttachmentMultiplier(theta,0,0),1);
near(sixfoldAttachmentMultiplier(theta,0,30),1);
for(const x of [0,15,30]){
  const m=sixfoldAttachmentMultiplier(theta,a,x);
  assert.ok(m>=0&&m<=1+a);
}

// γ(θ)=γ0(1+ε cos6θ) => γ+γ''=γ0(1-35ε cos6θ).
near(sixfoldSurfaceStiffness(0,.02,0),.30);
near(sixfoldSurfaceStiffness(0,.02,30),1.70);
near(sixfoldSurfaceStiffness(theta,0,0),1);
assert.throws(()=>sixfoldSurfaceStiffness(0,.03,0),RangeError);
assert.throws(()=>sixfoldSurfaceStiffness(0,-.03,0),RangeError);
assert.throws(()=>sixfoldSurfaceStiffness(0,NaN,0),RangeError);

// A faceted seed is actually put into PHYSICAL occupied cells, not a rendered
// fake outline. Rotating 30 degrees changes which face points along X.
const dx=.75e-6,r=4.5e-6;
assert.equal(pointInsideHexSeed(r,0,r,0),true);
assert.equal(pointInsideHexSeed(r,0,r,30),false);
assert.equal(pointInsideHexSeed(0,0,r,30),true);
const setup=seed=>makeHexGrowthState({
  N:97,dx,seedRadius:r,seedShape:seed,seedRotationDeg:30,
  Tinf:258.15,sigmaInf:.12
});
const circle=setup('circle'),hex=setup('hex');
const nIce=s=>s.ice.reduce((x,y)=>x+y,0);
assert.notEqual(nIce(circle),nIce(hex),'the shape of the actual ice seed must change');
assert.ok(nIce(hex)<nIce(circle),'same circumradius has smaller hexagonal area');

// Anisotropy phase cannot matter with kinetic amplitude=0 and isotropic
// surface energy. Same seed, vapor field, timesteps and budget exactly.
const kin=(amplitude,phase,epsilon=0)=>({
  prism:{mode:'nucleation',A:1,sigma0:.0035},
  basal:{mode:'nucleation',A:.004,sigma0:.01},
  anisotropy6:amplitude,kineticPhaseDeg:phase,
  surfaceEnergyEpsilon:epsilon,surfacePhaseDeg:phase,
  sdak:true,sdakStrength:.55,capillarity:true,gammaIv:.109
});
function run(amplitude,phase,epsilon=0,steps=36){
  const state=makeHexGrowthState({N:97,dx,seedRadius:1.5e-6,
    Tinf:258.15,sigmaInf:.10,kinetics:kin(amplitude,phase,epsilon)});
  for(let step=0;step<steps;step++){
    assert.equal(auditHexDomain(state).stop,false,'stop before domain-bias invalidity');
    resolvedHexGrowthStep(state,.5,{vaporIterations:24,heatIterations:30});
  }
  const morph=hexGrowthMorphology(state,720);
  // Distinguish the fixed-grid X direction from the actual orientation of
  // the *grown* sixfold morphology. Otherwise a 30° rotated intact arm
  // would be falsely scored as a split arm on the original axis.
  function spineScore(angularOrigin,outerOnly=false){
    let center=0,flanks=0,count=0;
    const limit=Math.min(morph.equivalentBasalRadius*2.2,(state.R-12)*dx);
    for(let arm=0;arm<6;arm++){
      const th=angularOrigin+arm*Math.PI/3;
      const armBin=((Math.round(((th%(2*Math.PI)+2*Math.PI)%(2*Math.PI))/
        (2*Math.PI)*morph.radius.length))%morph.radius.length);
      const armExtent=morph.radius[armBin];
      const end=outerOnly?Math.min(limit,armExtent*.95):limit;
      const start=outerOnly?Math.max(2.5e-6,armExtent*.55):2.5e-6;
      for(let radius=start;radius<=end;radius+=dx){
        const at=d=>sampleHexLinear(state,state.phi,radius*Math.cos(th+d),
          radius*Math.sin(th+d));
        const centerVal=at(0),left=at(-Math.PI/18),right=at(Math.PI/18);
        if([centerVal,left,right].every(Number.isFinite)){
          center+=centerVal;flanks+=(left+right)/2;count++;
        }
      }
    }
    return count?(center-flanks)/count:NaN;
  }
  const fixedSpineScore=spineScore(0);
  const alignedSpineScore=spineScore(morph.m6Phase);
  const outerAlignedSpineScore=spineScore(morph.m6Phase,true);
  assert.ok([fixedSpineScore,alignedSpineScore,outerAlignedSpineScore]
    .every(Number.isFinite),'spine/flank measurements require valid radial support');
  return {state,morph,diagnostics:{
    phase,anisotropy:amplitude,epsilon,time_s:state.t,
    mass_pg:state.totalDepositedMass*1e12,
    m6:morph.modes[6].relative,m4:morph.modes[4].relative,
    phase_m6_deg:morph.m6Phase*180/Math.PI,
    fixedAxisCenterMinusFlanks:fixedSpineScore,
    alignedCenterMinusFlanks:alignedSpineScore,
    outerAlignedCenterMinusFlanks:outerAlignedSpineScore,
    halo_cells:auditHexDomain(state).haloCells
  }};
}
const neutral0=run(0,0,0,16),neutral30=run(0,30,0,16);
assert.deepEqual([...neutral0.state.phi],[...neutral30.state.phi],
  'without angular physics, arbitrary phase selection must not change the ice');
assert.equal(neutral0.state.totalDepositedMass,neutral30.state.totalDepositedMass);

const p0=run(.50,0),p15=run(.50,15),p30=run(.50,30);
assert.ok(p0.state.totalDepositedMass>0);
assert.ok(p30.state.totalDepositedMass>0);
assert.ok(p0.state.lastBudget.interfaceCells>0);
assert.ok(p0.state.phi.some((v,k)=>v!==p30.state.phi[k]),
  'changing kinetic phase must actually change phase-field growth');
assert.ok(Math.abs(p0.diagnostics.fixedAxisCenterMinusFlanks-
  p30.diagnostics.fixedAxisCenterMinusFlanks)>1e-7,
  'fixed reference axes should reveal changed orientation');
assert.ok(Math.abs(p30.diagnostics.phase_m6_deg)>20,
  '30° kinetic phase should rotate the measured m6 direction in this test');
const surf0=run(.50,0,0,20),surf1=run(.50,0,.02,20);
assert.ok(surf0.state.phi.some((v,k)=>v!==surf1.state.phi[k]),
  'positive surface-energy anisotropy must modify the interface, not only UI');
console.log('anisotropy-phase.test.mjs PASS',JSON.stringify({
  study:[p0.diagnostics,p15.diagnostics,p30.diagnostics],
  neutral_identical:true,
  capillarityChangesIce:true,
  real_hex_seed_cells:nIce(hex),
  circle_seed_cells:nIce(circle)
},null,2));
