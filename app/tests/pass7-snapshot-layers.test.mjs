import assert from 'node:assert/strict';
import {makeHexGrowthState,resolvedHexGrowthStep,hexGrowthMorphology} from '../solver/hex25d.js';
import {makeGrowth25DState,resolvedGrowthStep,growthMorphology} from '../solver/growth25d.js';
import {versionIceSnapshot,verifyIceSnapshot,upgradeLegacyIceSnapshot,
  legacyHeightAt,ICE_SNAPSHOT_SCHEMA,ICE_SNAPSHOT_VERSION} from '../solver/ice-snapshot-v1.js';
import {makeUniformLayerPlan,FIVE_LAYER_COUNT,FIVE_LAYER_OPACITY} from '../render/uniform-five-layers.js';

/** Mirror the Observatory's read-only snapshot capture, not its rendering. */
function snapshot(state,hex=false){
  const morph=hex?hexGrowthMorphology(state,180):growthMorphology(state,180);
  return versionIceSnapshot({
    t:state.t,n:state.N,dx:state.dx,Tinf:state.Tinf,
    topology:hex?'triangular':'square',
    mask:hex?new Uint8Array(state.mask):null,
    u:new Float32Array(state.sigma),dT:new Float32Array(state.dT),
    phi:new Float32Array(state.phi),ice:new Uint8Array(state.ice),
    morph:{radius:morph.radius,halfThickness:morph.halfThickness},
    budget:state.lastBudget?{...state.lastBudget}:null,
    deposited:state.totalDepositedMass,energy:state.totalLatentHeat
  });
}
function comparePhysics(a,b){
  assert.equal(a.t,b.t);
  assert.equal(a.totalDepositedMass,b.totalDepositedMass);
  assert.equal(a.totalLatentHeat,b.totalLatentHeat);
  assert.equal(a.basalHalfThickness,b.basalHalfThickness);
  for(const field of ['phi','ice','sigma','dT','heatSource']){
    assert.deepEqual(a[field],b[field],'Preview must not alter '+field);
  }
}
for(const hex of [true,false]){
  const opts={N:41,dx:.75e-6,Tinf:258.15,sigmaInf:.075,seedRadius:1.5e-6};
  const create=hex?makeHexGrowthState:makeGrowth25DState;
  const step=hex?resolvedHexGrowthStep:resolvedGrowthStep;
  const model=create(opts),control=create(opts);
  for(let k=0;k<12;k++){
    step(model,.5,{vaporIterations:14,heatIterations:18});
    const frame=snapshot(model,hex);
    assert.equal(frame.snapshotSchema.name,ICE_SNAPSHOT_SCHEMA);
    assert.equal(frame.snapshotSchema.version,ICE_SNAPSHOT_VERSION);
    assert.equal(frame.iceStructure.representation,'legacy-uniform-slab');
    assert.equal(frame.iceStructure.localHeightSolved,false);
    assert.equal(frame.iceStructure.formationTimeSolved,false);
    assert.equal(frame.iceStructure.basalMassCoupled,false);
    assert.equal(frame.iceStructure.physicalTotalThicknessM,2*model.basalHalfThickness);
    assert.ok(verifyIceSnapshot(frame));
    const plan=makeUniformLayerPlan(frame);
    assert.equal(plan.sliceCount,FIVE_LAYER_COUNT);
    assert.equal(plan.layers.length,5);
    assert.equal(plan.mode,'uniform-five-layer-preview');
    assert.equal(plan.physicalLocalHeightSolved,false);
    assert.equal(plan.eachSliceOpacity,FIVE_LAYER_OPACITY);
    assert.ok(plan.layers.every((layer,i)=>
      layer.contourRadii===frame.morph.radius&&
      layer.zM===frame.iceStructure.physicalHalfThicknessM*([-1,-.5,0,.5,1][i])));
    assert.ok(plan.layers.every(layer=>layer.contourRadii.length===frame.morph.radius.length));
    assert.equal(frame.iceStructure.physicalTotalThicknessM,
      plan.layers[4].zM-plan.layers[0].zM);
    assert.equal(frame.iceStructure.localHeightSolved,false);
    const frozenIndex=frame.ice.findIndex(x=>x===1);
    assert.ok(frozenIndex>=0);
    const h=legacyHeightAt(frame,frozenIndex);
    assert.equal(h.topM,model.basalHalfThickness);
    assert.equal(h.bottomM,-model.basalHalfThickness);
    assert.equal(h.iceFraction,1);
    assert.throws(()=>versionIceSnapshot(frame),/already versioned/);
    assert.strictEqual(upgradeLegacyIceSnapshot(frame),frame);
    const legacy={...frame};delete legacy.snapshotSchema;delete legacy.iceStructure;
    assert.equal(upgradeLegacyIceSnapshot(legacy).snapshotSchema.version,1);
    const before=snapshot(model,hex);
    assert.deepEqual(frame.phi,before.phi);
    // Snapshot data must be copied, not aliased into the solver.
    frame.phi[frozenIndex]=.3;
    assert.equal(model.phi[frozenIndex],1);
    step(control,.5,{vaporIterations:14,heatIterations:18});
    comparePhysics(model,control);
  }
}
const invalid=makeHexGrowthState({N:41});
const raw=snapshot(invalid,true);
assert.throws(()=>verifyIceSnapshot({...raw,snapshotSchema:{...raw.snapshotSchema,version:99}}),/Unknown or unversioned/);
assert.throws(()=>verifyIceSnapshot({...raw,iceStructure:{...raw.iceStructure,localHeightSolved:true}}),/cannot claim/);
assert.throws(()=>legacyHeightAt(raw,-1),RangeError);
console.log('pass7-snapshot-layers.test.mjs PASS · 12 parity steps × hex and Cartesian, 5 slices, copied versioned snapshots');
