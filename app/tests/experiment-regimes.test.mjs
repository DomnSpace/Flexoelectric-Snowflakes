import assert from 'node:assert/strict';
import {EXPERIMENT_PRESETS,EXPERIMENT_VERSION,waterSaturationIcePercent,
  cloudCondition,validatePreset,sampleGrowthHistory,makeExperimentRun,
  makeCustomExperimentRun,applyGrowthEnvironment,cloneExperimentState,
  observedHabitTendency} from '../solver/experiment-regimes.js';
import {makeHexGrowthState,resolvedHexGrowthStep} from '../solver/hex25d.js';
import {saturationMassDensityIce} from '../solver/postseed.js';
import {iceIhProperties} from '../solver/ice-ih-eos.js';
import {localIceMass} from '../solver/ice-thickness.js';

const close=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<tol*Math.max(Math.abs(a),Math.abs(b),1e-22),'Expected '+a+' ~= '+b);
assert.equal(EXPERIMENT_VERSION,1);
assert.match(observedHabitTendency(-2),/plate/i);
assert.match(observedHabitTendency(-5),/needle/i);
assert.match(observedHabitTendency(-15),/dendrite/i);
assert.match(observedHabitTendency(-26),/vary/i);
close(waterSaturationIcePercent(-15),15.74,.001);
close(waterSaturationIcePercent(-5),4.98,.003);
assert.ok(cloudCondition(-15,14).relativeHumidityOverWaterPercent<100);
assert.ok(!cloudCondition(-15,14).aboveWaterSaturation);
assert.ok(cloudCondition(-5,18).aboveWaterSaturation);
assert.ok(!cloudCondition(-5,4).aboveWaterSaturation);
assert.equal(cloudCondition(-5,18).classification,'above-water-saturation: laboratory or transient only');
assert.equal(EXPERIMENT_PRESETS.needle.supported,false);
assert.equal(EXPERIMENT_PRESETS.cappedColumn.supported,false);
assert.throws(()=>makeExperimentRun('needle'),/Unsupported/);
for(const preset of Object.values(EXPERIMENT_PRESETS)){
  assert.ok(validatePreset(preset));
  const c=preset.controls;
  assert.ok([0,15,30].includes(c.kineticPhaseDeg));
  assert.ok(c.basalSigma0>0&&c.basalA>0);
  assert.ok(c.seedRadiusUm>0&&c.gridN%2===1);
  assert.ok(!c.electricEnabled&&!c.electricCouple);
  assert.match(preset.calibration,/uncalibrated/);
}
const timeline=makeExperimentRun('history');
assert.equal(timeline.stages.length,4);
assert.equal(timeline.stages[0].sigmaPercent,3);
assert.equal(sampleGrowthHistory(timeline.stages,14.9).sigmaPercent,3);
assert.equal(sampleGrowthHistory(timeline.stages,15).sigmaPercent,3);
close(sampleGrowthHistory(timeline.stages,16.5).sigmaPercent,8);
assert.equal(sampleGrowthHistory(timeline.stages,18).sigmaPercent,13);
assert.equal(sampleGrowthHistory(timeline.stages,45).sigmaPercent,13);
assert.equal(sampleGrowthHistory(timeline.stages,48).sigmaPercent,5);
assert.equal(sampleGrowthHistory(timeline.stages,68).sigmaPercent,14);
assert.equal(sampleGrowthHistory(timeline.stages,69).index,3);
assert.deepEqual(timeline.transitions,[]);
assert.equal(makeCustomExperimentRun(-13,8).id,'custom');

// Verify a live steady→ramp transition changes the real far-field BC and
// preserves the absolute interstitial c_v during the reference-temperature
// rebase. This does not claim conservation at an OPEN Dirichlet shell.
const state=makeHexGrowthState({
  N:41,dx:.75e-6,seedRadius:1.5e-6,localHeight:true,
  Tinf:258.15,sigmaInf:.03,
  pressurePa:101000
});
for(let i=0;i<2;i++)resolvedHexGrowthStep(state,.5,{
  vaporIterations:20,heatIterations:20
});
const originalTime=state.t,originalPhi=new Float64Array(state.phi),
  originalEnergy=state.totalLatentHeat,
  originalT=state.Tinf,oldC=saturationMassDensityIce(originalT),
  oldRho=iceIhProperties(originalT,101000).rho,
  originalMass=localIceMass(state,oldRho);
const center=state.R*state.N+state.R+7;
assert.ok(state.mask[center]&&!state.ice[center]);
const cOld=oldC*(1+state.sigma[center]);
const oldDeltaT=state.dT[center];
const customHistory={...makeCustomExperimentRun(-15,3),id:'explicit-two-stage',
  stages:[{t:0,label:'cold',temperatureC:-15,sigmaPercent:3},
    {t:5,label:'warmer',temperatureC:-12,sigmaPercent:6}]};
const before=applyGrowthEnvironment(state,customHistory,4);
assert.equal(before.changed,false);
const transition=applyGrowthEnvironment(state,customHistory,6);
assert.ok(transition.changed);
assert.ok(transition.rampFraction>0&&transition.rampFraction<1);
assert.equal(transition.index,1);
assert.equal(customHistory.transitions.length,1);
assert.equal(state.t,originalTime);
assert.deepEqual(state.phi,originalPhi);
assert.equal(state.totalLatentHeat,originalEnergy);
assert.equal(state.dT[center],oldDeltaT);
const newC=saturationMassDensityIce(state.Tinf);
close(newC*(1+state.sigma[center]),cOld,1e-12);
const nowRho=iceIhProperties(state.Tinf,101000).rho;
close(localIceMass(state,nowRho),originalMass,5e-12);
for(const k of state.edges)close(state.sigma[k],state.sigmaInf);
assert.equal(state.sigmaInf,transition.sigmaPercent/100);
const replayCopy=cloneExperimentState(customHistory);
replayCopy.stages[1].sigmaPercent=32;
assert.notEqual(customHistory.stages[1].sigmaPercent,32);
replayCopy.transitions[0].time=999;
assert.equal(customHistory.transitions[0].time,6);

// Constant conditions never rebase or restart the gas.
const steady=makeExperimentRun('dendrite');
const noChange=makeHexGrowthState({N:41,Tinf:258.15,sigmaInf:.14});
const initialU=new Float64Array(noChange.sigma);
assert.equal(applyGrowthEnvironment(noChange,steady,85).changed,false);
assert.deepEqual(noChange.sigma,initialU);
console.log('experiment-regimes.test.mjs PASS',{
  tempWaterCapMinus15:waterSaturationIcePercent(-15),
  stageSequence:timeline.stages.map(s=>[s.t,s.sigmaPercent]),
  absoluteVaporRebased:true,heightMassRebased:true,
  uncalibratedFacetCoefficients:true
});