/**
 * Milestone 1/15: deterministic, READ-ONLY legacy solver fingerprint.
 *
 * Run: node app/tests/legacy-baseline-capture.mjs
 *
 * Generates reference values for the 2026-10-15 pre-coupling baseline.
 * This module DOES NOT alter the imported solvers or render five fake fields.
 * CI node/runtime and frozen commit must be documented alongside digests.
 */
import {createHash} from 'node:crypto';
import {makeHexGrowthState,resolvedHexGrowthStep,hexGrowthMorphology}
  from '../solver/hex25d.js';
import {makeGrowth25DState,resolvedGrowthStep,growthMorphology}
  from '../solver/growth25d.js';
import {DEFAULT_ICE_KINETICS} from '../solver/growth25d.js';

const BASELINE_COMMIT='77f118c091a2fd1ae278c899e7d5f6a1715ebc54';
const options=Object.freeze({
  N:41,dx:0.75e-6,Tinf:258.15,sigmaInf:0.075,
  pressurePa:101000,seedRadius:1.5e-6,basalHalfThickness:1.5e-6
});
const fieldNames=['ice','phi','sigma','dT','heatSource'];

function fieldSHA256(state) {
  const hash=createHash('sha256');
  for(const name of fieldNames){
    const field=state[name];
    if(!ArrayBuffer.isView(field))throw new Error('Missing legacy typed-array field '+name);
    hash.update(name+'\0'+field.constructor.name+'\0'+field.length+'\0');
    hash.update(Buffer.from(field.buffer,field.byteOffset,field.byteLength));
  }
  return hash.digest('hex');
}

function runCase(mode){
  const hex=mode!=='cartesian';
  const kinetics=JSON.parse(JSON.stringify(DEFAULT_ICE_KINETICS));
  if(mode==='hex30'){
    kinetics.anisotropy6=0.5;
    kinetics.kineticPhaseDeg=30;
    kinetics.surfaceEnergyEpsilon=0;
    kinetics.surfacePhaseDeg=30;
  }
  const state=(hex?makeHexGrowthState:makeGrowth25DState)({
    ...options,kinetics, ...(hex?{seedShape:'circle',seedRotationDeg:0}:{})
  });
  const step=hex?resolvedHexGrowthStep:resolvedGrowthStep;
  const run=[];const snapshots=[];
  for(let iter=1;iter<=12;iter++){
    step(state,0.5,{vaporIterations:14,heatIterations:18});
    if([1,6,12].includes(iter)){
      const m=hex?hexGrowthMorphology(state,180):growthMorphology(state,180);
      run.push({
        step:iter,t:state.t,
        depositedMassKg:state.totalDepositedMass,
        totalLatentHeatJ:state.totalLatentHeat,
        halfThicknessM:state.basalHalfThickness,
        m6:m.modes?.[6]?.relative??null,
        interfaceCells:state.lastBudget?.interfaceCells??null
      });
      snapshots.push({step:iter,fieldsSHA256:fieldSHA256(state)});
    }
  }
  const finite=fieldNames.every(name=>Array.from(state[name]).every(Number.isFinite));
  if(!finite)throw new Error('Non-finite legacy reference field in '+mode);
  return {case:mode,topology:hex?'triangular':'square',
    gridN:state.N,gridSpacingM:state.dx,steps:12,dtSeconds:0.5,
    vaporRelaxationsPerStep:14,heatRelaxationsPerStep:18,
    snapshots,measurements:run};
}
export function captureLegacyBaseline(){
  return {
    schema:'snowflake/legacy-baseline-v1',
    legacyCommit:BASELINE_COMMIT,
    traceUnits:'SI',
    platformNote:'Float64 little-endian SHA256 fingerprints generated in Node 22 CI; compare within same runtime architecture',
    cases:['hex0','hex30','cartesian'].map(runCase)
  };
}
if(import.meta.url===new URL(process.argv[1]||'', 'file://').href){
  console.log('FIVE_LAYER_LEGACY_BASELINE_JSON_START');
  console.log(JSON.stringify(captureLegacyBaseline(),null,2));
  console.log('FIVE_LAYER_LEGACY_BASELINE_JSON_END');
}
