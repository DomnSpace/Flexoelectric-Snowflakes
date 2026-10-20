import assert from 'node:assert/strict';
import {makeGrowth25DState,resolvedGrowthStep} from '../solver/growth25d.js';
import {initializeElectrostatics,solveElectrostatics,dielectricEquilibriumLogShift,VACUUM_PERMITTIVITY} from '../solver/electric25d.js';
const approx=(a,b,tol,msg)=>assert.ok(Math.abs(a-b)<=tol,msg+': '+a+' vs '+b);

const baseline=makeGrowth25DState({N:31,dx:.75e-6,sigmaInf:.08,seedRadius:1.5e-6});
const zero=makeGrowth25DState({N:31,dx:.75e-6,sigmaInf:.08,seedRadius:1.5e-6,
  electric:{enabled:true,fieldVM:0,prescribedGradient:0,coupleGrowth:true}});
assert.ok(zero.electric?.params.enabled,'electric state instantiated');
for(let k=0;k<zero.electric.potential.length;k++){
  assert.ok(zero.electric.potential[k]===0,'no imposed potential');
  assert.ok(zero.electric.magnitude[k]===0,'no field without source');
}
for(let t=0;t<4;t++){
  resolvedGrowthStep(baseline,.5,{vaporIterations:15,heatIterations:15});
  resolvedGrowthStep(zero,.5,{vaporIterations:15,heatIterations:15});
}
approx(zero.totalDepositedMass,baseline.totalDepositedMass,1e-29,'zero-electric baseline mass must match');
approx(zero.basalHalfThickness,baseline.basalHalfThickness,1e-23,'zero-electric basal growth must match');
assert.equal(zero.lastBudget.meanElectricLogShift,0,'no fake forcing at zero E');

// In a uniform gas with a fixed voltage gradient, the discrete ∇·ε∇V
// equation admits the exact linear potential and a constant field.
const uniform=makeGrowth25DState({N:31,dx:.75e-6,seedRadius:0});
uniform.ice.fill(0);uniform.phi.fill(0);
initializeElectrostatics(uniform,{enabled:true,fieldVM:2e6,angleDeg:30,
  epsilonIce:3.2,prescribedGradient:0});
solveElectrostatics(uniform,200);
const c=(uniform.N**2-1)/2;
approx(uniform.electric.Ex[c],2e6*Math.cos(Math.PI/6),1e-5,'constant Ex');
approx(uniform.electric.Ey[c],2e6*Math.sin(Math.PI/6),1e-5,'constant Ey');
assert.ok(uniform.electric.diagnostics.residualRmsV<1e-11,'uniform field Poisson residual');

const e=makeGrowth25DState({N:33,seedRadius:2e-6,
  electric:{enabled:true,fieldVM:0,angleDeg:0,
    prescribedGradient:1500,gradientAngleDeg:0,
    flexoCoefficient:1.65e-9,coupleGrowth:false}});
assert.ok(e.electric.diagnostics.maxPolarization>0,'prescribed flexo source enabled');
assert.ok(e.electric.diagnostics.maxFieldVM>0,'bound polarization charge creates E');
assert.ok(e.electric.potential.every(Number.isFinite),'finite potential');
assert.ok(e.electric.magnitude.every(Number.isFinite),'finite field');
const y=(e.N-1)/2,left=e.electric.potential[y*e.N+3],right=e.electric.potential[y*e.N+e.N-4];
assert.ok(left*right<0,'flexoelectric bound charge creates opposing potential lobes');
assert.equal(dielectricEquilibriumLogShift(e,y*e.N+e.N/2|0,917),0,'diagnostic-only mode cannot alter growth');

// Controlled 10 MV/m laboratory test: enabling the optional fixed-field
// dielectric free-energy closure must generate a measurable signed response,
// not just light up a slider in the viewer.
function sample(coupleGrowth){
  return makeGrowth25DState({N:41,seedRadius:1.5e-6,sigmaInf:.08,
    electric:{enabled:true,fieldVM:1e7,angleDeg:20,epsilonIce:3.2,
      prescribedGradient:800,gradientAngleDeg:75,flexoCoefficient:1.65e-9,
      coupleGrowth}});
}
const passive=sample(false),active=sample(true);
for(let step=0;step<9;step++){
  resolvedGrowthStep(passive,.5,{vaporIterations:20,heatIterations:20});
  resolvedGrowthStep(active,.5,{vaporIterations:20,heatIterations:20});
}
assert.equal(passive.lastBudget.meanElectricLogShift,0);
assert.ok(active.lastBudget.meanElectricLogShift<0,'electric phase stabilization is negative');
assert.ok(active.lastBudget.maxElectricLogShift>0,'nonzero electric correction was evaluated');
assert.ok(active.totalDepositedMass>passive.totalDepositedMass,
  'coupled field must change physical mass deposition in the expected direction');
assert.ok(active.electric.diagnostics.maxFieldVM>0);
assert.ok(Number.isFinite(active.electric.diagnostics.residualRmsV));
assert.ok(active.electric.diagnostics.maxPolarization>0);
assert.ok(VACUUM_PERMITTIVITY>8e-12);
console.log('electric25d.test.mjs PASS', {
  neutral_mass_pg:baseline.totalDepositedMass*1e12,
  passive_mass_pg:passive.totalDepositedMass*1e12,
  active_mass_pg:active.totalDepositedMass*1e12,
  electric_exponent:active.lastBudget.meanElectricLogShift,
  Emax_MVm:active.electric.diagnostics.maxFieldVM/1e6,
  flexo_P_uCm2:active.electric.diagnostics.maxPolarization*1e6
});
