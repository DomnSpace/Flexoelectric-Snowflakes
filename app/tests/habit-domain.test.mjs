import assert from 'node:assert/strict';
import {GROWTH_HABITS,GRID_CHOICES,modelDomainUm,auditGrowthDomain} from '../solver/habits25d.js';
import {makeGrowth25DState,resolvedGrowthStep,growthMorphology} from '../solver/growth25d.js';
import {saturationVaporPressureIceMK,saturationVaporPressureWaterMK} from '../solver/thermo.js';

assert.equal(GROWTH_HABITS.needle.supported,false,
  'needle/column cannot be generated from this single-height basal-plane model');
for(const key of ['plate','sector','dendrite']){
  const h=GROWTH_HABITS[key];
  assert.equal(h.supported,true);
  assert.equal(h.category,'plate');
  assert.ok(GRID_CHOICES.includes(h.gridN),'select supported domain');
  assert.ok(h.basalA>0&&h.basalA<1,'basal accommodation in range');
  assert.ok(h.sigmaPercent>0&&h.temperatureC<0);
  const T=h.temperatureC+273.15;
  const sigmaW=100*(saturationVaporPressureWaterMK(T)/saturationVaporPressureIceMK(T)-1);
  assert.ok(h.sigmaPercent<sigmaW,
    key+' near-water-saturation cloud preset should not require above-water supersaturation');
}
assert.ok(GROWTH_HABITS.dendrite.sigmaPercent>GROWTH_HABITS.sector.sigmaPercent);
assert.ok(GROWTH_HABITS.sector.sigmaPercent>GROWTH_HABITS.plate.sigmaPercent);
assert.equal(modelDomainUm(73,.75e-6),54);
assert.equal(modelDomainUm(129,.75e-6),96);

const seed=makeGrowth25DState({N:73,dx:.75e-6,seedRadius:1.5e-6});
const start=auditGrowthDomain(seed);
assert.ok(start.haloCells>12,'initial seed far from Dirichlet boundary');
assert.equal(start.stop,false,'never stop on the seed');
assert.equal(start.domainUm,54);
assert.ok(start.tipRadiusUm>0);
const thin=auditGrowthDomain({...seed,basalHalfThickness:.2e-6,t:10});
assert.equal(thin.plateStop,false,'thin-plate model remains permitted');

// A sparse protrusion near one edge must trigger even if the equivalent-area
// radius remains tiny; this is the exact failure of old equivalent radius guard.
const N=seed.N;
seed.phi[36*N+(N-1-4)]=.90;
const near=auditGrowthDomain(seed);
assert.equal(near.haloCells,4);
assert.equal(near.stopHalo,true);
assert.equal(near.stop,true);
assert.ok(growthMorphology(seed).equivalentBasalRadius*1e6<8,
  'equivalent radius of a narrow arm must not hide a domain-edge incursion');
seed.phi[36*N+N-1-4]=0;

const bloated=auditGrowthDomain({...seed,basalHalfThickness:8e-6,t:25});
assert.equal(bloated.plateStop,true,'tall prism is not valid as a plate reconstruction');
assert.match(bloated.reason,/SHARED-THICKNESS/);

// The experimentally adjustable basal prefactor must ACTUALLY change the
// out-of-plane growth term, not simply change the HTML label.
function controlled(A){
  return makeGrowth25DState({N:41,Tinf:258.15,sigmaInf:.10,
    kinetics:{prism:{mode:'nucleation',A:1,sigma0:.003},
      basal:{mode:'nucleation',A,sigma0:.01},
      anisotropy6:.2,sdak:false,capillarity:true,gammaIv:.109}});
}
const low=controlled(.005),high=controlled(1);
for(let k=0;k<10;k++){
  resolvedGrowthStep(low,.5,{vaporIterations:24,heatIterations:25});
  resolvedGrowthStep(high,.5,{vaporIterations:24,heatIterations:25});
}
assert.ok(low.lastBudget.vBasal<high.lastBudget.vBasal/50,
  'low attachment prefactor physically suppresses vertical growth rate');
assert.ok(low.basalHalfThickness<high.basalHalfThickness,
  'thickness comes from basal kinetic integration');
for(const x of [low,high]){
  const d=auditGrowthDomain(x);
  assert.ok(Number.isFinite(d.plateRatio));
  assert.ok(d.haloUm>=0&&d.haloUm<=d.domainUm);
}
console.log('habit-domain.test.mjs PASS',{
  cloudCap_pct:100*(saturationVaporPressureWaterMK(258.15)/saturationVaporPressureIceMK(258.15)-1),
  seedHalo:start.haloCells,
  lowThicknessUm:low.basalHalfThickness*2e6,highThicknessUm:high.basalHalfThickness*2e6,
  nearEdgeHalo:near.haloCells
});
