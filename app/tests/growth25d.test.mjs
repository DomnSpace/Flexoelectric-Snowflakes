import assert from 'node:assert/strict';
import { makeGrowth25DState, resolvedGrowthStep, growthMorphology, fieldSample } from '../solver/growth25d.js';

const s=makeGrowth25DState({N:41,dx:0.5e-6,seedRadius:1.5e-6,Tinf:258.15,sigmaInf:0.08,basalHalfThickness:1.0e-6});
const m0=growthMorphology(s);
let lastMass=0;
for(let n=0;n<20;n++){
  resolvedGrowthStep(s,0.5,{vaporIterations:25,heatIterations:30});
  assert.ok(s.totalDepositedMass>=lastMass,'deposited mass monotone');
  lastMass=s.totalDepositedMass;
}
const m1=growthMorphology(s);
assert.ok(s.totalDepositedMass>0,'positive vapor deposition');
assert.ok(s.totalLatentHeat>0,'latent heat accumulated');
assert.ok(s.lastBudget.interfaceCells>0,'resolved interface exists');
assert.ok(s.lastBudget.meanPrismAlpha>=0&&s.lastBudget.meanPrismAlpha<=1);
assert.ok(s.lastBudget.alphaBasal>=0&&s.lastBudget.alphaBasal<=1);
assert.ok(m1.meanRadius>=m0.meanRadius,'basal-plane radius cannot shrink');
assert.ok(m1.halfThickness>m0.halfThickness,'basal thickness should grow');
assert.ok(Number.isFinite(m1.m6Amplitude));
assert.ok(Number.isFinite(m1.aspectRatio)&&m1.aspectRatio>0);

const c=(s.N-1)/2;
const center=fieldSample(s,c,c);
assert.equal(center.ice,true);
let hot=false,depleted=false;
for(let j=1;j<s.N-1;j++)for(let i=1;i<s.N-1;i++){
  const q=fieldSample(s,i,j);
  if(q.deltaT>0)hot=true;
  if(!q.ice&&q.sigma<s.sigmaInf)depleted=true;
}
assert.ok(hot,'latent heat should generate a positive temperature field');
assert.ok(depleted,'ice should deplete local supersaturation below far-field value');

console.log('growth25d.test.mjs PASS',{deposited_pg:s.totalDepositedMass*1e12,aspect:m1.aspectRatio,m6_um:m1.m6Amplitude*1e6});
