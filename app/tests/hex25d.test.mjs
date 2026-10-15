import assert from 'node:assert/strict';
import {makeHexGrowthState,HEX_NEIGHBORS,HEX_CELL_AREA_MULTIPLIER,
  axialToXY,axialRadius,axialRoundFromXY,
  hexCellAtXY,hexGrowthMorphology,relaxHexScalar,resolvedHexGrowthStep,
  auditHexDomain,enlargeHexDomain} from '../solver/hex25d.js';

const a=.75e-6,N=73,R=(N-1)/2;
const s=makeHexGrowthState({N,dx:a,Tinf:258.15,sigmaInf:.10});
assert.equal(s.topology,'triangular');
assert.equal(s.activeSites,1+3*R*(R+1));
assert.equal(HEX_NEIGHBORS.length,6);
assert.equal(new Set(HEX_NEIGHBORS.map(d=>d.join(':'))).size,6);
assert.ok(Math.abs(HEX_CELL_AREA_MULTIPLIER-Math.sqrt(3)/2)<1e-12);

const seen=[];
for(const [q,r] of [[0,0],[1,0],[0,1],[1,-1],[-3,8],[R,0],[-R,R]]){
  const pos=axialToXY(q,r,a),rounded=axialRoundFromXY(pos.x,pos.y,a);
  assert.deepEqual(rounded,{q,r});
  if(axialRadius(q,r)<=R)assert.ok(hexCellAtXY(s,pos.x,pos.y)>=0);
  seen.push(rounded);
}
assert.equal(hexCellAtXY(s,(R+3)*a,0),-1);

// The 6-neighbor triangular transport stencil is rotationally isotropic
// through second order. It must give exactly Δf=0 for linear f and
// Δ(x²+y²)=4 for quadratic f regardless of direction.
function laplacian(q,r,field){
  const p=axialToXY(q,r,a),sum=HEX_NEIGHBORS.reduce((tot,[dq,dr])=>{
    const k=axialToXY(q+dq,r+dr,a);return tot+field(k.x,k.y);
  },0);
  return 2/(3*a*a)*(sum-6*field(p.x,p.y));
}
for(const [q,r] of [[0,0],[2,-3],[-8,5]]){
  assert.ok(Math.abs(laplacian(q,r,(x,y)=>x+y*3))<1e-7);
  assert.ok(Math.abs(laplacian(q,r,(x,y)=>x*x+y*y)-4)<1e-10);
}
const start=auditHexDomain(s);
assert.equal(start.stop,false);
assert.equal(start.domainUm,54);
assert.ok(start.haloCells>12);
assert.ok(start.thicknessUm>0);
assert.ok(s.mask.reduce((a,b)=>a+b,0)<N*N,'domain corners are excluded');

// Hex crystal grows with actual mass and temperature source.
let mass=0;
for(let step=0;step<12;step++){
  resolvedHexGrowthStep(s,.5,{vaporIterations:32,heatIterations:38});
  assert.ok(Number.isFinite(s.totalDepositedMass));
  assert.ok(s.totalDepositedMass>=mass);
  mass=s.totalDepositedMass;
}
assert.ok(mass>0);
assert.ok(s.lastBudget.interfaceCells>0);
assert.ok(s.lastBudget.maxDeltaT>0);
assert.equal(s.lastBudget.stencil,'six-equidistant-triangular');
assert.equal(s.lastBudget.meanElectricLogShift,0);
const m=hexGrowthMorphology(s,180);
assert.ok(m.equivalentBasalRadius>0&&Number.isFinite(m.modes[6].relative));
assert.ok(Number.isFinite(m.modes[4].relative));
assert.ok(Number.isFinite(m.aspectRatio));
assert.ok(s.sigma.every(Number.isFinite));
assert.ok(s.dT.every(Number.isFinite));

// Boundary geometry is hexagonal, not the bounding square or equivalent area.
const c=R,idx=(i,j)=>j*N+i;
s.ice[idx(R+R-4,R)]=1;
const near=auditHexDomain(s);
assert.equal(near.haloCells,4);
assert.equal(near.stop,true);

// Conservatively preserve occupied physical cells, morphology mass, time,
// energy and phase through a domain increase, without restarting the flake.
const beforeArea=s.phi.reduce((t,x,i)=>t+(s.mask[i]?x:0),0);
const expanded=enlargeHexDomain(s,97);
assert.equal(expanded.N,97);
assert.equal(expanded.activeSites,1+3*48*49);
assert.equal(expanded.totalDepositedMass,s.totalDepositedMass);
assert.equal(expanded.totalLatentHeat,s.totalLatentHeat);
assert.equal(expanded.t,s.t);
assert.equal(expanded.domainExpansions,1);
const afterArea=expanded.phi.reduce((t,x,i)=>t+(expanded.mask[i]?x:0),0);
assert.ok(Math.abs(beforeArea-afterArea)<1e-12);
assert.ok(auditHexDomain(expanded).haloCells>near.haloCells);
assert.ok(expanded.sigma.every(Number.isFinite));
resolvedHexGrowthStep(expanded,.5,{vaporIterations:12,heatIterations:14});
assert.ok(expanded.totalDepositedMass>=mass);

// Square-free mesh still has SIXFOLD GRID SYMMETRY. This is a numerical
// property, not evidence that crystal sixfold morphology is physical.
console.log('hex25d.test.mjs PASS',{
  sites:s.activeSites,sourceMassPg:mass*1e12,
  m4:m.modes[4].relative,m6:m.modes[6].relative,
  haloBefore:near.haloCells,haloExpanded:auditHexDomain(expanded).haloCells
});
