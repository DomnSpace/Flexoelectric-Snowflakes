import assert from 'node:assert/strict';
import {makeHexGrowthState,resolvedHexGrowthStep} from '../solver/hex25d.js';
import {iceIhProperties} from '../solver/ice-ih-eos.js';
import {versionLocalIceSnapshot} from '../solver/ice-snapshot-v2.js';
import {computeActualHeightGeometry,computeActualHeightSlices,
  HEIGHT_SLICE_LEVELS,selectHeightMeshStride,PHASE_ISO} from '../render/actual-height-mesh.js';

const rho=iceIhProperties(258.15,101325).rho;
function frameOf(state){
  const N=state.N;
  return versionLocalIceSnapshot({
    n:N,dx:state.dx,t:state.t,topology:'triangular',
    mask:new Uint8Array(state.mask),ice:new Uint8Array(state.ice),
    phi:new Float32Array(state.phi),u:new Float32Array(state.sigma),
    dT:new Float32Array(state.dT),morph:{
      radius:Array.from({length:180},()=>3e-6),
      halfThickness:state.basalHalfThickness
    }
  },state,rho);
}
function physicalSignature(state){
  return [
    state.t,state.totalDepositedMass,state.totalLatentHeat,
    ...state.phi,...state.topHeight,...state.bottomHeight,
    ...state.sigma,...state.dT
  ];
}
function areaXY(positions){
  let sum=0;
  for(let i=0;i+8<positions.length;i+=9){
    const x0=positions[i],y0=positions[i+1],
      x1=positions[i+3],y1=positions[i+4],
      x2=positions[i+6],y2=positions[i+7];
    sum+=Math.abs((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0))/2;
  }
  return sum;
}
function isFiniteData(a){
  assert.equal(a.length%9,0,'a complete set of triangles');
  assert.ok(a.every(Number.isFinite),'finite triangular positions');
}
const state=makeHexGrowthState({N:41,dx:.75e-6,seedRadius:0,
  localHeight:true,Tinf:258.15,sigmaInf:.08});
const N=state.N,R=state.R,idx=(i,j)=>j*N+i;
for(let j=0;j<N;j++)for(let i=0;i<N;i++){
  const k=idx(i,j),q=i-R,r=j-R;
  if(!state.mask[k])continue;
  // Construct a known annular crystal with an ACTUAL central hole, and
  // higher central cross section through its thicker arm region.
  const radius=Math.hypot(q+r/2,r*Math.sqrt(3)/2);
  if(radius<3||radius>12)continue;
  state.phi[k]=1;state.ice[k]=1;
  const height=(radius>=6&&radius<9?4:2)*1e-6;
  state.topHeight[k]=height;state.bottomHeight[k]=-height;
}
const input=frameOf(state);
const before=physicalSignature(state);
assert.ok(input.iceStructure.maxLocalThicknessM>7e-6);
const full=computeActualHeightGeometry(input);
assert.equal(full.stride,1);
assert.ok(full.occupiedTriangles>50);
assert.ok(full.wallSegments>10,'inside hole and outside rim both have walls');
for(const arr of [full.top,full.bottom,full.walls])isFiniteData(arr);
assert.equal(full.volumeExactlyConserved,false);
assert.equal(full.source,'recorded snapshot-v2 top/bottom/phase fields');
assert.ok(full.top.every((z,i)=>i%3!==2||z>=0));
assert.ok(full.bottom.every((z,i)=>i%3!==2||z<=0));

const five=computeActualHeightSlices(input);
assert.equal(five.layers.length,5);
assert.deepEqual(HEIGHT_SLICE_LEVELS,[-.8,-.4,0,.4,.8]);
assert.ok(five.contourDerivedFromLocalHeight);
assert.equal(five.includesDisconnectedRegions,true);
assert.equal(five.preservesHolesWhereResolved,true);
for(const layer of five.layers){
  isFiniteData(layer.positions);
  assert.equal(layer.positions.length,layer.triangleCount*9);
  assert.ok(layer.positions.every((v,k)=>k%3!==2||
    Math.abs(v-layer.zM*1e6)<1e-4));
}
const areas=five.layers.map(l=>areaXY(l.positions));
assert.ok(areas[0]>0&&areas[2]>areas[0]*1.4,
  'middle plane includes both thin and thick annular material');
assert.ok(areas[4]>0&&areas[2]>areas[4]*1.4,
  'top plane includes only genuinely tall crystal regions');
assert.ok(five.layers.every(l=>l.triangleCount>0));
const areaTop=areaXY(five.layers[4].positions);
assert.ok(areaTop<areaXY(five.layers[2].positions));
assert.deepEqual(physicalSignature(state),before,
  'meshing must neither create ice mass nor mutate height or gas arrays');

// A thin uniform slab produces five identical XY cross-sections: no ridges
// are invented merely by using the more detailed mesh.
for(let k=0;k<state.phi.length;k++){
  if(state.phi[k]===0)continue;
  state.topHeight[k]=4e-6;state.bottomHeight[k]=-4e-6;
}
const uniform=computeActualHeightSlices(frameOf(state));
const uniformAreas=uniform.layers.map(l=>areaXY(l.positions));
assert.ok(uniformAreas.every(a=>Math.abs(a-uniformAreas[0])<1e-4),
  'equal heights cannot produce artificial ridges');
assert.equal(selectHeightMeshStride(frameOf(state)),1);

// Real running pass-8 state: physics is unchanged by meshing or changing view.
const grown=makeHexGrowthState({N:73,dx:.75e-6,Tinf:258.15,
  sigmaInf:.09,seedRadius:1.5e-6,localHeight:true});
for(let i=0;i<3;i++)resolvedHexGrowthStep(grown,.5,{
  vaporIterations:18,heatIterations:18
});
const gBefore=physicalSignature(grown),frame=frameOf(grown);
const geom=computeActualHeightGeometry(frame),layers=computeActualHeightSlices(frame);
assert.ok(geom.top.length>0&&geom.bottom.length>0);
assert.ok(layers.layers[2].triangleCount>0);
assert.deepEqual(physicalSignature(grown),gBefore);

// Large domain: graphics are generated at documented LOD only; mass and
// grid physics stay at full 0.75µm solver resolution.
const big=makeHexGrowthState({N:401,dx:.75e-6,localHeight:true});
assert.equal(selectHeightMeshStride(frameOf(big)),2);
const beforeLarge=physicalSignature(big);
const bigMesh=computeActualHeightGeometry(frameOf(big));
assert.equal(bigMesh.stride,2);
assert.ok(bigMesh.occupiedTriangles>0);
assert.deepEqual(physicalSignature(big),beforeLarge);
console.log('pass9-actual-height.test.mjs PASS',{
  fullTriangles:full.occupiedTriangles,wallSegments:full.wallSegments,
  fiveSliceAreasUm2:areas.map(a=>Number(a.toFixed(3))),
  equalHeightLayerArea:uniformAreas[0],bigRenderLOD:bigMesh.stride,
  seededLedgerPhysicalMassKg:input.iceStructure.iceMassKg
});