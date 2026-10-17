/**
 * Model B, Step 3 — five independently evolving ice-volume-fraction slabs.
 *
 * Only the five ice-eligible z slabs may solidify. Each slab has its own
 * phase field, stored ice mass and deposition history. Gas fractions shrink
 * as actual water-vapor molecules attach to exposed surfaces; the displaced
 * gas is transferred to adjacent accessible gas cells without mass creation.
 *
 * This uses a coarse voxel / subcell-volume closure, NOT molecular ice and
 * NOT a resolved smooth prism/basal interface. No new preference for the
 * middle layer is imposed. Initial seeds may be symmetric across all 5.
 *
 * The standard Observatory / hex25d legacy solver remains untouched.
 */
import {createLayeredGasFV,layeredGasIndex,axialHexRadius}
  from './layered-gas-fv.js';
import {layeredImplicitGasMassKg,layeredImplicitTrackedWaterKg,
  stepLayeredImplicit} from './layered-vapor-implicit.js';
import {iceIhProperties} from './ice-ih-eos.js';
import {attachLayeredThermal,stepLayeredThermal} from './layered-thermal.js';

const DIR=[[1,0],[-1,0],[0,1],[0,-1],[1,-1],[-1,1]];
const closeTo=(a,b)=>Math.abs(a-b)<=Math.max(1e-26,1e-8*Math.max(a,b));
function neighbors(s,k){
  const z=Math.floor(k/s.nxy),p=k-z*s.nxy;
  const i=p%s.N,j=Math.floor(p/s.N),result=[];
  for(const [dq,dr] of DIR){
    const x=i+dq,y=j+dr;
    if(x>=0&&y>=0&&x<s.N&&y<s.N&&s.hexMask[y*s.N+x])
      result.push(z*s.nxy+y*s.N+x);
  }
  for(const zz of [z-1,z+1])if(zz>=0&&zz<s.nz)
    result.push(zz*s.nxy+p);
  return result;
}
function inEligible(s,k){
  const z=Math.floor(k/s.nxy);
  return z>=s.gasHaloSlabs&&z<s.gasHaloSlabs+s.iceSlabs;
}
export function createFiveLayerCrystal({
  N=21,dx=.75e-6,dz=.75e-6,gasHaloSlabs=2,seedRadiusCells=1,
  seedLayerIndices=[0,1,2,3,4],
  Tinf=258.15,pressurePa=101325,sigmaInf=.14,
  thermalCoupling=false,
  prismAttachmentMPerS=5e-3,basalAttachmentMPerS=5e-3,
  boundaries={lateral:'reservoir',top:'reservoir',bottom:'reservoir'},
  ...rest
}={}){
  if(!Number.isInteger(seedRadiusCells)||seedRadiusCells<0||seedRadiusCells>5)
    throw new RangeError('Seed radius in axial sites must be 0..5');
  if(!Array.isArray(seedLayerIndices)||!seedLayerIndices.length||
    seedLayerIndices.some(l=>!Number.isInteger(l)||l<0||l>=5)||
    new Set(seedLayerIndices).size!==seedLayerIndices.length)
    throw new RangeError('Choose unique seed layers within the five ice slabs');
  const r=(N-1)/2,iceCells=[];
  for(const layer of seedLayerIndices){
    const z=gasHaloSlabs+layer;
    for(let j=0;j<N;j++)for(let i=0;i<N;i++){
      if(axialHexRadius(i-r,j-r)<=seedRadiusCells)
        iceCells.push({q:i-r,r:j-r,z});
    }
  }
  const s=createLayeredGasFV({N,dx,dz,iceSlabs:5,gasHaloSlabs,
    Tinf,pressurePa,sigmaInf,prismAttachmentMPerS,
    basalAttachmentMPerS,boundaries,iceCells,...rest});
  s.seedSolidMask=new Uint8Array(s.solid);
  s.iceFraction=new Float64Array(s.solid.length);
  s.iceMassAtCellKg=new Float64Array(s.solid.length);
  s.rhoIce=iceIhProperties(Tinf,Math.min(pressurePa,210e6)).rho;
  s.cellIceCapacityKg=s.rhoIce*s.volume;
  s.seedIceMassKg=iceCells.length*s.cellIceCapacityKg;
  s.trappedVaporKg=0;
  s.schema='snowflake/vertical-gas-fv-step2'; // retain transport schema for solver
  s.model='snowflake/model-b-five-independent-slabs-v1';
  s.interfaceMotionSolved=true; // voxel occupation, NOT smooth moving facets
  for(let k=0;k<s.solid.length;k++)if(s.solid[k]){
    s.iceFraction[k]=1;
    s.iceMassAtCellKg[k]=s.cellIceCapacityKg;
  }
  if(thermalCoupling)attachLayeredThermal(s);
  return s;
}
export function fiveLayerIceMassKg(s){
  let t=0;for(const v of s.iceMassAtCellKg)t+=v;return t;
}
export function fiveLayerOccupancy(s){
  const out=[];
  for(let l=0;l<s.iceSlabs;l++){
    const first=(s.gasHaloSlabs+l)*s.nxy;
    out.push(new Float64Array(s.iceFraction.subarray(first,first+s.nxy)));
  }
  return out;
}
export function fiveLayerDiagnostics(s){
  const massKg=[],occupiedSites=[],newGrowthSites=[],gasMassKg=[];
  for(let l=0;l<s.iceSlabs;l++){
    let iceMass=0,occupied=0,grown=0,gas=0;
    const z=s.gasHaloSlabs+l;
    for(let p=0;p<s.nxy;p++){
      if(!s.hexMask[p])continue;
      const k=z*s.nxy+p,ph=s.iceFraction[k];
      if(ph>0)occupied++;
      if(ph>0&&!s.seedSolidMask[k])grown++;
      iceMass+=s.iceMassAtCellKg[k];
      gas+=s.c[k]*s.volume*(1-ph);
    }
    massKg.push(iceMass);occupiedSites.push(occupied);
    newGrowthSites.push(grown);gasMassKg.push(gas);
  }
  return {massKg,occupiedSites,newGrowthSites,gasMassKg,
    depositedKg:s.totalDepositedKg,trappedVaporKg:s.trappedVaporKg,
    vaporInventoryKg:layeredImplicitGasMassKg(s),
    waterInventoryKg:layeredImplicitTrackedWaterKg(s)};
}
/** Reallocate gas displaced by growing solid fractions, simultaneously. */
export function settleFiveLayerFront(s,events){
  if(!s.iceFraction||!s.iceMassAtCellKg)
    throw new TypeError('Not a moving five-layer state');
  const prev=s.iceFraction.slice(),dMass=new Float64Array(s.c.length);
  const layerGains=new Float64Array(s.iceSlabs);
  let fullCellsAdded=0,partialCells=0,newTrapped=0;
  for(const e of events){
    const k=e.toIceIndex;
    if(!inEligible(s,k))throw new Error('Cannot nucleate ice inside gas-only halo');
    if(e.massKg<0||!Number.isFinite(e.massKg))throw new Error('Invalid condensed mass');
    s.iceMassAtCellKg[k]+=e.massKg;
    layerGains[e.iceLayer]+=e.massKg;
  }
  for(let k=0;k<s.c.length;k++){
    if(!inEligible(s,k)||!s.hexMask[k%s.nxy])continue;
    const target=s.iceMassAtCellKg[k]/s.cellIceCapacityKg;
    const next=Math.min(1,target);
    if(next+1e-12<prev[k])throw new Error('Ice occupancy spontaneously decreased');
    s.iceFraction[k]=next;
  }
  for(let k=0;k<s.c.length;k++){
    const delta=s.iceFraction[k]-prev[k];
    if(delta<=0)continue;
    // Gas pressure is not instantaneously raised when solid expands:
    // transfer the replaced volume's vapor to nearby exposed gas volumes.
    const displaced=Math.max(0,s.c[k]*s.volume*delta);
    const candidates=neighbors(s,k).filter(n=>!s.solid[n]&&s.iceFraction[n]<1);
    const weight=candidates.reduce((x,n)=>x+1-s.iceFraction[n],0);
    dMass[k]-=displaced;
    if(weight>0)for(const n of candidates)
      dMass[n]+=displaced*(1-s.iceFraction[n])/weight;
    else {s.trappedVaporKg+=displaced;newTrapped+=displaced;}
    if(s.iceFraction[k]>=1)fullCellsAdded++;
    else partialCells++;
  }
  for(let z=0;z<s.nz;z++)for(let p=0;p<s.nxy;p++){
    if(!s.hexMask[p])continue;
    const k=z*s.nxy+p;
    if(s.solid[k]&&prev[k]===1)continue;
    const initialMass=s.c[k]*s.volume*(1-prev[k]);
    const resulting=initialMass+dMass[k];
    if(resulting<-1e-25)throw new Error('Vapor displacement generated negative mass');
    const free=1-s.iceFraction[k];
    if(free<=0){
      if(resulting>1e-25)throw new Error('New solid retains unexplained gaseous mass');
      s.solid[k]=1;s.c[k]=0;
    }else{
      s.c[k]=Math.max(0,resulting)/(s.volume*free);
    }
    if(!Number.isFinite(s.c[k]))throw new Error('Vapor displacement is non-finite');
  }
  return {layerGainsKg:Array.from(layerGains),fullCellsAdded,partialCells,
    trappedVaporKg:newTrapped};
}
export function stepFiveLayerCrystal(s,dt,{tolerance=1e-8,maxIterations=1600}={}){
  if(s.model!=='snowflake/model-b-five-independent-slabs-v1')
    throw new TypeError('Expected independently moving five-layer crystal state');
  const before=layeredImplicitTrackedWaterKg(s);
  const startIceMass=fiveLayerIceMassKg(s);
  const report=stepLayeredImplicit(s,dt,{tolerance,maxIterations});
  const growth=settleFiveLayerFront(s,report.depositionEvents);
  // Solve thermodynamics only after the mass-conserving phase movement:
  // updated cell heat capacities use the NEW local ice fraction while
  // the prior enthalpy is carried unchanged through that volume change.
  const heat=s.thermal?stepLayeredThermal(s,dt,report.depositionEvents):null;
  const after=layeredImplicitTrackedWaterKg(s);
  const expected=report.boundaryMassInKg-report.boundaryMassOutKg;
  const difference=(after-before)-expected;
  if(!closeTo(after-before,expected))
    throw new Error('Model B water-budget error after voxel advancement: '+difference);
  const addedIce=fiveLayerIceMassKg(s)-startIceMass;
  if(!closeTo(addedIce,report.depositedKg))
    throw new Error('Front ice-mass increment differs from vapor uptake');
  if(!closeTo(fiveLayerIceMassKg(s)-s.seedIceMassKg,s.totalDepositedKg))
    throw new Error('Total grown ice mass differs from integrated deposition');
  return {...report,...growth,
    independentLayerMaps:true,
    layerDiagnostics:fiveLayerDiagnostics(s),
    postGrowthBudgetResidualKg:difference,
    heatTransportSolved:Boolean(heat),thermal:heat,
    iceVolumeClosure:'coarse voxel fraction + displaced-vapor exchange'};
}
