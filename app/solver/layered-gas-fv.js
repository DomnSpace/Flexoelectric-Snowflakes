/**
 * Step 2/15 — minimal CONSERVATIVE finite-volume vapor exchange on a
 * triangular basal lattice extruded into z-slabs.
 *
 * This is a separate opt-in numerical test kernel. It never edits or runs
 * hex25d.js, does not yet move the ice/gas boundary, resolve heat, or draw
 * five different snowflake layers. Ice is a FIXED cell mask; deposits are
 * separately inventoried for exact flux conservation.
 *
 * Each gas cell: c [kg/m^3], V=Ahex*dz [m^3], m=cV [kg].
 * Face transfer: D * faceArea/centreDistance * (cA-cB)*dt [kg].
 * Every transfer is posted once and subtracted/added on BOTH sides.
 */
import {diffusivityWaterVaporAir} from './thermo.js';
import {saturationMassDensityIce,sublimationEnthalpy} from './postseed.js';

const SQRT3=Math.sqrt(3);
const LATERAL=Object.freeze([[1,0],[0,1],[1,-1],[-1,0],[0,-1],[-1,1]]);
const positive=(v,name)=>{if(!Number.isFinite(v)||v<=0)throw new RangeError(name+' must be finite and positive');};
const nonnegative=(v,name)=>{if(!Number.isFinite(v)||v<0)throw new RangeError(name+' must be finite and nonnegative');};
export const FIVE_ICE_SLABS=5;
export const DEFAULT_HALO_SLABS=2;
export const phaseCellArea=dx=>SQRT3/2*dx*dx;
export const axialHexRadius=(q,r)=>Math.max(Math.abs(q),Math.abs(r),Math.abs(q+r));
const sumArray=a=>{let x=0;for(const v of a)x+=v;return x;};

function cellId(state,q,r,z){
  const i=q+state.R,j=r+state.R;
  if(i<0||i>=state.N||j<0||j>=state.N||z<0||z>=state.nz)return -1;
  if(!state.hexMask[j*state.N+i])return -1;
  return z*state.nxy+j*state.N+i;
}
export function layeredGasIndex(state,q,r,z){return cellId(state,q,r,z);}

function mode(x,name){
  if(!['sealed','reservoir'].includes(x))throw new RangeError(name+' boundary must be sealed or reservoir');
  return x;
}
export function createLayeredGasFV({
  N=21,dx=.75e-6,dz=.75e-6,iceSlabs=FIVE_ICE_SLABS,
  gasHaloSlabs=DEFAULT_HALO_SLABS,Tinf=258.15,pressurePa=101325,
  sigmaInf=.08,initialVaporKgM3=null,
  diffusivityM2s=null,prismAttachmentMPerS=0,basalAttachmentMPerS=0,
  surfaceEquilibriumKgM3=null,
  iceCells=[],
  boundaries={lateral:'sealed',top:'sealed',bottom:'sealed'}
}={}){
  if(!Number.isInteger(N)||N<7||N>97||N%2!==1)
    throw new RangeError('Step-2 FV prototype supports odd basal N=7..97 only');
  if(!Number.isInteger(iceSlabs)||iceSlabs<1||iceSlabs>9||
      !Number.isInteger(gasHaloSlabs)||gasHaloSlabs<1||gasHaloSlabs>5)
    throw new RangeError('Invalid number of ice-eligible / gas halo slabs');
  positive(dx,'dx');positive(dz,'dz');positive(Tinf,'temperature');
  positive(pressurePa,'pressure');nonnegative(sigmaInf,'ice supersaturation');
  nonnegative(prismAttachmentMPerS,'prism attachment rate');
  nonnegative(basalAttachmentMPerS,'basal attachment rate');
  const D=diffusivityM2s??diffusivityWaterVaporAir(Tinf,pressurePa);
  positive(D,'vapor diffusivity');
  const atmosphericIceSaturation=saturationMassDensityIce(Tinf);
  const ceq=surfaceEquilibriumKgM3??atmosphericIceSaturation;
  // A changed local attachment equilibrium does NOT change the far-field
  // gas inventory. It is an interface boundary condition, not the cloud
  // source concentration or a hidden overwrite of initial vapor mass.
  const c0=initialVaporKgM3??atmosphericIceSaturation*(1+sigmaInf);
  nonnegative(ceq,'surface equilibrium concentration');
  nonnegative(c0,'initial gas concentration');
  const bc={
    lateral:mode(boundaries.lateral??'sealed','lateral'),
    top:mode(boundaries.top??'sealed','top'),
    bottom:mode(boundaries.bottom??'sealed','bottom'),
    lateralConcentrationKgM3:boundaries.lateralConcentrationKgM3??c0,
    topConcentrationKgM3:boundaries.topConcentrationKgM3??c0,
    bottomConcentrationKgM3:boundaries.bottomConcentrationKgM3??c0
  };
  for(const side of ['lateral','top','bottom'])
    nonnegative(bc[side+'ConcentrationKgM3'],side+' reservoir concentration');
  const R=(N-1)/2,nxy=N*N,nz=iceSlabs+2*gasHaloSlabs,total=nxy*nz;
  const hexMask=new Uint8Array(nxy),solid=new Uint8Array(total);
  let basalSites=0;
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    if(axialHexRadius(i-R,j-R)<=R){hexMask[j*N+i]=1;basalSites++;}
  }
  const area=phaseCellArea(dx),volume=area*dz;
  const state={
    schema:'snowflake/vertical-gas-fv-step2',N,R,nxy,nz,dx,dz,iceSlabs,gasHaloSlabs,
    basalSites,hexMask,solid,
    c:new Float64Array(total),iceDepositedKg:new Float64Array(total),
    iceDepositionPerSlabKg:new Float64Array(nz),latentReleasedPerSlabJ:new Float64Array(nz),
    area,volume,sideArea:dx/SQRT3*dz,D,ceq,Tinf,pressurePa,
    kPrism:prismAttachmentMPerS,kBasal:basalAttachmentMPerS,
    boundaries:bc,Ls:sublimationEnthalpy(Tinf,pressurePa).Ls,
    timeSeconds:0,externalMassInKg:0,externalMassOutKg:0,
    totalDepositedKg:0,latentReleasedJ:0,
    lastStep:null,heatTransportSolved:false,interfaceMotionSolved:false
  };
  for(let z=0;z<nz;z++)for(let k=0;k<nxy;k++)if(hexMask[k])
    state.c[z*nxy+k]=c0;
  for(const p of iceCells){
    if(!Number.isInteger(p.q)||!Number.isInteger(p.r)||!Number.isInteger(p.z)||
      p.z<gasHaloSlabs||p.z>=gasHaloSlabs+iceSlabs)
      throw new RangeError('Ice sites must be integer axial positions inside five ice-eligible slabs');
    const k=cellId(state,p.q,p.r,p.z);
    if(k<0)throw new RangeError('Ice site outside hexagonal mesh');
    if(solid[k])throw new RangeError('Duplicate fixed ice site');
    solid[k]=1;state.c[k]=0;
  }
  return state;
}
function isGas(state,k){return k>=0&&!state.solid[k];}
export function layeredGasMassKg(state){
  let mass=0;for(let z=0;z<state.nz;z++){
    const off=z*state.nxy;
    for(let k=0;k<state.nxy;k++)if(state.hexMask[k]&&!state.solid[off+k])
      mass+=state.c[off+k]*state.volume;
  }
  return mass;
}
export function layeredTrackedWaterKg(state){
  return layeredGasMassKg(state)+state.totalDepositedKg;
}
/**
 * Conservative positive-flux explicit timestep recommendation, based on
 * worst-case gas/gas, reservoir half-cell and ice attachment conductances.
 * Taking <=this value avoids excessive transport CFL at small N; additional
 * inventory checking still detects nonlinear oversubscription.
 */
export function layeredExplicitDtLimit(state){
  const gs=state.D*state.sideArea/state.dx;
  const gz=state.D*state.area/state.dz;
  const ks=state.kPrism*state.sideArea,kz=state.kBasal*state.area;
  const maxConductance=12*gs+4*gz+6*ks+2*kz;
  return .45*state.volume/maxConductance;
}
/**
 * One physical step:
 *  1) evaluate ALL face requests from the unmodified previous concentration;
 *  2) compute total outward demand per gas volume (including all 5 ice slabs);
 *  3) either reject oversubscribed dt or scale all competing sinks together;
 *  4) apply opposite debit/credit for every internal flux, log boundary flow.
 *
 * No deposition is possible on buried solid-solid faces; the vertical
 * exposed orientation picks kBasal, lateral faces pick kPrism.
 */
export function stepLayeredGasFV(state,dt,{allowFluxLimiter=false}={}){
  positive(dt,'dt');
  const limit=layeredExplicitDtLimit(state);
  if(!allowFluxLimiter&&dt>limit*(1+1e-10))
    throw new RangeError('dt exceeds explicitly stable FV flux limit; use substeps / an implicit transport solver');
  const {N,R,nxy,nz,c,solid,hexMask,volume,D,area,sideArea,dx,dz}=state;
  const demands=new Float64Array(c.length),delta=new Float64Array(c.length),faces=[];
  const gas=k=>k>=0&&!solid[k];
  function add(from,to,mass,kind='gas'){
    if(!(mass>0))return;
    faces.push({from,to,mass,kind});
    if(from>=0)demands[from]+=mass;
  }
  function connect(a,b,conductance){
    if(!gas(a)||!gas(b))return;
    const dm=conductance*(c[a]-c[b])*dt;
    if(dm>0)add(a,b,dm);else if(dm<0)add(b,a,-dm);
  }
  function iceFace(gasIndex,solidIndex,faceArea,attachmentRate){
    if(!gas(gasIndex)||!solid[solidIndex]||attachmentRate<=0)return;
    const dm=attachmentRate*faceArea*Math.max(0,c[gasIndex]-state.ceq)*dt;
    add(gasIndex,solidIndex,dm,'deposit');
  }
  function outside(gasIndex,faceArea,halfSeparation,reservoirConcentration){
    if(!gas(gasIndex))return;
    const dm=D*faceArea/halfSeparation*(c[gasIndex]-reservoirConcentration)*dt;
    if(dm>0)add(gasIndex,-1,dm,'out');
    else if(dm<0)add(-1,gasIndex,-dm,'in');
  }
  const gxy=D*sideArea/dx,gz=D*area/dz;
  for(let z=0;z<nz;z++){
    const zoff=z*nxy;
    for(let j=0;j<N;j++)for(let i=0;i<N;i++){
      const base=j*N+i;if(!hexMask[base])continue;
      const k=zoff+base;
      if(!gas(k))continue;
      for(let direction=0;direction<6;direction++){
        const [dq,dr]=LATERAL[direction],ii=i+dq,jj=j+dr;
        const adj=(ii>=0&&ii<N&&jj>=0&&jj<N&&hexMask[jj*N+ii])?
          zoff+jj*N+ii:-1;
        if(adj<0){
          if(state.boundaries.lateral==='reservoir')
            outside(k,sideArea,dx/2,state.boundaries.lateralConcentrationKgM3);
        }else if(solid[adj]){
          iceFace(k,adj,sideArea,state.kPrism);
        }else if(direction<3)connect(k,adj,gxy);
      }
      for(const dd of [-1,1]){
        const zz=z+dd;
        if(zz<0||zz>=nz){
          const side=zz<0?'bottom':'top';
          if(state.boundaries[side]==='reservoir')
            outside(k,area,dz/2,state.boundaries[side+'ConcentrationKgM3']);
        }else{
          const adj=zz*nxy+base;
          if(solid[adj])iceFace(k,adj,area,state.kBasal);
          else if(dd>0)connect(k,adj,gz);
        }
      }
    }
  }
  let limitedFaces=0,limiterMinimum=1;
  const fraction=new Float64Array(c.length);
  for(let z=0;z<nz;z++)for(let p=0;p<nxy;p++){
    if(!hexMask[p])continue;
    const k=z*nxy+p;if(solid[k])continue;
    const available=c[k]*volume;
    const ratio=demands[k]>0?Math.min(1,available/demands[k]):1;
    fraction[k]=ratio;
    if(ratio<1-1e-12){limitedFaces++;limiterMinimum=Math.min(limiterMinimum,ratio);}
  }
  if(limitedFaces&&!allowFluxLimiter)
    throw new RangeError('Combined vapor demands exceed a gas-cell inventory; no state changed');
  let deposited=0,boundaryIn=0,boundaryOut=0,internalMoved=0;
  const layerDeposit=new Float64Array(nz);
  for(const f of faces){
    const amount=f.mass*(f.from<0?1:fraction[f.from]);
    if(f.from>=0)delta[f.from]-=amount/volume;
    if(f.kind==='gas'){delta[f.to]+=amount/volume;internalMoved+=amount;}
    else if(f.kind==='deposit'){
      const iceZ=Math.floor(f.to/nxy);
      state.iceDepositedKg[f.to]+=amount;
      state.iceDepositionPerSlabKg[iceZ]+=amount;
      layerDeposit[iceZ]+=amount;deposited+=amount;
      const heat=amount*state.Ls;
      state.latentReleasedPerSlabJ[iceZ]+=heat;
    }else if(f.kind==='out')boundaryOut+=amount;
    else if(f.kind==='in'){delta[f.to]+=amount/volume;boundaryIn+=amount;}
  }
  for(let z=0;z<nz;z++)for(let p=0;p<nxy;p++){
    if(!hexMask[p])continue;
    const k=z*nxy+p;if(solid[k])continue;
    c[k]+=delta[k];
    if(c[k]<0&&c[k]>-1e-12)c[k]=0;
    if(!Number.isFinite(c[k])||c[k]<0)throw new Error('Negative/non-finite vapor concentration');
  }
  state.totalDepositedKg+=deposited;
  state.latentReleasedJ+=deposited*state.Ls;
  state.externalMassInKg+=boundaryIn;
  state.externalMassOutKg+=boundaryOut;
  state.timeSeconds+=dt;
  const report={
    dtSeconds:dt,explicitLimitSeconds:limit,
    candidateFaces:faces.length,limitedCells:limitedFaces,limiterMinimum,
    internalMassTransferredKg:internalMoved,
    depositedKg:deposited,perIceSlabDepositedKg:Array.from(
      layerDeposit.slice(state.gasHaloSlabs,state.gasHaloSlabs+state.iceSlabs)),
    boundaryMassInKg:boundaryIn,boundaryMassOutKg:boundaryOut,
    gasMassKg:layeredGasMassKg(state),iceDepositTotalKg:state.totalDepositedKg,
    latentReleasedJ:deposited*state.Ls,
    heatTransportSolved:false,iceSurfaceMotionSolved:false
  };
  state.lastStep=report;
  return report;
}
