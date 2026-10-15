/**
 * Experimental 2.5-D ice growth on a triangular computational lattice.
 *
 * Six equidistant neighbor directions (60° apart), Voronoi control cells
 * shaped as regular hexagons. These are MICROMETRE transport volumes, NOT
 * individual hydrogen atoms or a molecular model of ice Ih.
 *
 * Axial (q,r) coordinates:
 *    x = a (q + r/2),  y = a sqrt(3) r/2
 * with six neighbors ±(1,0), ±(0,1), ±(1,-1).
 *
 * The six-neighbor Laplacian for isotropic diffusion is
 *    ∇²u ≈ [2/(3a²)] (Σ₆ u_neighbor - 6u).
 *
 * Hexagonal Dirichlet domain, NOT Cartesian image padding. No free charge,
 * elastic strain, vertical PDE, electrostatics, airflow, or calibrated SDAK.
 * Existing Cartesian solver remains independent for A/B comparison.
 */
import {diffusivityWaterVaporAir,thermalConductivityAir} from './thermo.js';
import {iceIhProperties} from './ice-ih-eos.js';
import {attachmentCoefficient,iceKineticVelocity,saturationMassDensityIce,sublimationEnthalpy} from './postseed.js';
import {DEFAULT_ICE_KINETICS,capillaryLengthIceVapor} from './growth25d.js';

const PI=Math.PI,RT3=Math.sqrt(3),SQRT3_OVER_2=RT3/2;
const DIR=Object.freeze([[1,0],[-1,0],[0,1],[0,-1],[1,-1],[-1,1]]);
const VECT=DIR.map(([q,r])=>[q+r/2,SQRT3_OVER_2*r]);
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const id=(i,j,N)=>j*N+i;
const valid=(i,j,N)=>i>=0&&j>=0&&i<N&&j<N;
export const HEX_NEIGHBORS=DIR;
export const HEX_CELL_AREA_MULTIPLIER=SQRT3_OVER_2;

export function axialToXY(q,r,a){return {x:a*(q+r/2),y:a*SQRT3_OVER_2*r};}
export function axialRadius(q,r){return Math.max(Math.abs(q),Math.abs(r),Math.abs(q+r));}
export function axialRoundFromXY(x,y,a){
  const q=x/a-y/(a*RT3),r=2*y/(a*RT3),s=-q-r;
  let Q=Math.round(q),R=Math.round(r),S=Math.round(s);
  const dq=Math.abs(Q-q),dr=Math.abs(R-r),ds=Math.abs(S-s);
  if(dq>dr&&dq>ds)Q=-R-S;
  else if(dr>ds)R=-Q-S;
  else S=-Q-R;
  return {q:Q,r:R};
}

export function makeHexGrowthState({
  N=73,dx=.75e-6,seedRadius=1.5e-6,Tinf=258.15,
  sigmaInf=.08,pressurePa=101325,basalHalfThickness=1.5e-6,
  kinetics=DEFAULT_ICE_KINETICS
}={}){
  if(!Number.isInteger(N)||N<13||N%2!==1)throw new RangeError('Hex grid N must be odd and >= 13');
  if(!(dx>0&&seedRadius>=0&&sigmaInf>=0&&Tinf>110&&pressurePa>0))throw new RangeError('Invalid hex physical grid/environment');
  const n=N*N,R=(N-1)/2,mask=new Uint8Array(n),ice=new Uint8Array(n);
  const phi=new Float64Array(n),sigma=new Float64Array(n),dT=new Float64Array(n),heatSource=new Float64Array(n);
  sigma.fill(sigmaInf);
  let activeSites=0;
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    const q=i-R,r=j-R,k=id(i,j,N);if(axialRadius(q,r)>R)continue;
    mask[k]=1;activeSites++;
    const {x,y}=axialToXY(q,r,dx);
    if(Math.hypot(x,y)<=seedRadius){ice[k]=1;phi[k]=1;sigma[k]=0;}
  }
  return {topology:'triangular',N,dx,R,mask,activeSites,
    ice,phi,sigma,dT,heatSource,Tinf,sigmaInf,pressurePa,
    basalHalfThickness,kinetics:JSON.parse(JSON.stringify(kinetics)),
    t:0,totalDepositedMass:0,totalLatentHeat:0,lastBudget:null,
    electric:null,domainExpansions:0};
}

export function hexCellAtXY(state,x,y){
  const a=axialRoundFromXY(x,y,state.dx);
  const i=a.q+state.R,j=a.r+state.R;
  if(!valid(i,j,state.N)||!state.mask[id(i,j,state.N)])return -1;
  return id(i,j,state.N);
}

function occupied(state,k){return state.ice[k]?1:state.phi[k];}
function forNeighbor(state,i,j,fn){
  for(let d=0;d<6;d++){
    const ii=i+DIR[d][0],jj=j+DIR[d][1];
    if(!valid(ii,jj,state.N))continue;
    const k=id(ii,jj,state.N);
    if(state.mask[k])fn(k,d);
  }
}
function atBoundary(state,i,j){return axialRadius(i-state.R,j-state.R)>=state.R;}

function phaseGradient(state,i,j){
  let gx=0,gy=0;let count=0;
  forNeighbor(state,i,j,(k,d)=>{
    const f=occupied(state,k);
    gx+=f*VECT[d][0];gy+=f*VECT[d][1];count++;
  });
  // centered six-site gradient; regular interior cells have 6 neighbors.
  return {gx:gx/(3*state.dx),gy:gy/(3*state.dx),count};
}
function outwardNormal(state,i,j){
  const grad=phaseGradient(state,i,j),g=Math.hypot(grad.gx,grad.gy);
  if(g>1e-10)return {nx:-grad.gx/g,ny:-grad.gy/g,g};
  const p=axialToXY(i-state.R,j-state.R,state.dx);
  const radial=Math.hypot(p.x,p.y);
  return radial>0?{nx:p.x/radial,ny:p.y/radial,g}: {nx:1,ny:0,g};
}

export function hexInterfaceGeometry(state,i,j){
  const normal=outwardNormal(state,i,j);
  let div=0;
  forNeighbor(state,i,j,(k,d)=>{
    const nx=outwardNormal(state,k%state.N,Math.floor(k/state.N));
    div+=nx.nx*VECT[d][0]+nx.ny*VECT[d][1];
  });
  const curvature=clamp(div/(3*state.dx),-2/state.dx,2/state.dx);
  return {...normal,theta:Math.atan2(normal.ny,normal.nx),curvature};
}

export function hexIsInterface(state,i,j){
  const k=id(i,j,state.N);
  if(!state.mask[k]||state.ice[k]||state.phi[k]>=.999||atBoundary(state,i,j))return false;
  let found=false;
  forNeighbor(state,i,j,k2=>{if(occupied(state,k2)>.15)found=true;});
  return found;
}

export function relaxHexScalar(state,field='sigma',iterations=40){
  const isVapor=field==='sigma';
  if(!isVapor&&field!=='dT')throw new Error('Unknown hex transport scalar');
  const {N,dx,Tinf,sigmaInf,mask,ice,heatSource}=state;
  const cref=isVapor?saturationMassDensityIce(Tinf):0;
  const kAir=isVapor?0:thermalConductivityAir(Tinf);
  let a=new Float64Array(state[field]),b=new Float64Array(a.length);
  for(let iter=0;iter<iterations;iter++){
    for(let j=0;j<N;j++)for(let i=0;i<N;i++){
      const k=id(i,j,N);
      if(!mask[k]){b[k]=isVapor?sigmaInf:0;continue;}
      if(atBoundary(state,i,j)){b[k]=isVapor?sigmaInf:0;continue;}
      if(isVapor&&ice[k]){
        b[k]=saturationMassDensityIce(Tinf+state.dT[k])/cref-1;continue;
      }
      let sum=0,n=0;
      forNeighbor(state,i,j,k2=>{sum+=a[k2];n++;});
      // All interior sites of the regular hex domain have six neighbors.
      if(n!==6)throw new Error('Broken six-neighbor hex interior at '+i+','+j);
      b[k]=sum/6+(isVapor?0:heatSource[k]*dx*dx/(4*kAir));
    }
    const temp=a;a=b;b=temp;
  }
  state[field].set(a);
  return state;
}

export function resolvedHexGrowthStep(state,dt,{vaporIterations=45,heatIterations=55}={}){
  if(!(dt>0&&Number.isFinite(dt)))throw new RangeError('Invalid dt');
  relaxHexScalar(state,'dT',heatIterations);
  relaxHexScalar(state,'sigma',vaporIterations);
  const {N,dx,Tinf,pressurePa,ice,phi,sigma,heatSource}=state;
  const rhoi=iceIhProperties(Tinf,Math.min(pressurePa,210e6)).rho;
  const cRef=saturationMassDensityIce(Tinf);
  const Dv=diffusivityWaterVaporAir(Tinf,pressurePa);
  const vkin=iceKineticVelocity(Tinf,rhoi);
  const Ls=sublimationEnthalpy(Tinf,pressurePa).Ls;
  const kParams=state.kinetics;
  const cellArea=HEX_CELL_AREA_MULTIPLIER*dx*dx;
  const faceLength=dx/RT3;
  const depth=2*state.basalHalfThickness;
  const nextHeat=new Float64Array(heatSource.length),updates=[];
  let totalMassRate=0,deposit=0,alphaSum=0,sigSum=0,curvatureSum=0,capillarySum=0,count=0;
  for(let j=1;j<N-1;j++)for(let i=1;i<N-1;i++)if(hexIsInterface(state,i,j)){
    const k=id(i,j,N);
    const g=hexInterfaceGeometry(state,i,j);
    let iceLinks=0,localTemp=0,iceWeight=0,filled=0;
    forNeighbor(state,i,j,k2=>{
      const occ=occupied(state,k2);
      if(occ>.15){iceLinks+=occ;filled++;localTemp+=occ*state.dT[k2];iceWeight+=occ;}
    });
    const Ts=Tinf+(iceWeight?localTemp/iceWeight:state.dT[k]);
    const cSat=saturationMassDensityIce(Ts);
    const delta=capillaryLengthIceVapor(Ts,rhoi,kParams.gammaIv);
    const exponent=kParams.capillarity===false?0:delta*g.curvature;
    const cEq=cSat*Math.exp(exponent);
    const cGas=cRef*Math.max(1e-9,1+sigma[k]);
    const sigmaDrive=Math.max(0,cGas/cEq-1);
    let model={...kParams.prism};
    if(kParams.sdak)model.sigma0*=Math.max(.08,1-(kParams.sdakStrength??.55)*(1-filled/6));
    const base=attachmentCoefficient(sigmaDrive,model);
    const alpha=clamp(base*(1+(kParams.anisotropy6??0)*Math.cos(6*g.theta)),0,1);
    const Rdiff=dx/(2*Math.max(1e-30,Dv));
    const Rkin=cEq/Math.max(1e-30,rhoi*alpha*vkin);
    const flux=Math.max(0,cGas-cEq)/(Rdiff+Rkin);
    const area=iceLinks*faceLength*depth;
    const rate=flux*area;
    const capacity=rhoi*cellArea*depth;
    const dphi=rate*dt/Math.max(1e-30,capacity);
    updates.push([k,dphi]);
    nextHeat[k]+=rate*Ls/Math.max(1e-30,cellArea*depth);
    totalMassRate+=rate;alphaSum+=alpha;sigSum+=sigmaDrive;
    curvatureSum+=g.curvature;capillarySum+=exponent;count++;
  }
  for(const [k,added] of updates){
    const prior=phi[k];phi[k]=Math.min(1,prior+added);
    if(phi[k]>=.999)ice[k]=1;
    deposit+=(phi[k]-prior)*rhoi*cellArea*depth;
  }
  heatSource.set(nextHeat);
  const meanSigma=count?sigSum/count:0;
  const alphaBasal=attachmentCoefficient(meanSigma,kParams.basal);
  const vBasal=alphaBasal*vkin*meanSigma;
  state.basalHalfThickness+=vBasal*dt;
  state.totalDepositedMass+=deposit;
  state.totalLatentHeat+=deposit*Ls;
  state.t+=dt;
  let maxDeltaT=0;for(let k=0;k<state.mask.length;k++)if(state.mask[k])maxDeltaT=Math.max(maxDeltaT,state.dT[k]);
  state.lastBudget={
    totalMassRate,deposited:deposit,Ls,rhoIce:rhoi,Dv,cRef,vkin,meanSigma,
    meanPrismAlpha:count?alphaSum/count:0,alphaBasal,vBasal,
    interfaceCells:count,maxDeltaT,
    meanCurvature:count?curvatureSum/count:0,
    meanCapillaryExponent:count?capillarySum/count:0,
    meanElectricLogShift:0,maxElectricLogShift:0,
    maxElectricFieldVM:0,maxFlexoPolarizationCm2:0,electricResidualRmsV:0,
    stencil:'six-equidistant-triangular',orientation:'six-site-phase-gradient',
    lattice:'mesoscopic transport mesh, NOT H-bond positions'
  };
  return state;
}

export function hexGrowthMorphology(state,rays=180){
  const radii=[],maxR=state.R*state.dx,step=state.dx*.3;
  let area=0;
  for(let k=0;k<state.mask.length;k++)if(state.mask[k])area+=occupied(state,k);
  for(let q=0;q<rays;q++){
    const theta=2*PI*q/rays;let radius=0;
    for(let r=0;r<=maxR;r+=step){
      const k=hexCellAtXY(state,r*Math.cos(theta),r*Math.sin(theta));
      if(k<0)break;
      if(occupied(state,k)>.15)radius=r;
    }
    radii.push(radius);
  }
  const meanRadius=radii.reduce((a,b)=>a+b,0)/rays,modes={};
  for(const mode of [1,2,3,4,6,8,12]){
    let a=0,b=0;
    for(let q=0;q<rays;q++){
      const theta=2*PI*q/rays,d=radii[q]-meanRadius;
      a+=d*Math.cos(mode*theta);b+=d*Math.sin(mode*theta);
    }
    a*=2/rays;b*=2/rays;
    modes[mode]={a,b,amplitude:Math.hypot(a,b),
      relative:meanRadius?Math.hypot(a,b)/meanRadius:0,phase:Math.atan2(-b,a)/mode};
  }
  const eqRadius=Math.sqrt(area*HEX_CELL_AREA_MULTIPLIER*state.dx*state.dx/PI);
  return {radius:radii,meanRadius,modes,m6Amplitude:modes[6].amplitude,
    m6Phase:modes[6].phase,equivalentBasalRadius:eqRadius,
    halfThickness:state.basalHalfThickness,
    aspectRatio:eqRadius/Math.max(1e-30,state.basalHalfThickness)};
}

export function auditHexDomain(state,{warningCells=12,stopCells=6,plateRatioWarning=.65,plateRatioStop=1.25}={}){
  let halo=state.R,tip=0,count=0;
  for(let j=0;j<state.N;j++)for(let i=0;i<state.N;i++){
    const k=id(i,j,state.N);
    if(!state.mask[k]||occupied(state,k)<.15)continue;
    const q=i-state.R,r=j-state.R;
    halo=Math.min(halo,state.R-axialRadius(q,r));
    const p=axialToXY(q,r,state.dx);tip=Math.max(tip,Math.hypot(p.x,p.y));
    count++;
  }
  const ratio=state.basalHalfThickness/Math.max(state.dx,tip);
  const stopHalo=halo<=stopCells;
  const plateWarning=ratio>plateRatioWarning;
  const plateStop=ratio>plateRatioStop&&state.t>5;
  const stop=stopHalo||plateStop;
  return {
    N:state.N,dx:state.dx,domainUm:(state.N-1)*state.dx*1e6,
    occupied:count,haloCells:halo,haloUm:halo*state.dx*1e6,
    tipRadiusUm:tip*1e6,thicknessUm:2*state.basalHalfThickness*1e6,
    plateRatio:ratio,warning:halo<=warningCells||plateWarning,
    stop,stopHalo,plateWarning,plateStop,
    reason:stopHalo?'HEXAGONAL FAR-FIELD BOUNDARY TOO CLOSE TO GROWTH FRONT':
      plateStop?'SHARED-THICKNESS PLATE MODEL INVALID: EXCESSIVE VERTICAL GROWTH':
      halo<=warningCells?'HEXAGONAL DIRICHLET BOUNDARY MAY SELECT M6':
      plateWarning?'SHARED-THICKNESS PLATE APPROXIMATION STRETCHED':
      'WITHIN EXPERIMENTAL HEX-MESH LIMITS'
  };
}

/**
 * Extend the finite transport domain without restarting the crystal.
 * Conserves ice, phase, heat-source arrays and cumulative mass/time.
 * Refilled outer vapor/thermal regions are relaxation initial guesses,
 * not an exact infinite-domain Green's function.
 */
export function enlargeHexDomain(state,newN){
  if(newN<=state.N||newN%2!==1)throw new RangeError('Hex grid must expand by an odd positive size');
  const next=makeHexGrowthState({
    N:newN,dx:state.dx,seedRadius:0,Tinf:state.Tinf,
    sigmaInf:state.sigmaInf,pressurePa:state.pressurePa,
    basalHalfThickness:state.basalHalfThickness,kinetics:state.kinetics
  });
  const offset=(newN-state.N)/2;
  for(let j=0;j<state.N;j++)for(let i=0;i<state.N;i++){
    const old=id(i,j,state.N);
    if(!state.mask[old])continue;
    const k=id(i+offset,j+offset,newN);
    next.ice[k]=state.ice[old];next.phi[k]=state.phi[old];
    next.sigma[k]=state.sigma[old];next.dT[k]=state.dT[old];
    next.heatSource[k]=state.heatSource[old];
  }
  next.t=state.t;
  next.totalDepositedMass=state.totalDepositedMass;
  next.totalLatentHeat=state.totalLatentHeat;
  next.lastBudget=state.lastBudget?{...state.lastBudget}:null;
  next.domainExpansions=(state.domainExpansions||0)+1;
  // Solve a new field in the enlarged gas domain before the next interface step.
  relaxHexScalar(next,'dT',100);
  relaxHexScalar(next,'sigma',100);
  return next;
}
