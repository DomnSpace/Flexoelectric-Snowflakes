import {
  CONSTANTS,
  diffusivityWaterVaporAir,
  saturationVaporPressureIceMK,
  thermalConductivityAir,
} from './thermo.js';
import { iceIhProperties } from './ice-ih-eos.js';
import { attachmentCoefficient, iceKineticVelocity, saturationMassDensityIce, sublimationEnthalpy } from './postseed.js';

/**
 * Milestone 2.5 resolved growth solver.
 *
 * sigma[] is a reference-normalized vapor-density excess
 *   u = c_v / c_sat,i(T_inf) - 1,
 * so ∇²u = 0 is the constant-D quasi-steady diffusion equation. Local physical
 * supersaturation is reconstructed against c_sat,i(T_inf + ΔT), including a
 * 2.5-D Gibbs-Thomson curvature correction at the interface.
 *
 * The scalar fields use a 9-point isotropic Laplacian to suppress artificial m=4
 * selection from the Cartesian mesh. Interface orientation comes from the local
 * phase-fraction gradient, never from radial angle about the domain centre.
 */

const dirs4=[[1,0],[-1,0],[0,1],[0,-1]];
const dirs8=[[1,0],[-1,0],[0,1],[0,-1],[1,1],[1,-1],[-1,1],[-1,-1]];
const id=(i,j,N)=>j*N+i;
const inb=(i,j,N)=>i>=0&&j>=0&&i<N&&j<N;
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));

export const ICE_VAPOR_CAPILLARITY=Object.freeze({
  gammaIv:0.109,
  units:'J m^-2',
  status:'reference-value',
  source:'Libbrecht 2005, Rep. Prog. Phys. 68, 855–895',
});

export const DEFAULT_ICE_KINETICS=Object.freeze({
  prism:{mode:'nucleation',A:1,sigma0:0.003},
  basal:{mode:'nucleation',A:1,sigma0:0.01},
  anisotropy6:0.18,
  sdak:false,
  sdakStrength:0.55,
  capillarity:true,
  gammaIv:ICE_VAPOR_CAPILLARITY.gammaIv,
});

export function makeGrowth25DState({
  N=73,
  dx=0.75e-6,
  seedRadius=1.5e-6,
  Tinf=258.15,
  sigmaInf=0.05,
  pressurePa=101325,
  basalHalfThickness=1.5e-6,
  kinetics=DEFAULT_ICE_KINETICS,
}={}){
  if(N%2===0)N+=1;
  const n=N*N,ice=new Uint8Array(n),phi=new Float64Array(n),sigma=new Float64Array(n),dT=new Float64Array(n),heatSource=new Float64Array(n);
  sigma.fill(sigmaInf);
  const c=(N-1)/2;
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    const x=(i-c)*dx,y=(j-c)*dx,r=Math.hypot(x,y),k=id(i,j,N);
    if(r<=seedRadius){ice[k]=1;phi[k]=1;sigma[k]=0;}
  }
  return {N,dx,Tinf,sigmaInf,pressurePa,ice,phi,sigma,dT,heatSource,basalHalfThickness,
    kinetics:JSON.parse(JSON.stringify(kinetics)),t:0,totalDepositedMass:0,totalLatentHeat:0,lastBudget:null};
}

function occ(state,i,j){
  if(!inb(i,j,state.N))return 0;
  const k=id(i,j,state.N);return state.ice[k]?1:state.phi[k];
}
function isBoundary(i,j,N){return i===0||j===0||i===N-1||j===N-1;}

export function isInterfaceGas(state,i,j){
  const {N,ice}=state,k=id(i,j,N);if(!inb(i,j,N)||ice[k]||state.phi[k]>=0.999)return false;
  return dirs8.some(([di,dj])=>inb(i+di,j+dj,N)&&occ(state,i+di,j+dj)>0.15);
}

/** Outward interface normal and in-plane curvature from the phase fraction. */
export function interfaceGeometry(state,i,j){
  const {dx,N}=state;
  const E=occ(state,i+1,j),W=occ(state,i-1,j),S=occ(state,i,j+1),Nn=occ(state,i,j-1);
  const NE=occ(state,i+1,j-1),NW=occ(state,i-1,j-1),SE=occ(state,i+1,j+1),SW=occ(state,i-1,j+1),C=occ(state,i,j);
  // Sobel first derivatives reduce pixel-axis bias.
  const fx=((NE+2*E+SE)-(NW+2*W+SW))/(8*dx);
  const fy=((SW+2*S+SE)-(NW+2*Nn+NE))/(8*dx);
  const fxx=(E-2*C+W)/(dx*dx),fyy=(S-2*C+Nn)/(dx*dx);
  const fxy=(SE-SW-NE+NW)/(4*dx*dx);
  const g=Math.hypot(fx,fy);
  let nx=0,ny=0;
  if(g>1e-12){nx=-fx/g;ny=-fy/g;}
  else{
    const c=(N-1)/2,x=(i-c)*dx,y=(j-c)*dx,r=Math.hypot(x,y)||1;nx=x/r;ny=y/r;
  }
  let curvature=0;
  if(g>1e-12){
    const num=fxx*fy*fy-2*fx*fy*fxy+fyy*fx*fx;
    curvature=-num/Math.max(1e-30,g**3);
  }
  curvature=clamp(curvature,-2/dx,2/dx);
  const interfaceLength=clamp(g*dx*dx,0.35*dx,Math.SQRT2*dx);
  return {nx,ny,theta:Math.atan2(ny,nx),curvature,gradientMagnitude:g,interfaceLength};
}

export function prismAttachmentAt(state,i,j,sigmaDrive){
  const th=interfaceGeometry(state,i,j).theta,k=state.kinetics;
  let model={...k.prism};
  if(k.sdak){
    let occupied=0;
    for(const[di,dj]of dirs8)if(inb(i+di,j+dj,state.N))occupied+=occ(state,i+di,j+dj)>0.15?1:0;
    const edgeFactor=Math.max(0,1-occupied/8);
    model.sigma0*=Math.max(0.08,1-(k.sdakStrength??0.55)*edgeFactor);
  }
  const base=attachmentCoefficient(Math.max(0,sigmaDrive),model);
  const six=1+(k.anisotropy6??0)*Math.cos(6*th);
  return Math.max(0,Math.min(1,base*six));
}

export function capillaryLengthIceVapor(T,rhoIce,gammaIv=ICE_VAPOR_CAPILLARITY.gammaIv){
  return gammaIv*CONSTANTS.MW/(rhoIce*CONSTANTS.R*T);
}

export function localInterfaceEquilibrium(state,i,j,iceProps=null){
  const {Tinf,pressurePa}=state;
  const props=iceProps??iceIhProperties(Tinf,Math.min(pressurePa,210e6));
  const geom=interfaceGeometry(state,i,j);
  let tempSum=0,wSum=0;
  for(const[di,dj]of dirs8){
    const ii=i+di,jj=j+dj;if(!inb(ii,jj,state.N))continue;
    const w=occ(state,ii,jj);if(w>0.15){tempSum+=w*state.dT[id(ii,jj,state.N)];wSum+=w;}
  }
  const Ts=Tinf+(wSum?tempSum/wSum:state.dT[id(i,j,state.N)]);
  const cSatLocal=saturationMassDensityIce(Ts);
  const delta=capillaryLengthIceVapor(Ts,props.rho,state.kinetics.gammaIv??ICE_VAPOR_CAPILLARITY.gammaIv);
  const capExponent=state.kinetics.capillarity===false?0:delta*geom.curvature;
  const cEq=cSatLocal*Math.exp(capExponent);
  return {Ts,cEq,cSatLocal,capillaryLength:delta,capillaryExponent:capExponent,...geom};
}

export function localIceSupersaturation(state,i,j){
  const k=id(i,j,state.N),cRef=saturationMassDensityIce(state.Tinf),cGas=cRef*Math.max(1e-9,1+state.sigma[k]);
  const eq=localInterfaceEquilibrium(state,i,j);
  return {sigmaLocal:cGas/eq.cEq-1,cGas,...eq};
}

export function relaxVaporField(state,iterations=80){
  const {N,ice,sigma,sigmaInf,Tinf}=state;
  const a=new Float64Array(sigma),b=new Float64Array(sigma.length),cRef=saturationMassDensityIce(Tinf);
  for(let it=0;it<iterations;it++){
    for(let j=0;j<N;j++)for(let i=0;i<N;i++){
      const k=id(i,j,N);
      if(isBoundary(i,j,N)){b[k]=sigmaInf;continue;}
      if(ice[k]){
        const Ts=Tinf+state.dT[k];
        b[k]=saturationMassDensityIce(Ts)/cRef-1;
        continue;
      }
      const card=a[id(i+1,j,N)]+a[id(i-1,j,N)]+a[id(i,j+1,N)]+a[id(i,j-1,N)];
      const diag=a[id(i+1,j+1,N)]+a[id(i+1,j-1,N)]+a[id(i-1,j+1,N)]+a[id(i-1,j-1,N)];
      b[k]=(4*card+diag)/20;
    }
    a.set(b);
  }
  sigma.set(a);return state;
}

export function relaxHeatField(state,iterations=100){
  const {N,dT,heatSource,dx,Tinf}=state,kAir=thermalConductivityAir(Tinf);
  const a=new Float64Array(dT),b=new Float64Array(dT.length);
  for(let it=0;it<iterations;it++){
    for(let j=0;j<N;j++)for(let i=0;i<N;i++){
      const k=id(i,j,N);
      if(isBoundary(i,j,N)){b[k]=0;continue;}
      const card=a[id(i+1,j,N)]+a[id(i-1,j,N)]+a[id(i,j+1,N)]+a[id(i,j-1,N)];
      const diag=a[id(i+1,j+1,N)]+a[id(i+1,j-1,N)]+a[id(i-1,j+1,N)]+a[id(i-1,j-1,N)];
      b[k]=(4*card+diag+6*heatSource[k]*dx*dx/kAir)/20;
    }
    a.set(b);
  }
  dT.set(a);return state;
}

export function resolvedGrowthStep(state,dt,{vaporIterations=45,heatIterations=55}={}){
  const {N,dx,Tinf,pressurePa,phi,sigma,heatSource}=state;
  relaxHeatField(state,heatIterations);relaxVaporField(state,vaporIterations);
  const iceProps=iceIhProperties(Tinf,Math.min(pressurePa,210e6)),rhoIce=iceProps.rho;
  const Dv=diffusivityWaterVaporAir(Tinf,pressurePa),cRef=saturationMassDensityIce(Tinf),vkin=iceKineticVelocity(Tinf,rhoIce),Ls=sublimationEnthalpy(Tinf,pressurePa).Ls;
  const nextHeat=new Float64Array(heatSource.length),adds=[];
  let totalMassRate=0,alphaSum=0,sigSum=0,count=0,curvatureSum=0,capillarySum=0;
  for(let j=1;j<N-1;j++)for(let i=1;i<N-1;i++)if(isInterfaceGas(state,i,j)){
    const k=id(i,j,N),cGas=cRef*Math.max(1e-9,1+sigma[k]);
    const eq=localInterfaceEquilibrium(state,i,j,iceProps);
    const sigmaDrive=Math.max(0,cGas/eq.cEq-1);
    const alpha=prismAttachmentAt(state,i,j,sigmaDrive);
    const Rdiff=dx/Math.max(1e-30,Dv);
    const Rkin=eq.cEq/Math.max(1e-30,rhoIce*alpha*vkin);
    const massFlux=Math.max(0,cGas-eq.cEq)/(Rdiff+Rkin);
    const interfaceArea=eq.interfaceLength*(2*state.basalHalfThickness);
    const massRate=massFlux*interfaceArea;
    const cellCapacity=rhoIce*dx*dx*(2*state.basalHalfThickness);
    const dphi=massRate*dt/Math.max(1e-30,cellCapacity);
    adds.push([k,dphi,massRate]);
    nextHeat[k]+=massRate*Ls/Math.max(1e-30,dx*dx*(2*state.basalHalfThickness));
    totalMassRate+=massRate;alphaSum+=alpha;sigSum+=sigmaDrive;curvatureSum+=eq.curvature;capillarySum+=eq.capillaryExponent;count++;
  }
  let deposited=0;
  for(const[k,dphi]of adds){
    const before=phi[k];phi[k]=Math.min(1,phi[k]+dphi);if(phi[k]>=0.999)state.ice[k]=1;
    deposited+=(phi[k]-before)*rhoIce*dx*dx*(2*state.basalHalfThickness);
  }
  heatSource.set(nextHeat);

  const meanSigma=count?sigSum/count:0;
  const alphaBasal=attachmentCoefficient(Math.max(0,meanSigma),state.kinetics.basal);
  const vBasal=alphaBasal*vkin*Math.max(0,meanSigma);
  state.basalHalfThickness+=vBasal*dt;
  state.t+=dt;state.totalDepositedMass+=deposited;state.totalLatentHeat+=deposited*Ls;
  let maxDeltaT=0;for(const x of state.dT)if(x>maxDeltaT)maxDeltaT=x;
  state.lastBudget={
    totalMassRate,deposited,Ls,rhoIce,Dv,cRef,vkin,meanSigma,
    meanPrismAlpha:count?alphaSum/count:0,alphaBasal,vBasal,interfaceCells:count,maxDeltaT,
    meanCurvature:count?curvatureSum/count:0,meanCapillaryExponent:count?capillarySum/count:0,
    stencil:'9-point-isotropic',orientation:'phase-gradient',
  };
  return state;
}

export function growthMorphology(state,rays=120){
  const {N,dx,ice,phi}=state,c=(N-1)/2;
  const radius=[];const maxR=c*dx;
  for(let q=0;q<rays;q++){
    const th=2*Math.PI*q/rays;let r=0;
    for(let rr=0;rr<=maxR;rr+=dx*0.35){
      const i=Math.round(c+rr*Math.cos(th)/dx),j=Math.round(c+rr*Math.sin(th)/dx);if(!inb(i,j,N))break;
      const k=id(i,j,N);if(ice[k]||phi[k]>0.15)r=rr;
    }
    radius.push(r);
  }
  const mean=radius.reduce((a,b)=>a+b,0)/radius.length;
  const modes={};
  for(const m of [1,2,3,4,6,8,12]){
    let a=0,b=0;for(let q=0;q<rays;q++){const th=2*Math.PI*q/rays,d=radius[q]-mean;a+=d*Math.cos(m*th);b+=d*Math.sin(m*th);}a*=2/rays;b*=2/rays;
    modes[m]={a,b,amplitude:Math.hypot(a,b),relative:mean?Math.hypot(a,b)/mean:0,phase:Math.atan2(-b,a)/m};
  }
  let areaCells=0;for(let k=0;k<ice.length;k++)areaCells+=ice[k]?1:phi[k];
  const eqRadius=Math.sqrt(areaCells*dx*dx/Math.PI);
  return {radius,meanRadius:mean,modes,m6Amplitude:modes[6].amplitude,m6Phase:modes[6].phase,equivalentBasalRadius:eqRadius,halfThickness:state.basalHalfThickness,aspectRatio:eqRadius/Math.max(1e-30,state.basalHalfThickness)};
}

export function fieldSample(state,i,j){
  const k=id(i,j,state.N),base={uRef:state.sigma[k],deltaT:state.dT[k],ice:!!state.ice[k],phi:state.phi[k],heatSource:state.heatSource[k]};
  if(!state.ice[k]){
    const cRef=saturationMassDensityIce(state.Tinf),cGas=cRef*Math.max(1e-9,1+state.sigma[k]),Tlocal=state.Tinf+state.dT[k],cEq=saturationMassDensityIce(Tlocal);
    base.localSupersaturation=cGas/cEq-1;base.growthPotential=Math.log(Math.max(1e-12,cGas/cEq));
  }
  return base;
}
