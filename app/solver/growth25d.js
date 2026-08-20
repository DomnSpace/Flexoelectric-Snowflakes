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
 * Basal plane: Cartesian finite-volume grid with an evolving ice volume fraction.
 *   ∇² sigma = 0                     (quasi-steady vapor field)
 *   ∇² ΔT = -q''' / k_air            (quasi-steady heat field)
 *
 * Interface deposition uses a series resistance:
 *   j = Δsigma / [dx/(D c_sat) + 1/(rho_i alpha v_kin)]
 *
 * and latent source q = j L_s. The grid is 2-D; an effective out-of-plane
 * thickness is carried separately and advanced with basal kinetics. Therefore
 * this is deliberately called 2.5-D, not a 3-D crystal solver.
 */

const dirs=[[1,0],[-1,0],[0,1],[0,-1]];
const id=(i,j,N)=>j*N+i;
const inb=(i,j,N)=>i>=0&&j>=0&&i<N&&j<N;

export const DEFAULT_ICE_KINETICS=Object.freeze({
  prism:{mode:'nucleation',A:1,sigma0:0.003},
  basal:{mode:'nucleation',A:1,sigma0:0.01},
  anisotropy6:0.18,
  sdak:false,
  sdakStrength:0.55,
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

function isBoundary(i,j,N){return i===0||j===0||i===N-1||j===N-1;}
function isInterfaceGas(state,i,j){
  const {N,ice}=state,k=id(i,j,N);if(ice[k])return false;
  return dirs.some(([di,dj])=>inb(i+di,j+dj,N)&&ice[id(i+di,j+dj,N)]);
}
function interfaceOrientation(state,i,j){
  const {N,dx}=state,c=(N-1)/2,x=(i-c)*dx,y=(j-c)*dx;
  return Math.atan2(y,x);
}
function localPrismAlpha(state,i,j,sigmaDrive){
  const th=interfaceOrientation(state,i,j),k=state.kinetics;
  let model={...k.prism};
  if(k.sdak){
    let iceN=0;for(const[di,dj]of dirs)if(inb(i+di,j+dj,state.N)&&state.ice[id(i+di,j+dj,state.N)])iceN++;
    const edgeFactor=Math.max(0,1-iceN/4);
    model.sigma0*=Math.max(0.08,1-(k.sdakStrength??0.55)*edgeFactor);
  }
  const base=attachmentCoefficient(Math.max(0,sigmaDrive),model);
  const six=1+(k.anisotropy6??0)*Math.cos(6*th);
  return Math.max(0,Math.min(1,base*six));
}

export function relaxVaporField(state,iterations=80){
  const {N,ice,sigma,sigmaInf,Tinf}=state;
  const a=new Float64Array(sigma),b=new Float64Array(sigma.length);
  const Ls=sublimationEnthalpy(Tinf,state.pressurePa).Ls,beta=Ls/(CONSTANTS.RV*Tinf*Tinf);
  for(let it=0;it<iterations;it++){
    for(let j=0;j<N;j++)for(let i=0;i<N;i++){
      const k=id(i,j,N);
      if(isBoundary(i,j,N)){b[k]=sigmaInf;continue;}
      if(ice[k]){b[k]=beta*state.dT[k];continue;}
      let s=0,n=0;for(const[di,dj]of dirs){s+=a[id(i+di,j+dj,N)];n++;}b[k]=s/n;
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
      let s=0;for(const[di,dj]of dirs)s+=a[id(i+di,j+dj,N)];
      b[k]=0.25*(s+heatSource[k]*dx*dx/kAir);
    }
    a.set(b);
  }
  dT.set(a);return state;
}

export function resolvedGrowthStep(state,dt,{vaporIterations=45,heatIterations=55}={}){
  const {N,dx,Tinf,pressurePa,ice,phi,sigma,heatSource}=state;
  relaxHeatField(state,heatIterations);relaxVaporField(state,vaporIterations);
  const iceProps=iceIhProperties(Tinf,Math.min(pressurePa,210e6)),rhoIce=iceProps.rho;
  const Dv=diffusivityWaterVaporAir(Tinf,pressurePa),cSat=saturationMassDensityIce(Tinf),vkin=iceKineticVelocity(Tinf,rhoIce),Ls=sublimationEnthalpy(Tinf,pressurePa).Ls;
  const beta=Ls/(CONSTANTS.RV*Tinf*Tinf);
  const nextHeat=new Float64Array(heatSource.length),adds=[];
  let totalMassRate=0,alphaSum=0,sigSum=0,count=0;
  for(let j=1;j<N-1;j++)for(let i=1;i<N-1;i++)if(isInterfaceGas(state,i,j)){
    const k=id(i,j,N),localSigma=sigma[k];
    let faces=0,sigmaEqSum=0;
    for(const[di,dj]of dirs){const kk=id(i+di,j+dj,N);if(ice[kk]){faces++;sigmaEqSum+=beta*state.dT[kk];}}
    const sigmaEq=faces?sigmaEqSum/faces:0;
    const sigmaDrive=Math.max(0,localSigma-sigmaEq);
    const alpha=localPrismAlpha(state,i,j,sigmaDrive);
    const Rdiff=dx/(Math.max(1e-30,Dv*cSat)),Rkin=1/(Math.max(1e-30,rhoIce*alpha*vkin));
    const massFlux=sigmaDrive/(Rdiff+Rkin); // kg m^-2 s^-1
    const interfaceArea=faces*dx*(2*state.basalHalfThickness);
    const massRate=massFlux*interfaceArea;
    const cellCapacity=rhoIce*dx*dx*(2*state.basalHalfThickness);
    const dphi=massRate*dt/Math.max(1e-30,cellCapacity);
    adds.push([k,dphi,massRate]);
    nextHeat[k]+=massRate*Ls/Math.max(1e-30,dx*dx*(2*state.basalHalfThickness));
    totalMassRate+=massRate;alphaSum+=alpha;sigSum+=sigmaDrive;count++;
  }
  let deposited=0;
  for(const[k,dphi]of adds){const before=phi[k];phi[k]=Math.min(1,phi[k]+dphi);if(phi[k]>=0.999)ice[k]=1;deposited+=(phi[k]-before)*rhoIce*dx*dx*(2*state.basalHalfThickness);}
  heatSource.set(nextHeat);

  const meanSigma=count?sigSum/count:0;
  const alphaBasal=attachmentCoefficient(Math.max(0,meanSigma),state.kinetics.basal);
  const vBasal=alphaBasal*vkin*Math.max(0,meanSigma);
  state.basalHalfThickness+=vBasal*dt;
  state.t+=dt;state.totalDepositedMass+=deposited;state.totalLatentHeat+=deposited*Ls;
  state.lastBudget={totalMassRate,deposited,Ls,rhoIce,Dv,cSat,vkin,meanSigma,meanPrismAlpha:count?alphaSum/count:0,alphaBasal,vBasal,interfaceCells:count,maxDeltaT:Math.max(...state.dT)};
  return state;
}

export function growthMorphology(state,rays=120){
  const {N,dx,ice,phi}=state,c=(N-1)/2;
  const radius=[];
  const maxR=c*dx;
  for(let q=0;q<rays;q++){
    const th=2*Math.PI*q/rays;let r=0;
    for(let rr=0;rr<=maxR;rr+=dx*0.35){const i=Math.round(c+rr*Math.cos(th)/dx),j=Math.round(c+rr*Math.sin(th)/dx);if(!inb(i,j,N))break;const k=id(i,j,N);if(ice[k]||phi[k]>0.15)r=rr;}
    radius.push(r);
  }
  const mean=radius.reduce((a,b)=>a+b,0)/radius.length;
  let a6=0,b6=0;for(let q=0;q<rays;q++){const th=2*Math.PI*q/rays,d=radius[q]-mean;a6+=d*Math.cos(6*th);b6+=d*Math.sin(6*th);}a6*=2/rays;b6*=2/rays;
  const eqRadius=Math.sqrt(state.ice.reduce((s,v,k)=>s+(v?1:state.phi[k]),0)*dx*dx/Math.PI);
  return {radius,meanRadius:mean,m6Amplitude:Math.hypot(a6,b6),m6Phase:Math.atan2(-b6,a6)/6,equivalentBasalRadius:eqRadius,halfThickness:state.basalHalfThickness,aspectRatio:eqRadius/Math.max(1e-30,state.basalHalfThickness)};
}

export function fieldSample(state,i,j){const k=id(i,j,state.N);return{sigma:state.sigma[k],deltaT:state.dT[k],ice:!!state.ice[k],phi:state.phi[k],heatSource:state.heatSource[k]};}
