import {
  CONSTANTS,
  diffusivityWaterVaporAir,
  thermalConductivityAir,
} from './thermo.js';
import { saturationMassDensityIce } from './postseed.js';

const id=(i,j,N)=>j*N+i;
const inb=(i,j,N)=>i>=0&&j>=0&&i<N&&j<N;
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));

function scalar(state,field,i,j){
  if(!inb(i,j,state.N))return 0;
  const k=id(i,j,state.N);
  if(field==='u')return state.sigma[k];
  if(field==='T')return state.dT[k];
  if(field==='chi'){
    const cRef=saturationMassDensityIce(state.Tinf);
    const cGas=cRef*Math.max(1e-12,1+state.sigma[k]);
    const cEq=saturationMassDensityIce(state.Tinf+state.dT[k]);
    return Math.log(Math.max(1e-14,cGas/cEq));
  }
  return 0;
}

export function scalarGradient(state,field,i,j){
  const h=state.dx;
  const gx=(scalar(state,field,i+1,j)-scalar(state,field,i-1,j))/(2*h);
  const gy=(scalar(state,field,i,j+1)-scalar(state,field,i,j-1))/(2*h);
  return {x:gx,y:gy,magnitude:Math.hypot(gx,gy)};
}

export function vaporMassFluxVector(state,i,j){
  const grad=scalarGradient(state,'u',i,j);
  const D=diffusivityWaterVaporAir(state.Tinf,state.pressurePa);
  const cRef=saturationMassDensityIce(state.Tinf);
  return {x:-D*cRef*grad.x,y:-D*cRef*grad.y,magnitude:D*cRef*grad.magnitude,units:'kg m^-2 s^-1'};
}

export function heatFluxVector(state,i,j){
  const grad=scalarGradient(state,'T',i,j);
  const k=thermalConductivityAir(state.Tinf);
  return {x:-k*grad.x,y:-k*grad.y,magnitude:k*grad.magnitude,units:'W m^-2'};
}

export function growthPotentialVector(state,i,j){
  const g=scalarGradient(state,'chi',i,j);
  return {x:-g.x,y:-g.y,magnitude:g.magnitude,units:'m^-1',meaning:'toward decreasing chemical growth potential'};
}

function bilinearScalar(state,field,x,y){
  const c=(state.N-1)/2;
  const gx=x/state.dx+c,gy=y/state.dx+c;
  const i0=Math.floor(gx),j0=Math.floor(gy),tx=gx-i0,ty=gy-j0;
  if(!inb(i0,j0,state.N)||!inb(i0+1,j0+1,state.N))return null;
  const s00=scalar(state,field,i0,j0),s10=scalar(state,field,i0+1,j0),s01=scalar(state,field,i0,j0+1),s11=scalar(state,field,i0+1,j0+1);
  return s00*(1-tx)*(1-ty)+s10*tx*(1-ty)+s01*(1-tx)*ty+s11*tx*ty;
}

function gradXY(state,field,x,y){
  const h=state.dx*0.5;
  const xp=bilinearScalar(state,field,x+h,y),xm=bilinearScalar(state,field,x-h,y),yp=bilinearScalar(state,field,x,y+h),ym=bilinearScalar(state,field,x,y-h);
  if([xp,xm,yp,ym].some(v=>v==null))return null;
  return {x:(xp-xm)/(2*h),y:(yp-ym)/(2*h)};
}

function occupiedXY(state,x,y){
  const c=(state.N-1)/2,i=Math.round(x/state.dx+c),j=Math.round(y/state.dx+c);
  if(!inb(i,j,state.N))return true;
  const k=id(i,j,state.N);return state.ice[k]||state.phi[k]>0.5;
}

export function traceStreamline(state,{x,y,field='chi',direction=1,step=0.45*state.dx,maxSteps=700,maxLength=null}={}){
  const points=[[x,y]],domain=((state.N-1)/2-1)*state.dx,Lmax=maxLength??domain*5;
  let length=0;
  for(let n=0;n<maxSteps;n++){
    const g=gradXY(state,field,x,y);if(!g)break;
    let vx,vy;
    if(field==='T'){vx=-g.x;vy=-g.y;}
    else{vx=-g.x;vy=-g.y;}
    const m=Math.hypot(vx,vy);if(!(m>1e-18))break;
    vx=direction*vx/m;vy=direction*vy/m;
    // RK2 midpoint
    const mx=x+0.5*step*vx,my=y+0.5*step*vy;
    const gm=gradXY(state,field,mx,my);if(!gm)break;
    let wx=-gm.x,wy=-gm.y,wm=Math.hypot(wx,wy);if(!(wm>1e-18))break;
    wx=direction*wx/wm;wy=direction*wy/wm;
    const xn=x+step*wx,yn=y+step*wy;
    if(Math.abs(xn)>domain||Math.abs(yn)>domain)break;
    points.push([xn,yn]);length+=Math.hypot(xn-x,yn-y);x=xn;y=yn;
    if(occupiedXY(state,x,y)||length>Lmax)break;
  }
  return {field,direction,points,length};
}

export function seedStreamlines(state,{field='chi',count=24,radiusFraction=0.92}={}){
  const domain=((state.N-1)/2-1)*state.dx,r=domain*radiusFraction,lines=[];
  for(let q=0;q<count;q++){
    const th=2*Math.PI*q/count;
    lines.push(traceStreamline(state,{x:r*Math.cos(th),y:r*Math.sin(th),field,direction:1}));
  }
  return lines;
}

function divergenceOfVector(state,vectorFn,i,j){
  const h=state.dx;
  const xp=vectorFn(state,i+1,j),xm=vectorFn(state,i-1,j),yp=vectorFn(state,i,j+1),ym=vectorFn(state,i,j-1);
  return (xp.x-xm.x)/(2*h)+(yp.y-ym.y)/(2*h);
}
function curlOfVector(state,vectorFn,i,j){
  const h=state.dx;
  const xp=vectorFn(state,i+1,j),xm=vectorFn(state,i-1,j),yp=vectorFn(state,i,j+1),ym=vectorFn(state,i,j-1);
  return (xp.y-xm.y)/(2*h)-(yp.x-ym.x)/(2*h);
}

export function fieldResidualDiagnostics(state,margin=2){
  let vaporDiv2=0,heatDiv2=0,chiCurl2=0,count=0,maxV=0,maxQ=0,maxChi=0;
  for(let j=margin;j<state.N-margin;j++)for(let i=margin;i<state.N-margin;i++){
    const k=id(i,j,state.N);if(state.ice[k]||state.phi[k]>0.15)continue;
    const v=vaporMassFluxVector(state,i,j),q=heatFluxVector(state,i,j),c=growthPotentialVector(state,i,j);
    maxV=Math.max(maxV,v.magnitude);maxQ=Math.max(maxQ,q.magnitude);maxChi=Math.max(maxChi,c.magnitude);
    const dv=divergenceOfVector(state,vaporMassFluxVector,i,j);
    const dq=divergenceOfVector(state,heatFluxVector,i,j);
    const cc=curlOfVector(state,growthPotentialVector,i,j);
    vaporDiv2+=dv*dv;heatDiv2+=dq*dq;chiCurl2+=cc*cc;count++;
  }
  return {
    count,
    vaporDivergenceRms:count?Math.sqrt(vaporDiv2/count):0,
    heatDivergenceRms:count?Math.sqrt(heatDiv2/count):0,
    growthPotentialCurlRms:count?Math.sqrt(chiCurl2/count):0,
    maxVaporFlux:maxV,maxHeatFlux:maxQ,maxGrowthPotentialGradient:maxChi,
  };
}

export function interfaceFluxBudget(state){
  let inwardVapor=0,outwardHeat=0,samples=0;
  const c=(state.N-1)/2;
  for(let j=1;j<state.N-1;j++)for(let i=1;i<state.N-1;i++){
    const k=id(i,j,state.N);if(state.ice[k]||state.phi[k]>=0.999)continue;
    let near=false;for(const[di,dj]of [[1,0],[-1,0],[0,1],[0,-1]]){const kk=id(i+di,j+dj,state.N);if(state.ice[kk]||state.phi[kk]>0.15)near=true;}
    if(!near)continue;
    const x=(i-c)*state.dx,y=(j-c)*state.dx,r=Math.hypot(x,y)||1,nx=x/r,ny=y/r;
    const v=vaporMassFluxVector(state,i,j),q=heatFluxVector(state,i,j);
    inwardVapor+=Math.max(0,-(v.x*nx+v.y*ny));
    outwardHeat+=Math.max(0,q.x*nx+q.y*ny);samples++;
  }
  return {samples,meanInwardVaporFlux:samples?inwardVapor/samples:0,meanOutwardHeatFlux:samples?outwardHeat/samples:0};
}
