/**
 * Model B / Step 3: backward-Euler, finite-volume vapor transport through
 * actual triangular x z gas control volumes. 5 independent ice-eligible slabs
 * compete for ONE vapor inventory. This is NOT the legacy hex25d solver.
 *
 * Implicit solve in u=c-c_eq >=0:
 * (I + dt/V * L_gas + dt/V * K_attach + dt/V * K_res) u_new
 *       = u_old + (dt/V) K_res u_res.
 * A is SPD with no artificial constant-mode pin. Jacobi-preconditioned CG
 * checks the algebraic residual before ANY physical state is committed.
 *
 * Model restrictions: uniform dz, uniform T, deposition only (u>=0),
 * constant positive diffusivity/attachment velocity, fixed geometry within
 * each implicit step; no heat PDE, no micro-terrace kinetics, no full 3-D ice.
 */
import {createLayeredGasFV,layeredGasMassKg,layeredTrackedWaterKg}
  from './layered-gas-fv.js';
import {iceIhProperties} from './ice-ih-eos.js';

const SQRT3=Math.sqrt(3);
const DIR=[[1,0],[-1,0],[0,1],[0,-1],[1,-1],[-1,1]];
const facePairs=[[1,0],[0,1],[1,-1]];
const nzSum=a=>a.reduce((s,v)=>s+v,0);
const validate=(value,label)=>{if(!(Number.isFinite(value)&&value>0))throw new RangeError(label+' must be positive');};
function indexed(s,i,j,z){return z*s.nxy+j*s.N+i;}
function neighbor(s,i,j,z,dq,dr){
  const q=i+dq,r=j+dr;if(q<0||q>=s.N||r<0||r>=s.N||!s.hexMask[r*s.N+q])return -1;
  return indexed(s,q,r,z);
}
function active(s,k){return k>=0&&!s.solid[k];}
function initialLayout(s){
  const nodes=[],map=new Int32Array(s.c.length).fill(-1);
  for(let z=0;z<s.nz;z++)for(let j=0;j<s.N;j++)for(let i=0;i<s.N;i++){
    const p=j*s.N+i;
    if(!s.hexMask[p])continue;
    const k=indexed(s,i,j,z);
    if(active(s,k)){map[k]=nodes.length;nodes.push(k);}
  }
  return {nodes,map};
}
export function assembleLayeredImplicit(s,dt){
  validate(dt,'dt');
  if(s.schema!=='snowflake/vertical-gas-fv-step2')throw new TypeError('Expected Step-2 gas FV state');
  const {nodes,map}=initialLayout(s);
  const n=nodes.length,diag=new Float64Array(n).fill(1),local=new Float64Array(n).fill(1),rhs=new Float64Array(n);
  const edges=[],sinks=[],reservoirs=[];
  const coef=dt/s.volume;
  const lateral=s.D*s.sideArea/s.dx*coef,vertical=s.D*s.area/s.dz*coef;
  for(let m=0;m<n;m++){
    const k=nodes[m],z=Math.floor(k/s.nxy),p=k-z*s.nxy;
    const i=p%s.N,j=Math.floor(p/s.N);
    rhs[m]=s.c[k]-s.ceq;
    if(rhs[m]<-1e-14)throw new RangeError('Condensation-only implicit solver requires c >= c_eq');
    rhs[m]=Math.max(0,rhs[m]);
    function coupled(other,weight){
      if(other<0||!s.solid[other]){
        if(other<0)return;
        const nb=map[other];
        if(nb<0)throw new Error('Invalid active neighbor mapping');
        if(nb>m){
          edges.push([m,nb,weight]);
          diag[m]+=weight;diag[nb]+=weight;
        }
      }else{
        const basal=Math.abs(other-k)===s.nxy;
        const attachment=(basal?s.kBasal*s.area:s.kPrism*s.sideArea)*coef;
        if(attachment>0){
          diag[m]+=attachment;local[m]+=attachment;
          sinks.push([m,other,attachment]);
        }
      }
    }
    function border(which,faceArea,halfWidth){
      if(s.boundaries[which]!=='reservoir')return;
      const weight=s.D*faceArea/halfWidth*coef;
      const cExternal=s.boundaries[which+'ConcentrationKgM3'];
      if(cExternal<s.ceq-1e-14)
        throw new RangeError('Condensation-only reservoir must be >= equilibrium concentration');
      diag[m]+=weight;local[m]+=weight;
      rhs[m]+=weight*Math.max(0,cExternal-s.ceq);
      reservoirs.push([m,which,weight,cExternal]);
    }
    for(const [dq,dr] of DIR){
      const other=neighbor(s,i,j,z,dq,dr);
      if(other<0)border('lateral',s.sideArea,s.dx/2);
      else coupled(other,lateral);
    }
    for(const dz of [-1,1]){
      const zz=z+dz;
      if(zz<0||zz>=s.nz)border(zz<0?'bottom':'top',s.area,s.dz/2);
      else coupled(indexed(s,i,j,zz),vertical);
    }
  }
  return {nodes,map,diag,local,rhs,edges,sinks,reservoirs,coef,dt,
    gasCells:n,meshFullCells:s.c.length};
}
function matvec(A,x,y){
  // Evaluate the diffusion operator as differences of neighboring values:
  //  y_i = (1+sink+reservoir)*x_i + Σ_ij w_ij*(x_i-x_j).
  // Computing diag*x_i - Σ w*x_j loses 8+ significant digits when
  // dt/dtExplicit is large and x varies very slowly. This flux form
  // preserves the exact constant nullspace of diffusion up to roundoff.
  for(let i=0;i<x.length;i++)y[i]=A.local[i]*x[i];
  for(const [i,j,w] of A.edges){
    const flux=w*(x[i]-x[j]);
    y[i]+=flux;y[j]-=flux;
  }
}
function dot(a,b){let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s;}
function norm(a){return Math.sqrt(dot(a,a));}
export function solveLayeredImplicit(s,dt,{
  tolerance=1e-8,maxIterations=1600
}={}){
  validate(tolerance,'tolerance');
  if(!Number.isInteger(maxIterations)||maxIterations<1)
    throw new RangeError('Invalid maximum CG iteration count');
  const A=assembleLayeredImplicit(s,dt);
  const n=A.nodes.length,x=new Float64Array(n),ap=new Float64Array(n);
  const residual=new Float64Array(n),z=new Float64Array(n),p=new Float64Array(n);
  for(let i=0;i<n;i++)x[i]=Math.max(0,s.c[A.nodes[i]]-s.ceq);
  matvec(A,x,ap);
  let rr0=0,bb=0;
  for(let i=0;i<n;i++){
    residual[i]=A.rhs[i]-ap[i];z[i]=residual[i]/A.diag[i];p[i]=z[i];
    rr0+=residual[i]*z[i];bb+=A.rhs[i]*A.rhs[i];
  }
  const bNorm=Math.sqrt(bb),target=tolerance*Math.max(bNorm,1e-30);
  let error=norm(residual),iterations=0,rho=rr0;
  while(error>target&&iterations<maxIterations){
    matvec(A,p,ap);
    const denom=dot(p,ap);
    if(!(denom>0&&Number.isFinite(denom)))throw new Error('Non-SPD or non-finite vapor operator');
    const alpha=rho/denom;
    for(let i=0;i<n;i++){x[i]+=alpha*p[i];residual[i]-=alpha*ap[i];}
    error=norm(residual);iterations++;
    if(error<=target)break;
    let next=0;
    for(let i=0;i<n;i++){z[i]=residual[i]/A.diag[i];next+=residual[i]*z[i];}
    const beta=next/rho;
    for(let i=0;i<n;i++)p[i]=z[i]+beta*p[i];
    rho=next;
  }
  // Project the common-mode numerical drift onto the exact finite-volume
  // water budget. Pairwise gas-gas terms sum to zero, so the mass-equation
  // sum is Σ_i local_i*x_i = Σ_i rhs_i. A constant shift changes that
  // sum by shift*Σ_i local_i without changing any diffusion gradient.
  // This is a conservation projection, NOT a substitute for checking the
  // full algebraic residual below.
  let rhsSum=0,lhsSum=0,localSum=0;
  for(let i=0;i<n;i++){
    rhsSum+=A.rhs[i];lhsSum+=A.local[i]*x[i];localSum+=A.local[i];
  }
  const correction=(rhsSum-lhsSum)/localSum;
  if(!Number.isFinite(correction))throw new Error('Non-finite vapor mass projection');
  for(let i=0;i<n;i++)x[i]+=correction;
  // Compute TRUE final residual (not only recursively propagated CG residual).
  matvec(A,x,ap);
  let true2=0;
  for(let i=0;i<n;i++){
    const d=A.rhs[i]-ap[i];true2+=d*d;
    if(!Number.isFinite(x[i])||x[i]<-1e-10)
      throw new Error('Nonpositive vapor field in implicit solver');
  }
  const trueResidual=Math.sqrt(true2)/Math.max(bNorm,1e-30);
  if(trueResidual>tolerance*1.25)
    throw new Error('Implicit vapor did not converge within '+maxIterations+
      ' CG iterations: relative residual '+trueResidual);
  return {A,x,iterations,relativeResidual:trueResidual,conservationShiftKgM3:correction};
}

/**
 * Atomically commits one converged implicit solve: no partial mutations if
 * CG fails or a temperature/diffusivity boundary is outside its assumptions.
 * Deposition goes to the FIXED solid ice site; the moving-front adapter uses
 * the separate resolved gas/ice face events to advance five independent
 * projected occupancy fields later (without inventing vapor).
 */
export function stepLayeredImplicit(s,dt,opts={}){
  const before=layeredTrackedWaterKg(s);
  const {A,x,iterations,relativeResidual,conservationShiftKgM3}=solveLayeredImplicit(s,dt,opts);
  const cvMass=s.volume,events=[],perIceSlab=new Float64Array(s.iceSlabs);
  let massToIce=0,boundaryInKg=0,boundaryOutKg=0;
  const nextC=new Float64Array(x.length);
  for(let i=0;i<x.length;i++)nextC[i]=Math.max(0,s.ceq+x[i]);
  for(const [gas,ice,w] of A.sinks){
    const deposited=dt/A.coef*w*Math.max(0,x[gas]);
    if(deposited<=0)continue;
    const z=Math.floor(ice/s.nxy),iceSlab=z-s.gasHaloSlabs;
    if(iceSlab<0||iceSlab>=s.iceSlabs)throw new Error('Deposit target outside eligible ice slabs');
    events.push({fromGasIndex:A.nodes[gas],toIceIndex:ice,iceLayer:iceSlab,
      massKg:deposited});
    perIceSlab[iceSlab]+=deposited;massToIce+=deposited;
  }
  for(const [gas,side,w,externalConcentration] of A.reservoirs){
    const mass=dt/A.coef*w*(externalConcentration-nextC[gas]);
    if(mass>0)boundaryInKg+=mass;
    else boundaryOutKg-=mass;
  }
  // Validate budget BEFORE changing any physical data.
  const newGasMass=x.reduce((total,u,i)=>total+(s.ceq+u)*cvMass,0);
  const expected=newGasMass+massToIce+boundaryOutKg-boundaryInKg;
  const residualKg=expected-layeredGasMassKg(s);
  const budgetScale=Math.max(before,1e-30);
  if(Math.abs(residualKg)>budgetScale*Math.max(opts.tolerance??1e-10,1e-10)*6)
    throw new Error('Implicit finite-volume mass residual too large: '+residualKg);
  for(let i=0;i<x.length;i++)s.c[A.nodes[i]]=nextC[i];
  for(const e of events){
    s.iceDepositedKg[e.toIceIndex]+=e.massKg;
    s.iceDepositionPerSlabKg[e.iceLayer+s.gasHaloSlabs]+=e.massKg;
    s.latentReleasedPerSlabJ[e.iceLayer+s.gasHaloSlabs]+=e.massKg*s.Ls;
  }
  s.totalDepositedKg+=massToIce;s.latentReleasedJ+=massToIce*s.Ls;
  s.externalMassInKg+=boundaryInKg;s.externalMassOutKg+=boundaryOutKg;
  s.timeSeconds+=dt;
  const report={dtSeconds:dt,solver:'backward-euler-pcg',iterations,
    relativeResidual,conservationShiftKgM3,massResidualKg:residualKg,faces:events.length,
    depositedKg:massToIce,perIceSlabDepositedKg:Array.from(perIceSlab),
    boundaryMassInKg:boundaryInKg,boundaryMassOutKg:boundaryOutKg,
    gasMassKg:layeredGasMassKg(s),latentReleasedJ:massToIce*s.Ls,
    heatTransportSolved:false,interfaceMotionSolved:false,
    depositionEvents:events};
  s.lastStep=report;
  return report;
}
export {createLayeredGasFV,layeredGasMassKg,layeredTrackedWaterKg};
