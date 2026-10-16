/**
 * Model B / Step 4 thermal pilot: coarse gas/ice voxel enthalpy conduction.
 *
 * A separate implicit finite-volume heat equation is solved on the SAME
 * 2D triangular x vertical slab geometry as the independent ice fields.
 * Deposition releases L_s*dm at the accepting ice voxel. Thermal energy
 * (J relative to Tinf) is the primary state: changing ice fraction does not
 * create/destroy heat merely by changing cell heat capacity.
 *
 * This is a homogenized mesoscale closure:
 * rhoAir*cpAir for gas, IAPWS Ice Ih cp and rho for ice,
 * volume-weighted conductivity (then harmonic face mean).
 * It does NOT resolve contact resistances, atmospheric convection or a
 * molecular ice surface. It uses fixed Tinf reference properties.
 */
import {thermalConductivityAir} from './thermo.js';
import {iceIhProperties} from './ice-ih-eos.js';
import {saturationMassDensityIce} from './postseed.js';

const SQRT3=Math.sqrt(3),LATERAL=[[1,0],[0,1],[1,-1]];
const safe=(v,n)=>{if(!Number.isFinite(v)||v<=0)throw new RangeError(n+' must be positive finite');};
const sum=a=>{let v=0;for(const x of a)v+=x;return v;};

export function attachLayeredThermal(s,{iceConductivityWmK=2.2,cpDryAirJkgK=1005,
  airGasConstantJkgK=287.05}={}){
  safe(iceConductivityWmK,'ice conductivity');
  safe(cpDryAirJkgK,'air heat capacity');
  safe(airGasConstantJkgK,'air specific gas constant');
  if(!s.iceFraction||!s.iceMassAtCellKg)throw new TypeError('Heat model needs the 5-slab volume fractions');
  const ice=iceIhProperties(s.Tinf,Math.min(s.pressurePa,210e6));
  s.thermal={
    kIce:iceConductivityWmK,kAir:thermalConductivityAir(s.Tinf),
    cpIce:ice.cp,cpAir:cpDryAirJkgK,
    rhoAir:s.pressurePa/(airGasConstantJkgK*s.Tinf)
  };
  s.temperatureDeltaK=new Float64Array(s.c.length);
  s.thermalEnthalpyJ=new Float64Array(s.c.length);
  s.localIceEqKgM3=new Float64Array(s.c.length);
  s.localIceEqKgM3.fill(s.ceq);
  s.thermalBoundaryHeatOutJ=0;
  s.thermalLatentInjectedJ=0;
  s.heatTransportSolved=true;
  return s;
}
export function layeredThermalEnergyJ(s){return sum(s.thermalEnthalpyJ||[]);}
export function layeredThermalCapacityAndK(s){
  if(!s.thermal)throw new TypeError('Thermal properties not initialized');
  const a=new Float64Array(s.c.length),kappa=new Float64Array(s.c.length);
  const t=s.thermal;
  for(let z=0;z<s.nz;z++)for(let p=0;p<s.nxy;p++){
    if(!s.hexMask[p])continue;
    const i=z*s.nxy+p,f=s.iceFraction[i],gas=1-f;
    a[i]=s.volume*(gas*t.rhoAir*t.cpAir+f*s.rhoIce*t.cpIce);
    kappa[i]=gas*t.kAir+f*t.kIce;
    safe(a[i],'heat capacity');
  }
  return {capacity:a,conductivity:kappa};
}
const harmonic=(a,b)=>2*a*b/(a+b);
function buildHeatMatrix(s,dt,source){
  safe(dt,'dt');
  const {capacity,conductivity}=layeredThermalCapacityAndK(s);
  const nodes=[],map=new Int32Array(s.c.length).fill(-1);
  for(let z=0;z<s.nz;z++)for(let p=0;p<s.nxy;p++)
    if(s.hexMask[p]){map[z*s.nxy+p]=nodes.length;nodes.push(z*s.nxy+p);}
  const n=nodes.length,local=new Float64Array(n),diag=new Float64Array(n);
  const rhs=new Float64Array(n),edges=[],reservoirs=[];
  for(let m=0;m<n;m++){
    const idx=nodes[m],z=Math.floor(idx/s.nxy),p=idx%s.nxy,
      i=p%s.N,j=Math.floor(p/s.N);
    local[m]=diag[m]=capacity[idx];
    rhs[m]=s.thermalEnthalpyJ[idx]+(source?.[idx]||0);
    function pair(other,area,distance){
      if(other<0)return;
      const q=map[other];if(q<=m)return;
      const w=dt*harmonic(conductivity[idx],conductivity[other])*area/distance;
      edges.push([m,q,w]);diag[m]+=w;diag[q]+=w;
    }
    function outer(face,area,dist){
      if(s.boundaries[face]!=='reservoir')return;
      const w=dt*conductivity[idx]*area/dist;
      diag[m]+=w;local[m]+=w;
      reservoirs.push([m,face,w]);
    }
    for(const [dq,dr] of LATERAL){
      const ii=i+dq,jj=j+dr;
      if(ii<0||jj<0||ii>=s.N||jj>=s.N||!s.hexMask[jj*s.N+ii])
        outer('lateral',s.sideArea,s.dx/2);
      else pair(z*s.nxy+jj*s.N+ii,s.sideArea,s.dx);
    }
    // For exterior faces in the three remaining lateral directions we
    // add boundary heat transfers separately (interior faces exactly once).
    for(const [dq,dr] of [[-1,0],[0,-1],[-1,1]]){
      const ii=i+dq,jj=j+dr;
      if(ii<0||jj<0||ii>=s.N||jj>=s.N||!s.hexMask[jj*s.N+ii])
        outer('lateral',s.sideArea,s.dx/2);
    }
    if(z>0)pair(idx-s.nxy,s.area,s.dz);
    else outer('bottom',s.area,s.dz/2);
    if(z+1<s.nz)pair(idx+s.nxy,s.area,s.dz);
    else outer('top',s.area,s.dz/2);
  }
  return {nodes,map,local,diag,rhs,edges,reservoirs,capacity,dt};
}
function dot(a,b){let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s;}
function apply(A,x,y){
  for(let i=0;i<x.length;i++)y[i]=A.local[i]*x[i];
  for(const [i,j,w] of A.edges){const flux=w*(x[i]-x[j]);y[i]+=flux;y[j]-=flux;}
}
export function stepLayeredThermal(s,dt,events,{
  tolerance=1e-9,maxIterations=2400
}={}){
  if(!s.thermal)throw new TypeError('Thermal model not attached');
  if(!Array.isArray(events))throw new TypeError('Supply measured vapor→ice deposition events');
  safe(dt,'dt');safe(tolerance,'tolerance');
  const source=new Float64Array(s.c.length);
  let latentJ=0;
  for(const e of events){
    if(!(e.massKg>=0&&Number.isFinite(e.massKg)))throw new RangeError('Invalid latent deposition mass');
    const energy=e.massKg*s.Ls;
    source[e.toIceIndex]+=energy;latentJ+=energy;
  }
  const A=buildHeatMatrix(s,dt,source),n=A.nodes.length;
  const x=new Float64Array(n),ax=new Float64Array(n),r=new Float64Array(n),z=new Float64Array(n),p=new Float64Array(n);
  let normB=0,rho=0;
  for(let i=0;i<n;i++){
    const k=A.nodes[i];
    x[i]=s.thermalEnthalpyJ[k]/A.capacity[k];
  }
  apply(A,x,ax);
  for(let i=0;i<n;i++){
    r[i]=A.rhs[i]-ax[i];z[i]=r[i]/A.diag[i];p[i]=z[i];
    normB+=A.rhs[i]**2;rho+=r[i]*z[i];
  }
  normB=Math.sqrt(normB);
  const target=tolerance*Math.max(normB,1e-30);
  let error=Math.sqrt(dot(r,r)),iterations=0;
  while(error>target&&iterations<maxIterations){
    apply(A,p,ax);
    const denom=dot(p,ax);
    if(!(denom>0&&Number.isFinite(denom)))throw new Error('Heat operator is not SPD');
    const alpha=rho/denom;
    for(let i=0;i<n;i++){x[i]+=alpha*p[i];r[i]-=alpha*ax[i];}
    error=Math.sqrt(dot(r,r));iterations++;
    if(error<=target)break;
    let next=0;
    for(let i=0;i<n;i++){z[i]=r[i]/A.diag[i];next+=r[i]*z[i];}
    const beta=next/rho;
    for(let i=0;i<n;i++)p[i]=z[i]+beta*p[i];
    rho=next;
  }
  // Separate common-mode heat-conservation projection and TRUE residual.
  let rsum=0,lsum=0,capacity=0;
  for(let i=0;i<n;i++){
    rsum+=A.rhs[i];lsum+=A.local[i]*x[i];capacity+=A.local[i];
  }
  const correction=(rsum-lsum)/capacity;
  for(let i=0;i<n;i++)x[i]+=correction;
  apply(A,x,ax);
  let true2=0,maxDeltaK=0,energy=0;
  for(let i=0;i<n;i++){
    const d=A.rhs[i]-ax[i];true2+=d*d;
    if(!Number.isFinite(x[i])||s.Tinf+x[i]<110||s.Tinf+x[i]>273.15)
      throw new Error('Thermal step outside ice equilibrium temperature range');
    energy+=x[i]*A.capacity[A.nodes[i]];
    maxDeltaK=Math.max(maxDeltaK,x[i]);
  }
  const relativeResidual=Math.sqrt(true2)/Math.max(normB,1e-30);
  if(relativeResidual>tolerance*1.3)
    throw new Error('Thermal implicit step did not converge: '+relativeResidual);
  let exported=0;
  for(const [m,face,w] of A.reservoirs)exported+=w*x[m];
  const oldEnergy=layeredThermalEnergyJ(s);
  const residual=energy+exported-oldEnergy-latentJ;
  const scale=Math.max(Math.abs(energy)+Math.abs(exported)+Math.abs(oldEnergy)+Math.abs(latentJ),1e-30);
  if(Math.abs(residual)>scale*Math.max(1e-7,5*tolerance))
    throw new Error('Heat-energy balance failure: '+residual);
  // Atomic commit only after finite-temperature, linear-residual and energy tests.
  for(let i=0;i<n;i++){
    const cell=A.nodes[i];
    s.temperatureDeltaK[cell]=x[i];
    s.thermalEnthalpyJ[cell]=x[i]*A.capacity[cell];
  }
  s.thermalBoundaryHeatOutJ+=exported;
  s.thermalLatentInjectedJ+=latentJ;
  // The local equilibrium condition is evaluated at the ICE temperature
  // on the absorbing facet, never simply at the global gas reservoir.
  for(let k=0;k<s.localIceEqKgM3.length;k++){
    if(!s.hexMask[k%s.nxy])continue;
    if(s.iceFraction[k]>0)
      s.localIceEqKgM3[k]=saturationMassDensityIce(s.Tinf+s.temperatureDeltaK[k]);
  }
  return {solver:'implicit-thermal-voxel-pcg',iterations,relativeResidual,
    depositedLatentJ:latentJ,heatBoundaryOutJ:exported,
    energyBalanceResidualJ:residual,energyStoredJ:energy,maxDeltaK,
    interfaceEquilibriumUpdated:true};
}
