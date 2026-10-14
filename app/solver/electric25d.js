/**
 * 2.5-D electrostatic extension for the existing resolved ice-growth solver.
 *
 * Solve ∇·D = 0 with D = ε0 εr E + P, E = -∇V:
 *   ∇·(εr ∇V) = ∇·P / ε0
 *
 * εr(φ) = ε_air + (ε_ice-ε_air)*occupancy. Harmonic face
 * permittivities enforce continuity of normal D across cells.
 *
 * The optional *prescribed* in-plane strain gradient G has units m^-1,
 * P = μ G [C/m²] within ice and vanishes in air. No elasticity equation
 * is solved; this is a controlled bending analogue, not a prediction
 * of strain gradients from snowflake surface curvature.
 *
 * Optional dielectric free energy is a clearly labeled approximate
 * fixed-field equilibrium closure, NOT a calibrated growth law.
 */
export const VACUUM_PERMITTIVITY = 8.8541878128e-12; // F/m
export const BOLTZMANN = 1.380649e-23;              // J/K
export const WATER_MOLECULE_MASS = 0.01801528/6.02214076e23; // kg
const PI = Math.PI;
const clamp=(v,lo,hi)=>Math.max(lo,Math.min(hi,v));
const id=(i,j,N)=>j*N+i;
const face=(a,b)=>2*a*b/(a+b);
const defaultExperiment={
  enabled:false,
  fieldVM:0,angleDeg:0,
  epsilonIce:3.2, epsilonAir:1.0006,
  prescribedGradient:0,gradientAngleDeg:0,
  flexoCoefficient:1.65e-9,
  coupleGrowth:false
};
export const DEFAULT_ELECTRIC_EXPERIMENT=Object.freeze({...defaultExperiment});

function occupancy(s,k){return s.ice[k]?1:clamp(s.phi[k],0,1);}
function validate(c){
  for(const name of ['fieldVM','angleDeg','epsilonIce','epsilonAir','prescribedGradient','gradientAngleDeg','flexoCoefficient']){
    if(!Number.isFinite(c[name]))throw new RangeError('Electric parameter '+name+' must be finite');
  }
  if(c.fieldVM<0||c.fieldVM>2e7)throw new RangeError('Electric magnitude must be 0..20 MV/m');
  if(c.epsilonIce<1||c.epsilonIce>120||c.epsilonAir<1||c.epsilonAir>2)throw new RangeError('Invalid permittivity');
  if(c.prescribedGradient<0||c.prescribedGradient>1e5)throw new RangeError('Invalid prescribed strain gradient');
  if(c.flexoCoefficient<0||c.flexoCoefficient>1e-7)throw new RangeError('Invalid flexoelectric coefficient');
}
export function initializeElectrostatics(state,params={}){
  const c={...defaultExperiment,...params};
  validate(c);
  const N=state.N,n=N*N,dx=state.dx,mid=(N-1)/2;
  const potential=new Float64Array(n),Ex=new Float64Array(n),Ey=new Float64Array(n);
  const magnitude=new Float64Array(n),epsilon=new Float64Array(n);
  const theta=c.angleDeg*PI/180,ExInf=c.fieldVM*Math.cos(theta),EyInf=c.fieldVM*Math.sin(theta);
  const gradientTheta=c.gradientAngleDeg*PI/180;
  const Px0=c.flexoCoefficient*c.prescribedGradient*Math.cos(gradientTheta);
  const Py0=c.flexoCoefficient*c.prescribedGradient*Math.sin(gradientTheta);
  const Px=new Float64Array(n),Py=new Float64Array(n);
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    const k=id(i,j,N),q=occupancy(state,k),x=(i-mid)*dx,y=(j-mid)*dx;
    potential[k]=-ExInf*x-EyInf*y;
    epsilon[k]=c.epsilonAir+(c.epsilonIce-c.epsilonAir)*q;
    Px[k]=Px0*q;Py[k]=Py0*q;
  }
  state.electric={params:c,potential,Ex,Ey,magnitude,epsilon,Px,Py,
    diagnostics:{residualRmsV:0,maxFieldVM:0,meanFieldVM:0,maxPolarization:Math.hypot(Px0,Py0),sweeps:0}};
  if(c.enabled)solveElectrostatics(state,160);
  return state.electric;
}

/** Recompute permittivity and imposed polarization after the ice mask grows. */
function updateMaterial(state){
  const e=state.electric,c=e.params,n=state.N*state.N;
  const a=c.gradientAngleDeg*PI/180;
  const px=c.flexoCoefficient*c.prescribedGradient*Math.cos(a),py=c.flexoCoefficient*c.prescribedGradient*Math.sin(a);
  for(let k=0;k<n;k++){
    const q=occupancy(state,k);
    e.epsilon[k]=c.epsilonAir+(c.epsilonIce-c.epsilonAir)*q;
    e.Px[k]=px*q;e.Py[k]=py*q;
  }
}
export function solveElectrostatics(state,iterations=80){
  const e=state.electric;
  if(!e||!e.params.enabled)return e;
  updateMaterial(state);
  const {potential:V,epsilon:eps,Px,Py}=e;
  const N=state.N,dx=state.dx,c=e.params;
  const th=c.angleDeg*PI/180,extX=c.fieldVM*Math.cos(th),extY=c.fieldVM*Math.sin(th);
  const mid=(N-1)/2,omega=1.55;
  // The solution has Dirichlet far-field voltage (linear potential).
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    if(i===0||j===0||i===N-1||j===N-1)
      V[id(i,j,N)]=-extX*(i-mid)*dx-extY*(j-mid)*dx;
  }
  function point(i,j){
    const k=id(i,j,N),eId=k+1,wId=k-1,nId=k-N,sId=k+N;
    const eFace=face(eps[k],eps[eId]),wFace=face(eps[k],eps[wId]);
    const nFace=face(eps[k],eps[nId]),sFace=face(eps[k],eps[sId]);
    const sum=eFace+wFace+nFace+sFace;
    const rhs=dx*((Px[eId]-Px[wId])+(Py[sId]-Py[nId]))/(2*VACUUM_PERMITTIVITY);
    return (eFace*V[eId]+wFace*V[wId]+nFace*V[nId]+sFace*V[sId]-rhs)/sum;
  }
  const sweeps=Math.max(0,Math.floor(iterations));
  for(let n=0;n<sweeps;n++){
    for(let j=1;j<N-1;j++)for(let i=1;i<N-1;i++){
      const k=id(i,j,N),next=point(i,j);
      V[k]+=omega*(next-V[k]);
    }
  }
  let residualSq=0,cellCount=0,maxField=0,sumField=0;
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    const k=id(i,j,N);
    if(i>0&&j>0&&i<N-1&&j<N-1){const err=point(i,j)-V[k];residualSq+=err*err;cellCount++;}
    const gx=(V[id(Math.min(N-1,i+1),j,N)]-V[id(Math.max(0,i-1),j,N)])/
      ((i===0||i===N-1)?dx:2*dx);
    const gy=(V[id(i,Math.min(N-1,j+1),N)]-V[id(i,Math.max(0,j-1),N)])/
      ((j===0||j===N-1)?dx:2*dx);
    e.Ex[k]=-gx;e.Ey[k]=-gy;
    e.magnitude[k]=Math.hypot(gx,gy);
    maxField=Math.max(maxField,e.magnitude[k]);sumField+=e.magnitude[k];
  }
  e.diagnostics={residualRmsV:Math.sqrt(residualSq/Math.max(1,cellCount)),
    maxFieldVM:maxField,meanFieldVM:sumField/(N*N),
    maxPolarization:c.flexoCoefficient*c.prescribedGradient,sweeps};
  return e;
}

/**
 * Approximate fixed-field dielectric stabilization of a deposited ice
 * molecule. The replacement of ambient vapor by bulk ice is modeled as
 * Δg = -½ ε0 (ε_ice-ε_air) |E|² [J/m³].
 *
 * Divide by ice molecule number density ρ_ice/m and by kBT to get the
 * dimensionless ln(c_eq,electric/c_eq,neutral). This neglects Maxwell
 * traction, interfacial molecular orientation, conductivity, and charge
 * transport. Keep opt-in and report the correction independently.
 */
export function dielectricEquilibriumLogShift(state,cellIndex,rhoIce){
  const e=state.electric;
  if(!e?.params.enabled||!e.params.coupleGrowth)return 0;
  if(!(rhoIce>0)||!(state.Tinf>0))return 0;
  const dEps=e.params.epsilonIce-e.params.epsilonAir;
  const E=e.magnitude[cellIndex];
  const exponent=-0.5*VACUUM_PERMITTIVITY*dEps*E*E*
    WATER_MOLECULE_MASS/(rhoIce*BOLTZMANN*state.Tinf);
  // Mathematical safeguard; flagged by tests and documentation.
  return clamp(exponent,-0.10,0.10);
}
