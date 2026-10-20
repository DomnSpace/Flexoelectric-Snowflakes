/**
 * Pass 8 -- conservative ice geometry ledger for the experimental hex solver.
 *
 * The gas/heat PDE is still the pre-existing quasi-steady 2-D solver.
 * Basal growth is supplied by an EXPLICIT UNRESOLVED VERTICAL VAPOR RESERVOIR.
 * Local ice mass is conserved against admitted prism + basal reservoir fluxes;
 * this is NOT a closed vapor+ice conservation law or a solved z transport PDE.
 */
export const LOCAL_HEIGHT_MODEL='local-height-open-reservoir-v1';
const SQRT3_2=Math.sqrt(3)/2;
const clamp=(x,a,b)=>Math.min(b,Math.max(a,x));
const fraction=(state,k)=>state.ice[k]?1:clamp(state.phi[k],0,1);
export function localCellArea(state){return SQRT3_2*state.dx*state.dx;}

export function seedLocalHeight(state,rhoIce){
  if(!state.mask||!state.phi||state.topHeight||state.bottomHeight)
    throw new Error('Only a fresh hex growth state may initialize the local-height model');
  if(!(rhoIce>0))throw new RangeError('Invalid density');
  const n=state.N*state.N,top=new Float64Array(n),bottom=new Float64Array(n);
  const h=state.basalHalfThickness;
  for(let k=0;k<n;k++)if(state.mask[k]&&fraction(state,k)>0){
    top[k]=h;bottom[k]=-h;
  }
  state.topHeight=top;state.bottomHeight=bottom;
  const initial=localIceMass(state,rhoIce);
  state.heightLedger={
    model:LOCAL_HEIGHT_MODEL,
    initialMassKg:initial,
    prismInKg:0,basalVerticalInKg:0,
    lastPrismKg:0,lastBasalKg:0,
    lastResidualKg:0,relativeResidual:0,
    source:'open quasi-steady lateral vapor + unresolved vertical basal reservoir'
  };
  return state;
}

export function localIceMass(state,rhoIce){
  if(!state.topHeight||!state.bottomHeight)throw new Error('Local thickness not initialized');
  const area=localCellArea(state);
  let volume=0;
  for(let k=0;k<state.phi.length;k++)if(state.mask[k]){
    const h=state.topHeight[k]-state.bottomHeight[k];
    if(h<0||!Number.isFinite(h))throw new RangeError('Negative or nonfinite ice thickness');
    volume+=fraction(state,k)*h;
  }
  return volume*area*rhoIce;
}

export function localHeightStats(state,rhoIce){
  const area=localCellArea(state),top=state.topHeight,bottom=state.bottomHeight;
  if(!top||!bottom)throw new Error('Local thickness not initialized');
  let maxThickness=0,minThickness=Infinity,weighted=0,fractions=0,sites=0;
  for(let k=0;k<state.phi.length;k++){
    if(!state.mask[k])continue;
    const phi=fraction(state,k);
    if(phi<=0)continue;
    const depth=top[k]-bottom[k];
    if(!Number.isFinite(depth)||depth<0)throw new RangeError('Invalid height');
    maxThickness=Math.max(maxThickness,depth);
    minThickness=Math.min(minThickness,depth);
    weighted+=phi*depth;fractions+=phi;sites++;
  }
  return {
    physicalMassKg:localIceMass(state,rhoIce),
    minThicknessM:sites?minThickness:0,
    maxThicknessM:maxThickness,
    meanThicknessM:fractions?weighted/fractions:0,
    occupiedCells:sites,iceFootprintM2:fractions*area,
    localHeightSolved:true,verticalVaporTransportSolved:false
  };
}

/**
 * Prism accretion adds occupied in-plane volume at the inherited physical
 * local depth, capped to available φ. Incoming mass is never credited twice.
 */
export function acceptPrismMass(state,k,offeredKg,rhoIce,inheritTop,inheritBottom){
  if(!state.mask[k]||state.ice[k]||state.phi[k]>=1)return 0;
  if(!(offeredKg>=0)||!Number.isFinite(offeredKg))
    throw new RangeError('Invalid prism mass flux');
  const top=state.topHeight,bottom=state.bottomHeight;
  const prior=clamp(state.phi[k],0,1);
  if(prior===0){
    if(!(inheritTop>inheritBottom))throw new RangeError('Missing local facet thickness');
    top[k]=inheritTop;bottom[k]=inheritBottom;
  }
  const depth=top[k]-bottom[k];
  if(!(depth>0))throw new RangeError('Growth requires positive inherited thickness');
  const maxMass=rhoIce*localCellArea(state)*depth*(1-prior);
  const admitted=Math.min(offeredKg,maxMass);
  const next=prior+admitted/(rhoIce*localCellArea(state)*depth);
  state.phi[k]=next>=1-1e-12?1:next;
  if(state.phi[k]===1)state.ice[k]=1;
  return admitted;
}

/** Exposed basal surfaces are supplied by an explicit external z reservoir. */
export function acceptBasalMass(state,k,fluxKgM2S,dt,rhoIce){
  const phi=fraction(state,k);
  if(!state.mask[k]||phi===0)return 0;
  if(!Number.isFinite(fluxKgM2S)||fluxKgM2S<0)throw new RangeError('Invalid basal flux');
  // Numerical guard, not an ice-material upper bound. If reached, the
  // admitted mass is reduced consistently with the actual geometry change.
  const eachFaceStep=Math.min(fluxKgM2S*dt/rhoIce,state.dx*.1);
  state.topHeight[k]+=eachFaceStep;
  state.bottomHeight[k]-=eachFaceStep;
  return rhoIce*localCellArea(state)*phi*(2*eachFaceStep);
}

/**
 * Strict ice-side finite-volume accounting:
 * M(t)-M(seed) = M_from_xy_reservoir + M_from_unresolved_z_reservoir.
 * Neither reservoir has its own finite-volume conservation equation yet.
 */
export function finishLocalHeightStep(state,{rhoIce,Ls,dt,prismKg,basalKg}){
  const ledger=state.heightLedger;
  if(!ledger||!(rhoIce>0&&Ls>0&&dt>0))
    throw new Error('Missing conservative ice-side ledger');
  ledger.prismInKg+=prismKg;
  ledger.basalVerticalInKg+=basalKg;
  ledger.lastPrismKg=prismKg;ledger.lastBasalKg=basalKg;
  const actual=localIceMass(state,rhoIce);
  const predicted=ledger.initialMassKg+ledger.prismInKg+ledger.basalVerticalInKg;
  const error=actual-predicted;
  ledger.lastResidualKg=error;
  ledger.relativeResidual=Math.abs(error)/Math.max(1e-24,predicted);
  if(ledger.relativeResidual>2e-10)
    throw new Error('Ice mass ledger drift: '+ledger.relativeResidual);
  const result=localHeightStats(state,rhoIce);
  state.basalHalfThickness=result.meanThicknessM/2; // diagnostic proxy ONLY
  state.totalDepositedMass+=prismKg+basalKg;
  state.totalLatentHeat+=(prismKg+basalKg)*Ls;
  return {...result,prismKg,basalKg,totalAddedKg:prismKg+basalKg,
    lastLatentJ:(prismKg+basalKg)*Ls,ledgerResidualKg:error,
    relativeLedgerResidual:ledger.relativeResidual};
}

/** Transfer the whole reduced-order mass state to a larger same-dx hex domain. */
export function extendLocalHeight(oldState,nextState,offset){
  if(!oldState.topHeight)throw new Error('Old state has no local-height physics');
  if(nextState.topHeight)throw new Error('Destination should be uninitialized');
  const n=nextState.N*nextState.N,top=new Float64Array(n),bottom=new Float64Array(n);
  for(let j=0;j<oldState.N;j++)for(let i=0;i<oldState.N;i++){
    const old=j*oldState.N+i;
    if(!oldState.mask[old])continue;
    const k=(j+offset)*nextState.N+i+offset;
    top[k]=oldState.topHeight[old];bottom[k]=oldState.bottomHeight[old];
  }
  nextState.topHeight=top;nextState.bottomHeight=bottom;
  nextState.heightLedger={...oldState.heightLedger};
  nextState.localHeight=true;
  return nextState;
}
