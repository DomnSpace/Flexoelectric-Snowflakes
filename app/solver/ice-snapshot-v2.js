/**
 * Version 2 Observatory snapshot: an OPT-IN solver-owned local ice thickness
 * field, not a 3-D transport model. Legacy v1 snapshots remain unchanged.
 */
import {ICE_SNAPSHOT_SCHEMA} from './ice-snapshot-v1.js';
import {localHeightStats} from './ice-thickness.js';
export const ICE_SNAPSHOT_VERSION_LOCAL=2;

export function versionLocalIceSnapshot(frame,state,rhoIce){
  if(frame.snapshotSchema||frame.iceStructure)
    throw new RangeError('Local snapshot already versioned');
  if(!state.localHeight||!state.topHeight||state.topHeight.length!==frame.n*frame.n)
    throw new Error('Local height snapshot needs a running height-resolved solver');
  if(state.N!==frame.n||frame.topology!=='triangular')
    throw new Error('Local height coordinates must agree with triangular physical grid');
  const stats=localHeightStats(state,rhoIce);
  const upper=new Float32Array(state.topHeight),
    lower=new Float32Array(state.bottomHeight);
  const schema=Object.freeze({
    name:ICE_SNAPSHOT_SCHEMA,
    version:ICE_SNAPSHOT_VERSION_LOCAL,
    topology:'triangular',unitSystem:'SI',
    fields:Object.freeze({
      phase:'phi',vapor:'u',temperature:'dT',frozenMask:'ice',
      localUpperHeight:'hTop',localLowerHeight:'hBottom',
      localThickness:'hTop - hBottom',firstFreezeTime:null
    })
  });
  const structure=Object.freeze({
    representation:'locally-resolved-height-open-z-reservoir',
    physicalHalfThicknessM:stats.meanThicknessM/2,
    physicalTotalThicknessM:stats.meanThicknessM,
    maxLocalThicknessM:stats.maxThicknessM,
    minLocalThicknessM:stats.minThicknessM,
    localHeightSolved:true,verticalVaporTransportSolved:false,
    basalMassCoupled:true,formationTimeSolved:false,
    basalSupplyModel:'unresolved external vertical vapor reservoir',
    iceMassKg:stats.physicalMassKg,
    ledgerInitialKg:state.heightLedger.initialMassKg,
    ledgerPrismKg:state.heightLedger.prismInKg,
    ledgerBasalKg:state.heightLedger.basalVerticalInKg,
    ledgerResidualKg:state.heightLedger.lastResidualKg
  });
  return {...frame,hTop:upper,hBottom:lower,
    snapshotSchema:schema,iceStructure:structure};
}
export function verifyLocalIceSnapshot(frame){
  if(frame.snapshotSchema?.name!==ICE_SNAPSHOT_SCHEMA||
    frame.snapshotSchema.version!==ICE_SNAPSHOT_VERSION_LOCAL||
    frame.snapshotSchema.topology!=='triangular')
    throw new Error('Invalid local-height snapshot schema');
  if(frame.hTop?.length!==frame.n*frame.n||
    frame.hBottom?.length!==frame.n*frame.n||
    frame.phi?.length!==frame.n*frame.n)
    throw new Error('Local height arrays do not match the replay grid');
  if(frame.iceStructure?.localHeightSolved!==true||
    frame.iceStructure.verticalVaporTransportSolved!==false||
    frame.iceStructure.basalMassCoupled!==true)
    throw new Error('Local-height model capabilities are inconsistent');
  if(!Number.isFinite(frame.iceStructure.iceMassKg)||
    !Number.isFinite(frame.iceStructure.ledgerResidualKg))
    throw new Error('Nonfinite ice-side mass ledger');
  return true;
}
export function localThicknessAt(frame,index){
  verifyLocalIceSnapshot(frame);
  if(!Number.isInteger(index)||index<0||index>=frame.n*frame.n)
    throw new RangeError('Invalid local-height index');
  return frame.hTop[index]-frame.hBottom[index];
}
