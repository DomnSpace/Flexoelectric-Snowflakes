/**
 * Pass 7: versioned, non-mutating 2.5-D Observatory snapshot envelope.
 *
 * The solver owns phase/transport fields. This module describes them for
 * replay/rendering without inventing a locally varying thickness or changing
 * basal mass/heat accounts (deferred to Pass 8).
 */
export const ICE_SNAPSHOT_SCHEMA='snowflake-observatory/snapshot';
export const ICE_SNAPSHOT_VERSION=1;
export const UNIFORM_THICKNESS_MODE='legacy-uniform-slab';

function finiteNonnegative(value,label){
  if(!Number.isFinite(value)||value<0)throw new RangeError(label+' must be finite and nonnegative');
}
function validateFields(frame){
  const n=frame.n;
  if(!Number.isInteger(n)||n<3||n%2!==1)throw new RangeError('Invalid grid width');
  if(!Number.isFinite(frame.dx)||frame.dx<=0)throw new RangeError('Invalid grid spacing');
  finiteNonnegative(frame.t,'Model time');
  for(const field of ['phi','u','dT','ice']){
    if(frame[field]?.length!==n*n)throw new RangeError('Field '+field+' does not match grid');
  }
  if(frame.topology==='triangular'&&frame.mask?.length!==n*n)
    throw new RangeError('Hex snapshot must carry physical domain mask');
  if(!frame.morph?.radius?.length)throw new RangeError('Snapshot requires solved phase-field contour');
  finiteNonnegative(frame.morph.halfThickness,'Physical uniform half-thickness');
}

/** Capture a compact, explicit geometry descriptor; no n² height allocation. */
export function versionIceSnapshot(frame){
  validateFields(frame);
  if(frame.snapshotSchema||frame.iceStructure)throw new RangeError('Snapshot was already versioned');
  const half=frame.morph.halfThickness;
  return {
    ...frame,
    snapshotSchema:Object.freeze({
      name:ICE_SNAPSHOT_SCHEMA,
      version:ICE_SNAPSHOT_VERSION,
      topology:frame.topology,
      unitSystem:'SI',
      fields:Object.freeze({
        phase:'phi',vapor:'u',temperature:'dT',frozenMask:'ice',
        localThickness:null,firstFreezeTime:null
      })
    }),
    iceStructure:Object.freeze({
      representation:UNIFORM_THICKNESS_MODE,
      physicalHalfThicknessM:half,
      physicalTotalThicknessM:2*half,
      localHeightSolved:false,
      verticalVaporTransportSolved:false,
      basalMassCoupled:false,
      formationTimeSolved:false
    })
  };
}
export function verifyIceSnapshot(frame){
  validateFields(frame);
  if(frame.snapshotSchema?.name!==ICE_SNAPSHOT_SCHEMA||
    frame.snapshotSchema?.version!==ICE_SNAPSHOT_VERSION)
    throw new RangeError('Unknown or unversioned ice snapshot');
  if(frame.iceStructure?.representation!==UNIFORM_THICKNESS_MODE)
    throw new RangeError('Unknown ice thickness representation');
  if(frame.iceStructure.localHeightSolved!==false)
    throw new RangeError('Pass 7 cannot claim solved local thickness');
  finiteNonnegative(frame.iceStructure.physicalTotalThicknessM,'Uniform thickness');
  if(Math.abs(frame.iceStructure.physicalTotalThicknessM-
    2*frame.morph.halfThickness)>1e-15)
    throw new RangeError('Snapshot thickness differs from physical legacy solver');
  return true;
}
/** Explicit adapter for an old in-memory snap when replay compatibility is needed. */
export function upgradeLegacyIceSnapshot(frame){
  if(frame.snapshotSchema)return (verifyIceSnapshot(frame),frame);
  return versionIceSnapshot(frame);
}
/**
 * Read-only compatibility interface for Pass 8's eventual height arrays.
 * Ice phase occupancy is already modeled by phi. The original height remains
 * spatially uniform in occupied cells; neither a ridge nor height mass is
 * inferred from the interface gradient or formation time.
 */
export function legacyHeightAt(frame,index){
  verifyIceSnapshot(frame);
  if(!Number.isInteger(index)||index<0||index>=frame.n*frame.n)
    throw new RangeError('Invalid physical cell index');
  if(frame.mask&&!frame.mask[index])return {topM:0,bottomM:0,iceFraction:0};
  const fraction=Math.max(0,Math.min(1,frame.ice[index]?1:frame.phi[index]));
  const half=frame.iceStructure.physicalHalfThicknessM;
  return {topM:fraction>0?half:0,bottomM:fraction>0?-half:0,iceFraction:fraction};
}
