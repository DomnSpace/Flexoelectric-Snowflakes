/**
 * Pass 7 five-sheet optical preview. This is NOT a resolved z field.
 * Every sheet has the same basal-plane phase-derived contour and is placed
 * at a physically sampled elevation in the current uniform-height slab.
 * Neither fake concentric ridges nor displaced XY arm tips are introduced.
 *
 * The pure plan is tested without Three.js; the WebGL adapter is render-only.
 */
import {verifyIceSnapshot} from '../solver/ice-snapshot-v1.js';
import {verifyLocalIceSnapshot} from '../solver/ice-snapshot-v2.js';

export const FIVE_LAYER_COUNT=5;
export const FIVE_LAYER_OPACITY=0.245;
export const FIVE_LAYER_Z_FRACTIONS=Object.freeze([-1,-0.5,0,0.5,1]);
const finite=(v)=>Number.isFinite(v);

/** Five planes: z=-h,-h/2,0,+h/2,+h. No invented height gradient. */
export function makeUniformLayerPlan(frame){
  const local=frame.snapshotSchema?.version===2;
  if(local)verifyLocalIceSnapshot(frame);
  else verifyIceSnapshot(frame);
  // Pass 9 will replace these five identical outlines by true iso-z loops.
  // When local heights exist, this is explicitly a MEAN-height optical guide.
  const h=frame.iceStructure.physicalHalfThicknessM;
  const radius=frame.morph.radius;
  if(radius.length<6||!Array.from(radius).every(v=>finite(v)&&v>=0))
    throw new RangeError('Invalid modeled basal-plane contour');
  return Object.freeze({
    mode:local?'mean-height-guide-not-local-slices':'uniform-five-layer-preview',
    representation:frame.iceStructure.representation,
    physicalThicknessM:2*h,
    sliceCount:FIVE_LAYER_COUNT,
    eachSliceOpacity:FIVE_LAYER_OPACITY,
    contourSource:'modeled-phase-isocontour',
    physicalLocalHeightSolved:local,
    sliceGeometryFromLocalHeight:false,
    layers:Object.freeze(FIVE_LAYER_Z_FRACTIONS.map((fraction,index)=>Object.freeze({
      index,zM:fraction*h,relativeZ:fraction,
      opacity:FIVE_LAYER_OPACITY,
      // Same input for every layer: no claims of simulated core ridges.
      contourRadii:radius
    })))
  });
}

export function buildUniformLayerGroup(THREE,frame,scaleMetersToScene=1e6){
  const plan=makeUniformLayerPlan(frame);
  const group=new THREE.Group();
  group.name='uniform-five-slice-preview';
  const shape=new THREE.Shape(),outline=[];
  for(let i=0;i<frame.morph.radius.length;i++){
    const angle=i*2*Math.PI/frame.morph.radius.length;
    const radius=scaleMetersToScene*frame.morph.radius[i];
    const x=radius*Math.cos(angle),y=radius*Math.sin(angle);
    if(i===0)shape.moveTo(x,y);else shape.lineTo(x,y);
    outline.push(new THREE.Vector3(x,y,0));
  }
  shape.closePath();
  const geometry=new THREE.ShapeGeometry(shape,1);
  const lineGeometry=new THREE.BufferGeometry().setFromPoints(outline);
  for(const layer of plan.layers){
    const material=new THREE.MeshPhysicalMaterial({
      color:0xc1eaff,roughness:.24,metalness:0,
      clearcoat:.7,clearcoatRoughness:.16,
      transparent:true,opacity:layer.opacity,depthWrite:false,
      side:THREE.DoubleSide,polygonOffset:true,
      polygonOffsetFactor:1,polygonOffsetUnits:1
    });
    const mesh=new THREE.Mesh(geometry,material);
    mesh.name='ice-slice-'+(layer.index+1);
    mesh.position.z=layer.zM*scaleMetersToScene;
    group.add(mesh);
    const line=new THREE.LineLoop(lineGeometry,new THREE.LineBasicMaterial({
      color:0x93dcfa,transparent:true,opacity:.18,depthWrite:false
    }));
    line.name='ice-slice-outline-'+(layer.index+1);
    line.position.z=mesh.position.z;
    group.add(line);
  }
  return {group,plan};
}

/** Dispose unique GPU geometries/materials; shared slice meshes are freed once. */
export function disposeLayerGroup(group){
  if(!group)return;
  const geometries=new Set(),materials=new Set();
  group.traverse(object=>{
    if(object.geometry)geometries.add(object.geometry);
    if(object.material){
      for(const m of Array.isArray(object.material)?object.material:[object.material])materials.add(m);
    }
  });
  for(const geo of geometries)geo.dispose();
  for(const mat of materials)mat.dispose();
  group.clear();
}
