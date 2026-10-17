/**
 * Pass 9 — true local-height surfaces on the existing TRIANGULAR transport grid.
 * Entirely render-side: no crystals, mass, heat, vapor or steps are modified.
 *
 * Intersect each equilateral triangle with phi>=0.5, then either triangulate
 * interpolated top/bottom faces and construct footprint boundary walls, or
 * intersect it with the five horizontal planes hBottom<=z<=hTop.
 *
 * Interpolated phase is an optical interface approximation. Its triangulated
 * volume is NOT identical to the solver's rho*A_hex*sum(phi*h) finite-volume
 * mass; that ledger remains the sole physical ice-mass measure.
 *
 * Full triangles at 401² are expensive on phones; for N>161 use stride=2
 * render LOD. No smoothing and no ridge physics or generated ribs.
 */
import {verifyLocalIceSnapshot} from '../solver/ice-snapshot-v2.js';

export const PHASE_ISO=0.5;
export const HEIGHT_SLICE_LEVELS=Object.freeze([-.8,-.4,0,.4,.8]);
const SQ3=Math.sqrt(3);
const EPS=1e-18;
const interpolate=(a,b,t,plane)=>({
  x:a.x+(b.x-a.x)*t,y:a.y+(b.y-a.y)*t,
  phi:a.phi+(b.phi-a.phi)*t,
  top:a.top+(b.top-a.top)*t,
  bottom:a.bottom+(b.bottom-a.bottom)*t,
  cut:plane
});
function clipPolygon(poly,name,threshold,keepGreater){
  if(!poly.length)return [];
  const value=p=>p[name];
  const inside=p=>keepGreater?value(p)>=threshold:value(p)<=threshold;
  const out=[];
  for(let k=0;k<poly.length;k++){
    const a=poly[k],b=poly[(k+1)%poly.length];
    const ai=inside(a),bi=inside(b);
    if(ai)out.push(a);
    if(ai!==bi){
      const delta=value(b)-value(a);
      if(Math.abs(delta)<EPS)continue;
      const t=Math.max(0,Math.min(1,(threshold-value(a))/delta));
      // When clipping phase, z values outside ice have no defined physical
      // meaning; inherit the solid-side height at the cut boundary.
      const c=interpolate(a,b,t,name);
      if(name==='phi'){
        const solid=ai?a:b;
        c.top=solid.top;c.bottom=solid.bottom;
      }
      out.push(c);
    }
  }
  return out;
}
function area2(a,b,c){return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);}
function emit(out,a,b,c,scale=1e6,flip=false,flatZ=null){
  // Use signed orientation for TOP surfaces; bottom faces face down.
  let v1=b,v2=c;
  if((area2(a,b,c)<0)!==flip){v1=c;v2=b;}
  if(Math.abs(area2(a,b,c))<1e-24)return;
  for(const v of [a,v1,v2]){
    out.push(v.x*scale,v.y*scale,(flatZ===null?v.top:flatZ)*scale);
  }
}
function emitFace(out,poly,isTop,scale){
  if(poly.length<3)return 0;
  const a=poly[0],oriented=[];
  for(let i=1;i+1<poly.length;i++){
    const b=poly[i],c=poly[i+1];
    if(Math.abs(area2(a,b,c))<1e-24)continue;
    // The face emitter copies physical z from its top slot.
    const x={...a,top:isTop?a.top:a.bottom},
      y={...b,top:isTop?b.top:b.bottom},
      z={...c,top:isTop?c.top:c.bottom};
    emit(out,x,y,z,scale,!isTop);
    oriented.push(i);
  }
  return oriented.length;
}
function emitWall(out,a,b,scale){
  if(Math.hypot(b.x-a.x,b.y-a.y)<1e-12)return;
  if(Math.min(a.top-a.bottom,b.top-b.bottom)<=0)return;
  const aa={...a,top:a.top},bb={...b,top:b.top};
  const ab={...a,top:a.bottom},ba={...b,top:b.bottom};
  for(const [p,q,r] of [[aa,ab,bb],[bb,ab,ba]]){
    for(const v of [p,q,r])out.push(v.x*scale,v.y*scale,v.top*scale);
  }
}
function triangleCells(frame,stride,fn){
  const N=frame.n,center=(N-1)/2,a=frame.dx,root=SQ3/2;
  const node=(i,j)=>{
    const k=j*N+i;
    return {
      x:a*((i-center)+(j-center)/2),
      y:a*root*(j-center),phi:frame.ice[k]?1:frame.phi[k],
      top:frame.hTop[k],bottom:frame.hBottom[k],cut:null
    };
  };
  let baseTriangles=0;
  for(let j=0;j+stride<N;j+=stride)for(let i=0;i+stride<N;i+=stride){
    const k=j*N+i,kr=j*N+i+stride,
      ku=(j+stride)*N+i,kru=(j+stride)*N+i+stride;
    // A clipped triangle must never bridge outside the physical hex mask.
    const pairs=[[k,kr,ku],[kru,ku,kr]];
    let nodes=null;
    for(const [a0,b0,c0] of pairs){
      if(!frame.mask[a0]||!frame.mask[b0]||!frame.mask[c0])continue;
      nodes??=[node(i,j),node(i+stride,j),node(i,j+stride),
               node(i+stride,j+stride)];
      const tri=a0===k?[nodes[0],nodes[1],nodes[2]]:
        [nodes[3],nodes[2],nodes[1]];
      fn(tri);baseTriangles++;
    }
  }
  return baseTriangles;
}
export function selectHeightMeshStride(frame){
  return frame.n>161?2:1;
}
export function computeActualHeightGeometry(frame,{stride=selectHeightMeshStride(frame)}={}){
  verifyLocalIceSnapshot(frame);
  if(!Number.isInteger(stride)||stride<1||stride>4)throw new RangeError('Invalid rendering stride');
  const top=[],bottom=[],walls=[];
  let occupiedTriangles=0,wallSegments=0;
  const inspected=triangleCells(frame,stride,tri=>{
    const poly=clipPolygon(tri,'phi',PHASE_ISO,true);
    if(poly.length<3)return;
    occupiedTriangles+=emitFace(top,poly,true,1e6);
    emitFace(bottom,poly,false,1e6);
    // Only phi=iso boundary segments of a clipped triangle become side walls.
    // Internal triangle edges never create false "seams".
    for(let k=0;k<poly.length;k++){
      const a=poly[k],b=poly[(k+1)%poly.length];
      if(a.cut==='phi'&&b.cut==='phi'){
        emitWall(walls,a,b,1e6);wallSegments++;
      }
    }
  });
  return {
    representation:'local-height-triangulated-v1',stride,
    phaseIso:PHASE_ISO,
    top:new Float32Array(top),bottom:new Float32Array(bottom),
    walls:new Float32Array(walls),
    occupiedTriangles,wallSegments,inspectedTriangles:inspected,
    source:'recorded snapshot-v2 top/bottom/phase fields',
    volumeExactlyConserved:false
  };
}
export function computeActualHeightSlices(frame,{stride=selectHeightMeshStride(frame)}={}){
  verifyLocalIceSnapshot(frame);
  const half=Math.max(frame.iceStructure.maxLocalThicknessM/2,1e-12);
  const heights=HEIGHT_SLICE_LEVELS.map(f=>f*half);
  const vertices=heights.map(()=>[]),triangles=heights.map(()=>0),
    boundarySegments=heights.map(()=>0);
  const inspected=triangleCells(frame,stride,tri=>{
    let poly=clipPolygon(tri,'phi',PHASE_ISO,true);
    if(poly.length<3)return;
    for(let z=0;z<heights.length;z++){
      const level=heights[z];
      let p=clipPolygon(poly,'top',level,true);
      p=clipPolygon(p,'bottom',level,false);
      if(p.length<3)continue;
      const base=p[0],arr=vertices[z];
      for(let i=1;i+1<p.length;i++){
        const x={...base,top:level},y={...p[i],top:level},
          w={...p[i+1],top:level};
        if(Math.abs(area2(x,y,w))<1e-24)continue;
        emit(arr,x,y,w);
        triangles[z]++;
      }
      // Segment count is a diagnostic for clipped level sets, not a
      // topological loop count; no fake radial outline is generated.
      for(let i=0;i<p.length;i++){
        const a=p[i],b=p[(i+1)%p.length];
        if(a.cut&&a.cut===b.cut)boundarySegments[z]++;
      }
    }
  });
  return {
    representation:'five-physical-iso-z-cross-sections',
    levelM:heights,stride,inspectedTriangles:inspected,
    layers:vertices.map((v,i)=>({
      index:i,zM:heights[i],positions:new Float32Array(v),
      triangleCount:triangles[i],boundarySegments:boundarySegments[i]
    })),
    contourDerivedFromLocalHeight:true,
    includesDisconnectedRegions:true,
    preservesHolesWhereResolved:true,
    volumeExactlyConserved:false
  };
}
export function buildActualIceMeshes(THREE,frame,{slices=false}={}){
  const group=new THREE.Group();
  group.name=slices?'actual-five-isoheight-slices':'actual-local-ice-surfaces';
  const mkGeo=positions=>{
    const geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.BufferAttribute(positions,3));
    geo.computeVertexNormals();
    return geo;
  };
  if(slices){
    const plan=computeActualHeightSlices(frame);
    for(const layer of plan.layers){
      const mat=new THREE.MeshPhysicalMaterial({
        color:0xb5e4ff,roughness:.23,metalness:0,
        clearcoat:.68,transparent:true,opacity:.245,
        depthWrite:false,side:THREE.DoubleSide
      });
      const mesh=new THREE.Mesh(mkGeo(layer.positions),mat);
      mesh.name='actual-z-slice-'+(layer.index+1);
      mesh.userData.heightUm=layer.zM*1e6;
      mesh.userData.triangleCount=layer.triangleCount;
      group.add(mesh);
    }
    return {group,plan};
  }
  const plan=computeActualHeightGeometry(frame);
  for(const [name,vertices,color,opacity] of [
    ['top',plan.top,0xd5f0ff,1],
    ['bottom',plan.bottom,0x9ad7f5,1],
    ['rim',plan.walls,0x78b9dc,1]
  ]){
    const mat=new THREE.MeshPhysicalMaterial({
      color,roughness:.18,metalness:.03,clearcoat:.83,
      emissive:0x163449,emissiveIntensity:.06,
      side:THREE.DoubleSide,transparent:opacity<1,opacity
    });
    const mesh=new THREE.Mesh(mkGeo(vertices),mat);
    mesh.name='local-ice-'+name;
    group.add(mesh);
  }
  return {group,plan};
}
