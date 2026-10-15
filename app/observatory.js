import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import {
  makeGrowth25DState, resolvedGrowthStep, growthMorphology,
  DEFAULT_ICE_KINETICS
} from './solver/growth25d.js';
import {GROWTH_HABITS, GRID_CHOICES, auditGrowthDomain, modelDomainUm} from './solver/habits25d.js';
import {saturationVaporPressureWaterMK, saturationVaporPressureIceMK} from './solver/thermo.js';
import {makeHexGrowthState,resolvedHexGrowthStep,hexGrowthMorphology,
  auditHexDomain,enlargeHexDomain,hexCellAtXY} from './solver/hex25d.js';

/*
 * Snowflake Observatory
 * Truth: app/solver/growth25d.js. Every curve, slice and crystal outline comes
 * from that same resolved phase field. The scene is a 3-D extrusion of a 2-D
 * contour, not a claim of 3-D vapor or heat transport.
 *
 * The fixed timestep, deterministic initial state and stored snapshots make
 * visual replay independent from frame rate.
 */
const $ = id => document.getElementById(id);
const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
const num = (x, n = 2) => Number.isFinite(x) ? x.toFixed(n) : '—';
const um = x => x * 1e6;
const pg = x => x * 1e12;
const STEPTIME = 0.5;
const BASE_STEPS_PER_SECOND = 3;
const MAX_SNAPSHOTS = 360;
const MAX_STORED_GRID_CELLS=360*73*73; // bounded replay on larger grids
const DX=0.75e-6;
const ids = ['temp','sigma','pressure','anisotropy','prism','basalA','seed','gridSize','sdak','capillarity','electricEnabled','electricCouple','eField','eAngle','epsIce','bendGradient','bendAngle','muFlexo'];
const labels = [
  ['temp','tempValue',v => num(v,0)+' °C'],
  ['sigma','sigmaValue',v => num(v,1)+' %'],
  ['pressure','pressureValue',v => num(v,0)+' kPa'],
  ['anisotropy','anisotropyValue',v => num(v,2)],
  ['prism','prismValue',v => num(v,4)],
  ['seed','seedValue',v => num(v,2)+' µm'],
  ['basalA','basalAValue',v => num(v,3)],
  ['gridSize','gridSizeValue',v => v+' × '+v+' · '+num(modelDomainUm(v,DX),1)+' µm'],
  ['eField','eFieldValue',v => num(v,2)+' MV/m'],
  ['eAngle','eAngleValue',v => num(v,0)+'°'],
  ['epsIce','epsIceValue',v => num(v,1)],
  ['bendGradient','bendGradientValue',v => num(v,0)+' m⁻¹'],
  ['bendAngle','bendAngleValue',v => num(v,0)+'°'],
  ['muFlexo','muFlexoValue',v => num(v,2)+' nC/m']
];
const palette = ['#8de0ff','#b7a9ff','#f9cb8b','#72d8b7','#ff9da5','#e8efff'];
let state, frames=[], selectedIndex=0, replay=false, playing=true, pending=false;
let activeHabit='custom';
const isHex=()=>state?.topology==='triangular';
const auditState=()=>isHex()?auditHexDomain(state):auditGrowthDomain(state);
let lastFrame=performance.now(), accumulator=0, mesh=null, edge=null, worldScale=1;
let scene, camera, controls, renderer, crystal, snowLight;
let drawCanvas;

function config() {
  const v = id => Number($(id).value);
  return {
    N: v('gridSize'), dx: DX, Tinf: 273.15 + v('temp'),
    sigmaInf: v('sigma')/100, pressurePa: v('pressure')*1000,
    seedRadius: v('seed')*1e-6, basalHalfThickness: 1.5e-6,
    kinetics: {
      prism: {mode:'nucleation',A:1,sigma0:v('prism')},
      basal: {...DEFAULT_ICE_KINETICS.basal,A:v('basalA')},
      anisotropy6: v('anisotropy'), sdak:$('sdak').checked,
      sdakStrength:0.55, capillarity:$('capillarity').checked,
      gammaIv:DEFAULT_ICE_KINETICS.gammaIv
    },
    electric:{
      enabled:$('electricEnabled').checked,
      fieldVM:v('eField')*1e6,angleDeg:v('eAngle'),
      epsilonIce:v('epsIce'),epsilonAir:1.0006,
      prescribedGradient:v('bendGradient'),gradientAngleDeg:v('bendAngle'),
      flexoCoefficient:v('muFlexo')*1e-9,
      coupleGrowth:$('electricCouple').checked
    }
  };
}
function cloudSaturationSigma(tC){
  const T=273.15+tC;
  return saturationVaporPressureWaterMK(T)/saturationVaporPressureIceMK(T)-1;
}
function formatControls() {
  for(const [id,out,fmt] of labels) $(out).textContent=fmt(Number($(id).value));
  const cap=cloudSaturationSigma(Number($('temp').value));
  const above=Number($('sigma').value)/100>cap+1e-7;
  $('waterLimit').textContent=num(cap*100,1)+' % above ice'+(above?' · ABOVE WATER SATURATION':'');
  $('waterLimit').style.color=above?'#ffcb98':'';
}
function setStatus(message, error=false) {
  $('status').textContent=message;
  $('status').style.color=error?'#ffb6a7':'';
}
function setPending(v) {
  pending=v; $('pending').hidden=!v;
}
function snap() {
  const m=isHex()?hexGrowthMorphology(state,180):growthMorphology(state,180);
  return {
    t:state.t, n:state.N, dx:state.dx, Tinf:state.Tinf,
    topology:state.topology||'square',mask:isHex()?new Uint8Array(state.mask):null,
    expansions:state.domainExpansions||0,
    u:new Float32Array(state.sigma), dT:new Float32Array(state.dT),
    phi:new Float32Array(state.phi), ice:new Uint8Array(state.ice),
    ePotential:state.electric?new Float32Array(state.electric.potential):null,
    eMagnitude:state.electric?new Float32Array(state.electric.magnitude):null,
    eDiagnostics:state.electric?{...state.electric.diagnostics}:null,
    eEnabled:!!state.electric?.params.enabled,
    eCoupled:!!state.electric?.params.coupleGrowth,
    morph:{
      radius:m.radius, meanRadius:m.meanRadius, equivalentBasalRadius:m.equivalentBasalRadius,
      halfThickness:m.halfThickness, aspectRatio:m.aspectRatio,
      modes:{4:{...m.modes[4]},6:{...m.modes[6]},12:{...m.modes[12]}}
    },
    budget:state.lastBudget?{...state.lastBudget}:null,
    audit:auditState(),habit:activeHabit,
    waterSaturationPercent:cloudSaturationSigma(state.Tinf-273.15)*100,
    deposited:state.totalDepositedMass, energy:state.totalLatentHeat
  };
}
function current() {return frames[clamp(selectedIndex,0,frames.length-1)];}
function capture() {
  frames.push(snap());
  const historyLimit=Math.max(40,Math.min(MAX_SNAPSHOTS,Math.floor(MAX_STORED_GRID_CELLS/(state.N*state.N))));
  if(frames.length>historyLimit)frames.shift();
  selectedIndex=frames.length-1;
  $('scrub').max=String(frames.length-1);
  $('scrub').value=String(selectedIndex);
}
function reset(start=true) {
  try {
    const next=config();
    if(!GRID_CHOICES.includes(next.N))throw new RangeError('Grid size is not an approved resolution');
    if($('topology').value==='triangular'){
      if(next.electric.enabled||next.electric.coupleGrowth)
        throw new Error('Hexagonal experimental mesh does not support the electric Poisson model yet. Switch to Cartesian.');
      state=makeHexGrowthState(next);
    }else state=makeGrowth25DState(next);
    frames=[];replay=false;playing=start;accumulator=0;
    $('errorOverlay').hidden=true;
    setPending(false);
    capture();sync();setStatus(start?'SOLVER RUNNING':'SOLVER PAUSED');
    const selected=GROWTH_HABITS[activeHabit];
    $('scenario').textContent=(selected?.supported?selected.name.toUpperCase():'CUSTOM EXPERIMENT')+
      ' · '+num(Number($('temp').value),0)+' °C / '+num(Number($('sigma').value),1)+' %';
    $('electricReadout').textContent=isHex()?'NOT COUPLED TO TRIANGULAR MESH':
      state.electric?(state.electric.params.coupleGrowth?'POISSON + DIELECTRIC GROWTH':'POISSON / DIAGNOSTIC ONLY'):'DISABLED';
    $('play').textContent=start?'❚❚ Pause':'▶ Play';
    $('live').classList.add('active');
  } catch(e) { fail(e); }
}
function fail(e) {
  playing=false; $('play').textContent='▶ Play';
  setStatus('SOLVER ERROR',true);
  $('errorOverlay').hidden=false;
  $('errorOverlay').textContent='Physics or rendering failed: '+(e?.message||String(e));
  console.error(e);
}
function advance() {
  if(!state)return;
  try{
    // Preserve the grown ice, time, vapor fields and cumulative budgets when
    // expanding the experimental hex mesh (no restart or fake second seed).
    if(isHex()&&$('autoExpand').checked){
      const options=GRID_CHOICES, index=options.indexOf(state.N);
      if(index>=0&&index<options.length-1&&auditHexDomain(state).haloCells<=12){
        state=enlargeHexDomain(state,options[index+1]);
        $('gridSize').value=String(state.N);
        formatControls();
        $('fieldInspect').textContent='HEX DOMAIN GROWN · '+state.domainExpansions+' EXPANSIONS';
      }
    }
    const audit=auditState();
    if(audit.stop){
    playing=false;$('play').textContent='▶ Play';
    setStatus('MODEL DOMAIN STOP · '+audit.reason);
    return;
    }
    // Each model frame advances exactly one resolved transport + deposition step.
    if(isHex())resolvedHexGrowthStep(state,STEPTIME,{vaporIterations:45,heatIterations:55});
    else resolvedGrowthStep(state,STEPTIME,{vaporIterations:45,heatIterations:55});
    capture();sync();
    const d=current().audit;
    if(d.stop){
      playing=false;$('play').textContent='▶ Play';
      setStatus('MODEL DOMAIN STOP · '+d.reason);
    }else if(d.warning){
      setStatus('MODEL VALIDITY WARNING · '+d.reason);
    }else setStatus(playing?'SOLVER RUNNING':'SOLVER PAUSED');
  } catch(e) {fail(e);}
}
function play() {
  if(state&&auditState().stop&&!(isHex()&&$('autoExpand').checked&&GRID_CHOICES.indexOf(state.N)<GRID_CHOICES.length-1)){
    playing=false;setStatus('MODEL DOMAIN STOP · CHANGE DOMAIN OR PARAMETERS AND RESTART');return;
  }
  if(replay){replay=false;selectedIndex=frames.length-1;$('scrub').value=String(selectedIndex);}
  playing=!playing;accumulator=0;
  $('play').textContent=playing?'❚❚ Pause':'▶ Play';
  setStatus(playing?'SOLVER RUNNING':'SOLVER PAUSED');
  $('live').classList.toggle('active',!replay);
  sync();
}
function sync() {
  const s=current();if(!s)return;
  $('kTime').innerHTML=num(s.t,1)+' <small>s</small>';
  $('kRadius').innerHTML=num(um(s.morph.equivalentBasalRadius),2)+' <small>µm</small>';
  $('kMass').innerHTML=num(pg(s.deposited),4)+' <small>pg</small>';
  $('kHeat').innerHTML=num((s.budget?.maxDeltaT||0)*1e3,2)+' <small>mK</small>';
  $('kSix').innerHTML=num((s.morph.modes[6].relative||0)*100,2)+' <small>%</small>';
  $('sizeLabel').textContent=num((s.n-1)*s.dx*1e6,1)+' µm DOMAIN · '+s.n+'²';
  $('stepLabel').textContent='STEP '+String(Math.round(s.t/STEPTIME)).padStart(4,'0');
  $('scrubTime').textContent='t = '+num(s.t,1)+' s';
  $('gridDim').textContent=s.topology==='triangular'
    ?Math.round(1+3*((s.n-1)/2)*((s.n+1)/2))+' HEX CELLS':s.n+'×'+s.n;
  $('fieldGeometry').textContent=s.topology==='triangular'
    ?'TRIANGULAR TRANSPORT / HEX VORONOI · FINITE HEX BOUNDARY SELECTS m6':
    'SQUARE DIRICHLET FAR-FIELD · ARTIFACTS NEAR EDGES';
  $('topologyValue').textContent=s.topology==='triangular'?'6-neighbor triangular':'9-point Cartesian';
  $('live').classList.toggle('active',!replay);
  updateBudget(s);
  updateDomain(s);
  updateMesh(s);
  drawField(s);
  drawAllGraphs(s);
}
function updateDomain(s) {
  const d=s.audit;
  const label=d.stop?'STOPPED — OUTSIDE VALIDITY':d.warning?'CAUTION — APPROXIMATION': 'VALIDITY MARGIN AVAILABLE';
  $('domainStatus').textContent=label;
  $('domainDetails').textContent='TIP '+num(d.tipRadiusUm,1)+' µm · HALO '+num(d.haloUm,1)+
    ' µm ('+d.haloCells+' cells) · THICKNESS '+num(d.thicknessUm,1)+
    ' µm · h/r '+num(d.plateRatio,2)+' · '+d.reason+(s.topology==='triangular'?' · EXTENSIONS '+s.expansions:'');
  $('domainReport').classList.toggle('domain-stop',d.stop);
  $('domainReport').classList.toggle('domain-warning',d.warning&&!d.stop);
  $('validity').textContent=d.stop?'STOP · INVALID':d.warning?'WARNING':'MONITORING';
  $('validity').style.color=d.stop?'#ffa0a0':d.warning?'#ffce91':'';
}
function updateBudget(s) {
  const b=s.budget;
  const fields={
    bSigma:b?num(100*b.meanSigma,3)+' %':'—',
    bAlpha:b?num(b.meanPrismAlpha,4):'—',
    bBasal:b?num(um(b.vBasal),5)+' µm/s':'—',
    bCells:b?String(b.interfaceCells):'—',
    bCurve:b?num(b.meanCurvature*1e-6,4)+' µm⁻¹':'—',
    bCap:b?num(b.meanCapillaryExponent,6):'—',
    bFlux:b?num(pg(b.totalMassRate),5)+' pg/s':'—',
    bEField:s.eEnabled?num((s.eDiagnostics?.maxFieldVM||0)*1e-6,4)+' MV/m':'OFF',
    bPolar:s.eEnabled?num((s.eDiagnostics?.maxPolarization||0)*1e6,4)+' µC/m²':'OFF',
    bElectricShift:b?num((b.meanElectricLogShift||0)*1e6,5)+' ppm':'—',
    bEResidual:s.eEnabled?num(s.eDiagnostics?.residualRmsV||0,7)+' V':'OFF'
  };
  for(const [id,value] of Object.entries(fields)) $(id).textContent=value;
}

function initializeScene() {
  const viewport=$('viewport');
  renderer=new THREE.WebGLRenderer({antialias:true,preserveDrawingBuffer:true,alpha:true,powerPreference:'high-performance'});
  renderer.setPixelRatio(Math.min(window.devicePixelRatio||1,2));
  renderer.outputColorSpace=THREE.SRGBColorSpace;
  renderer.toneMapping=THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure=1.55;
  viewport.insertBefore(renderer.domElement,viewport.firstChild);
  scene=new THREE.Scene();
  scene.fog=new THREE.FogExp2(0x0b1d2b,0.013);
  camera=new THREE.PerspectiveCamera(39,1,0.05,300);
  camera.up.set(0,0,1);
  camera.position.set(5.8,-9.1,9.4);
  camera.lookAt(0,0,0);
  controls=new OrbitControls(camera,renderer.domElement);
  controls.enableDamping=true;controls.dampingFactor=.07;
  controls.enablePan=false;controls.minDistance=3;controls.maxDistance=32;
  controls.autoRotate=$('rotate').checked;controls.autoRotateSpeed=.55;
  scene.add(new THREE.HemisphereLight(0xd2edff,0x123247,2.2));
  snowLight=new THREE.DirectionalLight(0xd8f2ff,3.3);
  snowLight.position.set(-5,3,8);scene.add(snowLight);
  const fill=new THREE.PointLight(0x6cc8ff,95,27,2);
  fill.position.set(4,-5,-3);scene.add(fill);
  const rim=new THREE.DirectionalLight(0xa38dff,2);
  rim.position.set(3,5,-4);scene.add(rim);
  crystal=new THREE.Group();scene.add(crystal);
  const guide=new THREE.GridHelper(8,16,0x3b708c,0x1d3a4e);
  guide.rotation.x=Math.PI/2;guide.position.z=-0.12;
  guide.material.transparent=true;guide.material.opacity=.19;
  scene.add(guide);
  const ringPoints=[];
  for(let j=1;j<=3;j++){
    for(let i=0;i<=160;i++){
      const a=i/160*Math.PI*2,r=j*1.06;
      ringPoints.push(new THREE.Vector3(r*Math.cos(a),r*Math.sin(a),-.07));
    }
    const ring=new THREE.Line(new THREE.BufferGeometry().setFromPoints(ringPoints.splice(0)),new THREE.LineBasicMaterial({color:0x396882,transparent:true,opacity:.24}));
    scene.add(ring);
  }
  // Subtle non-data star field solely for spatial contrast.
  const particles=[];
  let seed=53117;
  const rand=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;};
  for(let i=0;i<450;i++){
    const th=rand()*Math.PI*2,zz=(rand()-.5)*2,r=15+rand()*24,rr=Math.sqrt(1-zz*zz);
    particles.push(r*rr*Math.cos(th),r*rr*Math.sin(th),r*zz);
  }
  const stars=new THREE.BufferGeometry();
  stars.setAttribute('position',new THREE.Float32BufferAttribute(particles,3));
  scene.add(new THREE.Points(stars,new THREE.PointsMaterial({color:0x7ca9c2,size:.035,transparent:true,opacity:.55,sizeAttenuation:true})));
  function resize(){
    const w=viewport.clientWidth,h=viewport.clientHeight;
    if(!w||!h)return;
    camera.aspect=w/h;camera.updateProjectionMatrix();renderer.setSize(w,h,false);
  }
  new ResizeObserver(resize).observe(viewport);
  resize();
  viewport.addEventListener('dblclick',viewIso);
}
function viewIso(){
  camera.position.set(5.8,-9.1,9.4);
  camera.up.set(0,0,1);
  controls.target.set(0,0,0);controls.update();
  $('cameraLabel').textContent='ISOMETRIC';
}
function viewTop(){
  camera.position.set(.001,-.001,11);
  camera.up.set(0,1,0);controls.target.set(0,0,0);controls.update();
  $('cameraLabel').textContent='TOP-DOWN';
}
function updateMesh(s) {
  if(!crystal)return;
  if(mesh){crystal.remove(mesh);mesh.geometry.dispose();mesh.material.dispose();}
  if(edge){crystal.remove(edge);edge.geometry.dispose();edge.material.dispose();}
  const radii=s.morph.radius;
  const n=radii.length;
  const shape=new THREE.Shape();
  for(let k=0;k<n;k++){
    const a=2*Math.PI*k/n, r=um(radii[k]);
    const x=r*Math.cos(a),y=r*Math.sin(a);
    if(k===0)shape.moveTo(x,y);else shape.lineTo(x,y);
  }
  shape.closePath();
  const depth=Math.max(.15,um(2*s.morph.halfThickness));
  const bevel=Math.min(.18,Math.max(.035,depth*.08));
  const geometry=new THREE.ExtrudeGeometry(shape,{
    depth:Math.max(.05,depth-2*bevel),steps:1,bevelEnabled:true,
    bevelThickness:bevel,bevelSize:Math.min(.12,bevel),bevelSegments:2,
    curveSegments:1
  });
  geometry.computeVertexNormals();
  const material=new THREE.MeshPhysicalMaterial({
    color:0xc6ecff,metalness:.16,roughness:.24,
    clearcoat:.85,clearcoatRoughness:.15,
    emissive:0x163649,emissiveIntensity:.22,
    side:THREE.DoubleSide,flatShading:false
  });
  mesh=new THREE.Mesh(geometry,material);
  mesh.position.z=-depth/2;
  crystal.add(mesh);
  const lineGeo=new THREE.EdgesGeometry(geometry,27);
  edge=new THREE.LineSegments(lineGeo,new THREE.LineBasicMaterial({color:0x8ddfff,transparent:true,opacity:.38}));
  edge.position.z=mesh.position.z;crystal.add(edge);
  // A common scale for the grid; adapt slowly enough to avoid zoomed-out seed.
  const maxRadius=Math.max(3.2,um(s.morph.equivalentBasalRadius)*1.5);
  worldScale=2.5/maxRadius;
  crystal.scale.setScalar(worldScale);
  $('viewLabel').textContent=(s.topology==='triangular'?'TRIANGULAR-LATTICE · ':'')+
    '2.5-D EXTRUSION · '+num(depth,2)+' µm';
}

function fieldColor(t,ice=false) {
  if(ice)return [201,238,249];
  const stops=[[14,35,61],[22,81,121],[40,145,187],[112,207,225],[244,248,247]];
  const p=clamp(t,0,1)*(stops.length-1),a=Math.min(stops.length-2,Math.floor(p)),f=p-a;
  return [0,1,2].map(c=>Math.round(stops[a][c]*(1-f)+stops[a+1][c]*f));
}
function fieldValues(s,kind) {
  if(kind==='temperature')return {values:s.dT,unit:'mK',scale:1000};
  if(kind==='occupancy')return {values:s.phi,unit:'φ',scale:1};
  if(kind==='potential')return {values:s.ePotential||new Float32Array(s.n*s.n),unit:'V',scale:1};
  if(kind==='electric')return {values:s.eMagnitude||new Float32Array(s.n*s.n),unit:'MV/m',scale:1e-6};
  return {values:s.u,unit:'%',scale:100};
}
function hexRasterGeometry(s,can){
  const R=(s.n-1)/2;
  const units=Math.min((can.width-32)/(2*R+2),(can.height-32)/(Math.sqrt(3)*R+2));
  return {unit:units,cx:can.width/2,cy:can.height/2,R};
}
function drawHexField(s,ctx,can,kind,values,scale,unit){
  can.width=600;can.height=600;
  const g=hexRasterGeometry(s,can),N=s.n;
  const visible=[];
  for(let k=0;k<N*N;k++)if(s.mask[k]&&
    (kind==='occupancy'||!s.ice[k]))visible.push(kind==='occupancy'&&s.ice[k]?1:values[k]);
  let lo=kind==='occupancy'?0:Math.min(...visible);
  let hi=kind==='occupancy'?1:Math.max(...visible);
  if(!Number.isFinite(lo)||!Number.isFinite(hi)){lo=0;hi=1;}
  const den=Math.max(1e-10,hi-lo);
  ctx.fillStyle='#071a28';ctx.fillRect(0,0,can.width,can.height);
  // Draw actual Voronoi hexagons of triangular transport nodes, not a
  // square heatmap with sixfold color applied after the solve.
  const rad=g.unit/Math.sqrt(3)+.28;
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    const k=j*N+i;if(!s.mask[k])continue;
    const q=i-g.R,r=j-g.R;
    const x=g.cx+g.unit*(q+r/2),y=g.cy-g.unit*(Math.sqrt(3)/2*r);
    const value=kind==='occupancy'&&s.ice[k]?1:values[k];
    const rgb=fieldColor((value-lo)/den,['vapor','temperature'].includes(kind)&&s.ice[k]);
    ctx.fillStyle='rgb('+rgb.join(',')+')';
    ctx.beginPath();
    for(let t=0;t<6;t++){
      const th=Math.PI/6+t*Math.PI/3,px=x+rad*Math.cos(th),py=y+rad*Math.sin(th);
      if(t===0)ctx.moveTo(px,py);else ctx.lineTo(px,py);
    }
    ctx.closePath();ctx.fill();
  }
  $('fieldMin').textContent=num(lo*scale,3)+' '+unit;
  $('fieldMax').textContent=num(hi*scale,3)+' '+unit;
}
function drawField(s) {
  const can=$('fieldCanvas'),ctx=can.getContext('2d'),kind=$('field').value;
  if(!ctx)return;
  const N=s.n,{values,scale,unit}=fieldValues(s,kind);
  if(s.topology==='triangular'){drawHexField(s,ctx,can,kind,values,scale,unit);return;}
  const data=[];
  for(let k=0;k<N*N;k++){
    if(['occupancy','electric','potential'].includes(kind))data.push(kind==='occupancy'&&s.ice[k]?1:values[k]);
    else if(!s.ice[k])data.push(values[k]);
  }
  let lo=kind==='occupancy'?0:Math.min(...data),hi=kind==='occupancy'?1:Math.max(...data);
  if(!Number.isFinite(lo)||!Number.isFinite(hi)){lo=0;hi=1;}
  const den=Math.max(1e-10,hi-lo);
  can.width=N;can.height=N;
  const image=ctx.createImageData(N,N);
  for(let k=0;k<N*N;k++){
    const fraction=kind==='occupancy'?(s.ice[k]?1:values[k]):values[k];
    const rgb=fieldColor((fraction-lo)/den,['vapor','temperature'].includes(kind)&&s.ice[k]);
    const o=k*4;image.data[o]=rgb[0];image.data[o+1]=rgb[1];image.data[o+2]=rgb[2];image.data[o+3]=255;
  }
  ctx.putImageData(image,0,0);
  const decimal=['electric','temperature'].includes(kind)?3:2;
  $('fieldMin').textContent=num(lo*scale,decimal)+' '+unit;
  $('fieldMax').textContent=num(hi*scale,decimal)+' '+unit;
}
function onFieldMove(evt) {
  const s=current();if(!s)return;
  const r=evt.currentTarget.getBoundingClientRect(),N=s.n;
  let i,j,k,x,y;
  if(s.topology==='triangular'){
    const can=evt.currentTarget,g=hexRasterGeometry(s,can);
    const px=(evt.clientX-r.left)/r.width*can.width,py=(evt.clientY-r.top)/r.height*can.height;
    x=(px-g.cx)/g.unit*s.dx;y=-(py-g.cy)/g.unit*s.dx;
    k=hexCellAtXY({N:s.n,R:g.R,dx:s.dx,mask:s.mask},x,y);
    if(k<0){$('fieldInspect').textContent='OUTSIDE HEXAGONAL TRANSPORT DOMAIN';return;}
    i=k%N;j=Math.floor(k/N);
  }else{
    i=clamp(Math.floor((evt.clientX-r.left)/r.width*N),0,N-1);
    j=clamp(Math.floor((evt.clientY-r.top)/r.height*N),0,N-1);
    k=j*N+i;const c=(N-1)/2;
    x=(i-c)*s.dx;y=(j-c)*s.dx;
  }
  $('fieldInspect').textContent='x '+num(um(x),1)+' µm · y '+num(um(y),1)+' µm · u '+num(s.u[k]*100,3)+' % · ΔT '+num(s.dT[k]*1000,2)+' mK · φ '+num(s.ice[k]?1:s.phi[k],2)+
  (s.eEnabled?' · |E| '+num(s.eMagnitude[k]*1e-6,3)+' MV/m · V '+num(s.ePotential[k],3)+' V':'');
}

// Canvas plots are pixel-perfect and use the same stored physical samples.
function context2D(id) {
  const canvas=$(id),box=canvas.getBoundingClientRect(),dpr=Math.min(window.devicePixelRatio||1,2);
  if(box.width<1)return null;
  const width=Math.round(box.width*dpr),height=Math.round(box.height*dpr);
  if(canvas.width!==width||canvas.height!==height){canvas.width=width;canvas.height=height;}
  const ctx=canvas.getContext('2d');
  ctx.setTransform(dpr,0,0,dpr,0,0);
  return {ctx,w:box.width,h:box.height};
}
function chart(id,series,xTitle,yTitle,opts={}) {
  const p=context2D(id);if(!p)return;
  const {ctx,w,h}=p;
  ctx.clearRect(0,0,w,h);
  const L=47,R=12,T=12,B=29,W=w-L-R,H=h-T-B;
  if(W<10||H<10)return;
  const points=series.flatMap(s=>s.points).filter(q=>Number.isFinite(q[0])&&Number.isFinite(q[1]));
  if(!points.length)return;
  let x0=opts.x0??Math.min(...points.map(p=>p[0])),x1=opts.x1??Math.max(...points.map(p=>p[0]));
  let y0=opts.y0??Math.min(0,...points.map(p=>p[1])),y1=opts.y1??Math.max(...points.map(p=>p[1]));
  if(x1===x0)x1=x0+1;
  if(y1-y0<1e-7)y1=y0+1e-4;
  const paddingY=(y1-y0)*.08;y0-=paddingY;y1+=paddingY;
  const X=x=>L+(x-x0)/(x1-x0)*W,Y=y=>T+H-(y-y0)/(y1-y0)*H;
  ctx.strokeStyle='#274557';ctx.lineWidth=1;ctx.font='9px Inter,system-ui,sans-serif';ctx.textAlign='right';ctx.fillStyle='#7392a5';
  for(let i=0;i<=4;i++){
    const yy=T+H*i/4,val=y1-(y1-y0)*i/4;
    ctx.beginPath();ctx.moveTo(L,yy);ctx.lineTo(w-R,yy);ctx.stroke();
    ctx.fillText((Math.abs(val)<.1?val.toFixed(3):val.toFixed(2)),L-6,yy+3);
  }
  ctx.textAlign='center';
  for(let i=0;i<=4;i++){
    const v=x0+(x1-x0)*i/4,x=X(v);
    ctx.fillText((x1-x0<5?v.toFixed(1):v.toFixed(0)),x,h-13);
  }
  ctx.fillStyle='#789eb1';ctx.font='9px Inter,system-ui,sans-serif';
  ctx.textAlign='left';ctx.fillText(yTitle,L+4,T+8);
  ctx.textAlign='right';ctx.fillText(xTitle,w-R,h-2);
  ctx.save();ctx.beginPath();ctx.rect(L,T,W,H);ctx.clip();
  for(const [i,s] of series.entries()){
    if(!s.points.length)continue;
    ctx.beginPath();
    s.points.forEach(([x,y],k)=>{if(k===0)ctx.moveTo(X(x),Y(y));else ctx.lineTo(X(x),Y(y));});
    ctx.strokeStyle=s.color||palette[i%palette.length];ctx.globalAlpha=.9;
    ctx.lineWidth=s.width||1.7;ctx.lineJoin='round';ctx.stroke();
  }
  ctx.globalAlpha=1;
  if(opts.markerX!==undefined){
    const x=X(opts.markerX);ctx.beginPath();ctx.moveTo(x,T);ctx.lineTo(x,T+H);
    ctx.strokeStyle='#ffffffb0';ctx.setLineDash([3,4]);ctx.stroke();ctx.setLineDash([]);
  }
  ctx.restore();
}
function radialSeries(s){
  const rays=6,radials=[];
  const center=(s.n-1)/2;
  const maxR=(center-1)*s.dx;
  const points=64;
  for(let arm=0;arm<rays;arm++){
    const theta=arm*Math.PI/3,out=[];
    for(let k=0;k<points;k++){
      const r=maxR*k/(points-1),i=clamp(Math.round(center+r*Math.cos(theta)/s.dx),0,s.n-1);
      const j=clamp(Math.round(center+r*Math.sin(theta)/s.dx),0,s.n-1);
      const sampleIndex=s.topology==='triangular'?
        hexCellAtXY({N:s.n,R:center,dx:s.dx,mask:s.mask},
          r*Math.cos(theta),r*Math.sin(theta)):j*s.n+i;
      if(sampleIndex>=0)out.push([um(r),100*s.u[sampleIndex]]);
    }
    radials.push({color:palette[arm],points:out,width:1.55});
  }
  return radials;
}
function drawAllGraphs(s) {
  chart('radialGraph',radialSeries(s),'radius [µm]','u [%]');
  const time=frames.map(f=>f.t);
  chart('growthGraph',[
    {color:palette[0],width:2.3,points:frames.map((f,i)=>[time[i],um(f.morph.equivalentBasalRadius)])},
    {color:palette[2],points:frames.map((f,i)=>[time[i],um(f.morph.halfThickness)])}
  ],'time [s]','µm',{markerX:replay?s.t:undefined});
  chart('modeGraph',[
    {color:palette[0],width:2.3,points:frames.map((f,i)=>[time[i],(f.morph.modes[6].relative||0)*100])},
    {color:palette[1],points:frames.map((f,i)=>[time[i],(f.morph.modes[4].relative||0)*100])}
  ],'time [s]','rel. %',{markerX:replay?s.t:undefined});
}

function download(name,type,data) {
  const url=URL.createObjectURL(new Blob([data],{type}));
  const a=document.createElement('a');
  a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();
  setTimeout(()=>URL.revokeObjectURL(url),2000);
}
function downloadCsv(){
  const columns=['time_s','radius_um','mean_radius_um','half_thickness_um','deposited_mass_pg','latent_energy_J','max_dT_mK','m4_relative','m6_relative','m12_relative','mean_interface_sigma','mean_prism_alpha','alpha_basal','interface_cells','mass_rate_pg_s','basal_velocity_um_s','mean_curvature_per_um','mean_capillary_exponent',
    'max_electric_field_MV_m','flexoelectric_P_uC_m2','mean_dielectric_log_shift','poisson_residual_V','electric_coupled',
    'domain_um','tip_radius_um','halo_cells','thickness_um','halfthickness_over_tipradius','domain_stop','water_saturation_ice_percent','habit',
    'topology','domain_expansions'];
  const rows=[columns.join(',')];
  for(const s of frames){
    const b=s.budget||{},m=s.morph;
    const values=[s.t,um(m.equivalentBasalRadius),um(m.meanRadius),um(m.halfThickness),pg(s.deposited),s.energy,(b.maxDeltaT||0)*1000,
      m.modes[4].relative,m.modes[6].relative,m.modes[12].relative,b.meanSigma,b.meanPrismAlpha,b.alphaBasal,b.interfaceCells,
      b.totalMassRate===undefined?'':pg(b.totalMassRate),b.vBasal===undefined?'':um(b.vBasal),
      b.meanCurvature===undefined?'':b.meanCurvature*1e-6,b.meanCapillaryExponent,
      s.eDiagnostics?s.eDiagnostics.maxFieldVM*1e-6:'',s.eDiagnostics?s.eDiagnostics.maxPolarization*1e6:'',
      b.meanElectricLogShift??0,s.eDiagnostics?.residualRmsV??'',s.eCoupled?1:0,
      s.audit.domainUm,s.audit.tipRadiusUm,s.audit.haloCells,s.audit.thicknessUm,s.audit.plateRatio,
      s.audit.stop?1:0,s.waterSaturationPercent,s.habit,s.topology,s.expansions];
    rows.push(values.map(v=>v??'').join(','));
  }
  download('snowflake-observatory-'+new Date().toISOString().slice(0,10)+'.csv','text/csv;charset=utf-8',rows.join('\r\n')+'\r\n');
}
function downloadPng(){
  if(!renderer)return;
  renderer.render(scene,camera);
  const a=document.createElement('a');
  a.href=renderer.domElement.toDataURL('image/png');a.download='snowflake-observatory-view.png';a.click();
}
function showReplay(index){
  if(!frames.length)return;
  selectedIndex=clamp(index,0,frames.length-1);
  replay=selectedIndex!==frames.length-1;
  if(replay)playing=false;
  $('play').textContent=playing?'❚❚ Pause':'▶ Play';
  $('scrub').value=String(selectedIndex);
  setStatus(replay?'REPLAYING STORED SOLVER STATE':playing?'SOLVER RUNNING':'SOLVER PAUSED');
  sync();
}
function describeSelectedHabit(){
  const habit=GROWTH_HABITS[$('habit').value];
  if(!habit)return;
  $('habitNote').textContent=(habit.supported?'OBSERVED HABIT TARGET · ':'NOT SIMULATED · ')+
    habit.expectation+' '+(habit.limitation||habit.calibration)+'. '+habit.evidence+'.';
  $('preset').disabled=!habit.supported;
  $('preset').title=habit.supported?'Apply exploratory scenario; no guaranteed morphology':habit.limitation;
}
function loadSelectedHabit(){
  const key=$('habit').value,habit=GROWTH_HABITS[key];
  if(!habit?.supported){
    describeSelectedHabit();
    setStatus('UNSUPPORTED HABIT · NEEDLE/COLUMN REQUIRES RESOLVED 3-D',true);
    return;
  }
  const values={
    temp:habit.temperatureC,sigma:habit.sigmaPercent,pressure:habit.pressureKPa,
    anisotropy:habit.anisotropy6,prism:habit.prismSigma0,basalA:habit.basalA,
    seed:habit.seedRadiusUm,gridSize:habit.gridN
  };
  for(const [id,value] of Object.entries(values))$(id).value=String(value);
  $('sdak').checked=habit.sdak;$('capillarity').checked=habit.capillarity;
  // Weather scenarios are neutral; electric forcing is a separate experiment.
  $('electricEnabled').checked=false;$('electricCouple').checked=false;
  $('eField').value='2';$('bendGradient').value='0';
  $('topology').value='triangular';topologyControls();
  activeHabit=key;formatControls();reset(true);describeSelectedHabit();
}
function topologyControls(){
  const hex=$('topology').value==='triangular';
  // The Poisson module is currently discretized on the Cartesian grid.
  // Do not pretend it solves on the triangular mesh.
  for(const name of ['electricEnabled','electricCouple']){
    if(hex)$(name).checked=false;
    $(name).disabled=hex;
  }
  for(const opt of $('field').options){
    if(['electric','potential'].includes(opt.value))opt.disabled=hex;
  }
  if(hex&&['electric','potential'].includes($('field').value))$('field').value='vapor';
  $('topologyValue').textContent=hex?'6-neighbor triangular':'9-point Cartesian';
}
function setupEvents(){
  formatControls();
  describeSelectedHabit();
  topologyControls();
  $('topology').addEventListener('change',()=>{
    activeHabit='custom';topologyControls();setPending(true);
  });
  for(const [id] of labels)$(id).addEventListener('input',()=>{
    activeHabit='custom';formatControls();setPending(true);
  });
  for(const id of ['sdak','capillarity','electricEnabled','electricCouple','epsIce'])$(id).addEventListener('change',()=>{
    activeHabit='custom';formatControls();setPending(true);
  });
  $('apply').addEventListener('click',()=>reset(true));
  $('reset').addEventListener('click',()=>reset(false));
  $('step').addEventListener('click',()=>{
    if(replay){replay=false;selectedIndex=frames.length-1;}
    playing=false;$('play').textContent='▶ Play';advance();
    const audit=current()?.audit;
    setStatus(audit?.stop?'MODEL DOMAIN STOP · '+audit.reason:
      audit?.warning?'MODEL VALIDITY WARNING · '+audit.reason:'SOLVER PAUSED');
  });
  $('play').addEventListener('click',play);
  $('habit').addEventListener('change',describeSelectedHabit);
  $('preset').addEventListener('click',loadSelectedHabit);
  $('field').addEventListener('change',()=>{if(current())drawField(current());});
  $('fieldCanvas').addEventListener('pointermove',onFieldMove);
  $('fieldCanvas').addEventListener('pointerleave',()=>$('fieldInspect').textContent='HOVER FIELD TO INSPECT A CELL');
  $('rotate').addEventListener('change',()=>{if(controls)controls.autoRotate=$('rotate').checked;});
  $('viewTop').addEventListener('click',viewTop);
  $('viewIso').addEventListener('click',viewIso);
  $('csv').addEventListener('click',downloadCsv);
  $('png').addEventListener('click',downloadPng);
  $('speed').addEventListener('change',()=>accumulator=0);
  $('scrub').addEventListener('input',()=>showReplay(Number($('scrub').value)));
  $('live').addEventListener('click',()=>{
    replay=false;selectedIndex=frames.length-1;$('scrub').value=String(selectedIndex);
    setStatus(playing?'SOLVER RUNNING':'SOLVER PAUSED');sync();
  });
  window.addEventListener('keydown',event=>{
    if(event.target instanceof HTMLInputElement || event.target instanceof HTMLSelectElement)return;
    if(event.code==='Space'){event.preventDefault();play();}
    if(event.key.toLowerCase()==='s'){event.preventDefault();$('step').click();}
    if(event.key.toLowerCase()==='r'){event.preventDefault();reset(false);}
  });
  new ResizeObserver(()=>{if(current())drawAllGraphs(current());}).observe(document.querySelector('.graphs'));
}

function loop(now){
  requestAnimationFrame(loop);
  const elapsed=Math.min(.25,Math.max(0,(now-lastFrame)/1000));lastFrame=now;
  if(playing&&!replay&&state){
    accumulator+=elapsed*BASE_STEPS_PER_SECOND*Number($('speed').value);
    let count=0;
    while(accumulator>=1&&count<2&&playing){
      accumulator--;advance();count++;
    }
    // Prevent a stalled tab from demanding an unbounded backlog of solver steps.
    accumulator=Math.min(accumulator,2.5);
  }
  if(controls){controls.update();renderer.render(scene,camera);}
}
setupEvents();
try{
  initializeScene();
  reset(true);
}catch(e){fail(e);}
requestAnimationFrame(loop);
