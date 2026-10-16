import {
  createFiveLayerCrystal,stepFiveLayerCrystal,
  fiveLayerOccupancy,fiveLayerDiagnostics
} from './solver/layered-ice-front.js';
import {layeredImplicitTrackedWaterKg} from './solver/layered-vapor-implicit.js';
import {layeredThermalEnergyJ} from './solver/layered-thermal.js';

const $=id=>document.getElementById(id);
const nformat=x=>Number.isFinite(x)?x.toExponential(3):'—';
let state=null,initialWater=0,last=null,busy=false;
const N=()=>Number($('grid').value),dt=()=>Number($('dt').value);
function reset(){
  if(busy)return;
  state=createFiveLayerCrystal({N:N(),seedRadiusCells:0,
    prismAttachmentMPerS:5e-3,basalAttachmentMPerS:6e-3,
    Tinf:258.15,sigmaInf:.14,pressurePa:101325,
    thermalCoupling:$('thermalMode').value==='coupled'});
  initialWater=layeredImplicitTrackedWaterKg(state);
  last=null;
  $('status').style.color='';
  $('status').textContent='READY · INDEPENDENT ICE PLANES · '+
    (state.thermal?'THERMAL COUPLED':'NO HEAT CONTROL');
  draw();
}
function hexPath(ctx,x,y,r){
  ctx.beginPath();
  for(let k=0;k<6;k++){
    const a=Math.PI/6+k*Math.PI/3;
    const u=x+r*Math.cos(a),v=y+r*Math.sin(a);
    if(k===0)ctx.moveTo(u,v);else ctx.lineTo(u,v);
  }
  ctx.closePath();
}
function drawLayer(canvas,phi,seed,temperature,maxTemperature){
  const ctx=canvas.getContext('2d'),w=canvas.width,h=canvas.height;
  ctx.clearRect(0,0,w,h);
  ctx.fillStyle='#061827';ctx.fillRect(0,0,w,h);
  const N=state.N,R=state.R;
  const unit=Math.min((w-36)/(2*R+2),(h-36)/(Math.sqrt(3)*R+2));
  const hexR=unit/Math.sqrt(3)+.28;
  const gain=Number($('gain').value);
  const heatMode=$('mapMode').value==='temperature';
  let hasNew=0;
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    const k=j*N+i;if(!state.hexMask[k])continue;
    const q=i-R,r=j-R,x=w/2+unit*(q+r/2);
    const y=h/2-unit*Math.sqrt(3)/2*r;
    const ice=phi[k],isSeed=seed[k];
    const visible=isSeed?1:Math.max(0,Math.min(1,ice*gain));
    const shade=Math.sqrt(visible);
    if(ice>0&&!isSeed)hasNew++;
    if(heatMode){
      const warmth=Math.min(1,Math.max(0,(temperature?.[k]??0)/Math.max(maxTemperature,1e-12)));
      ctx.fillStyle='rgba('+Math.round(35+205*warmth)+','+
        Math.round(80+107*(1-warmth))+','+
        Math.round(160-100*warmth)+','+(.14+.86*Math.sqrt(warmth))+')';
    }else{
      ctx.fillStyle=isSeed?'rgb(230,247,255)':
        ice>0?'rgba(135,219,249,'+(0.09+0.90*shade)+')':'rgba(35,82,105,.27)';
    }
    hexPath(ctx,x,y,hexR);ctx.fill();
    if(!heatMode&&ice>0&&!isSeed&&visible>.5){
      ctx.strokeStyle='rgba(213,250,255,.35)';ctx.lineWidth=.5;ctx.stroke();
    }
  }
  return hasNew;
}
function draw(){
  if(!state)return;
  const maps=fiveLayerOccupancy(state),diag=fiveLayerDiagnostics(state);
  const k=state.gasHaloSlabs;
  const hottest=state.temperatureDeltaK?Math.max(...state.temperatureDeltaK):0;
  for(let l=0;l<5;l++){
    const seed=state.seedSolidMask.subarray((k+l)*state.nxy,(k+l+1)*state.nxy);
    const temp=state.temperatureDeltaK?.subarray((k+l)*state.nxy,(k+l+1)*state.nxy);
    const active=drawLayer($('layer'+l),maps[l],seed,temp,hottest);
    const newlyDeposited=state.iceDepositionPerSlabKg[k+l];
    $('mass'+l).textContent=nformat(diag.massKg[l])+
      ' kg total · +'+nformat(newlyDeposited)+' kg deposited · '+
      active+' newly occupied sites'+
      (temp?' · peak ΔT '+Math.max(...temp).toFixed(5)+' K':'');
  }
  const gains=Array.from(state.iceDepositionPerSlabKg.slice(k,k+5));
  const max=Math.max(...gains,1e-99);
  $('bars').replaceChildren(...gains.map((mass,i)=>{
    const row=document.createElement('div');row.className='bar-row';
    const lab=document.createElement('span');lab.textContent='Layer '+(i+1);
    const track=document.createElement('div');track.className='track';
    const fill=document.createElement('div');fill.className='fill';
    fill.style.width=(mass/max*100).toFixed(3)+'%';track.append(fill);
    const amount=document.createElement('span');amount.textContent=nformat(mass)+' kg';
    row.append(lab,track,amount);return row;
  }));
  $('time').textContent=state.timeSeconds.toFixed(3)+' s';
  $('ice').textContent=nformat(state.totalDepositedKg)+' kg';
  $('residual').textContent=last?nformat(last.relativeResidual):'—';
  $('balance').textContent=nformat(
    layeredImplicitTrackedWaterKg(state)-initialWater-
    (state.externalMassInKg-state.externalMassOutKg))+' kg';
  $('maxTemp').textContent=state.thermal?hottest.toFixed(6)+' K':'OFF';
  $('heatBalance').textContent=state.thermal?
    nformat(state.thermalLatentInjectedJ-state.thermalBoundaryHeatOutJ-
      layeredThermalEnergyJ(state))+' J':'OFF';
  if(last){
    $('status').textContent='COMPUTED · '+last.iterations+
      ' CG ITERATIONS · 5 DISTINCT PHYSICAL FIELDS'+
      (last.thermal?' · LATENT HEAT COUPLED':' · NO-HEAT CONTROL');
  }
}
async function advance(steps){
  if(busy||!state)return;
  busy=true;
  for(const b of ['step','step10','reset'])$(b).disabled=true;
  try{
    for(let i=0;i<steps;i++){
      last=stepFiveLayerCrystal(state,dt(),{tolerance:1e-8,maxIterations:1600});
      draw();
      // Yield to browser painting and input between steps.
      if(i<steps-1)await new Promise(resolve=>setTimeout(resolve,0));
    }
  }catch(e){
    $('status').textContent='NUMERICAL STOP · '+e.message;
    $('status').style.color='#ffa69b';
    console.error(e);
  }finally{
    busy=false;
    for(const b of ['step','step10','reset'])$(b).disabled=false;
  }
}
$('step').addEventListener('click',()=>advance(1));
$('step10').addEventListener('click',()=>advance(10));
$('reset').addEventListener('click',reset);
$('grid').addEventListener('change',reset);
$('dt').addEventListener('change',()=>{});
$('gain').addEventListener('change',draw);
$('mapMode').addEventListener('change',draw);
$('thermalMode').addEventListener('change',reset);
reset();
