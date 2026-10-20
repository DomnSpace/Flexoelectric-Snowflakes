/**
 * Snowflake Observatory experiment definitions.
 *
 * A habit name is an OBSERVATIONAL TARGET, never a solver output guarantee.
 * Temperature, water/ice saturation and pressure have observational units;
 * the facet coefficients below are explicit UNCALIBRATED DEMONSTRATOR values.
 * A regime timeline changes the far-field environment during one run. It
 * never resets the seed or paints a target shape into the phase field.
 */
import {GROWTH_HABITS} from './habits25d.js';
import {saturationVaporPressureWaterMK,saturationVaporPressureIceMK} from './thermo.js';
import {saturationMassDensityIce} from './postseed.js';
import {iceIhProperties} from './ice-ih-eos.js';

export const EXPERIMENT_SCHEMA='snowflake-observatory/experiment';
export const EXPERIMENT_VERSION=1;
export const RAMP_SECONDS=3;

export function waterSaturationIcePercent(tC){
  const T=273.15+tC;
  if(!(T>110))throw new RangeError('Invalid cloud temperature');
  return 100*(saturationVaporPressureWaterMK(T)/saturationVaporPressureIceMK(T)-1);
}
export function cloudCondition(tC,sigmaPercent){
  if(!Number.isFinite(tC)||!Number.isFinite(sigmaPercent)||
    sigmaPercent<0)throw new RangeError('Invalid far-field cloud condition');
  const limit=waterSaturationIcePercent(tC);
  const overWaterPercent=100*(1+sigmaPercent/100)/(1+limit/100);
  return {
    temperatureC:tC,iceSupersaturationPercent:sigmaPercent,
    waterSaturationIcePercent:limit,
    relativeHumidityOverWaterPercent:overWaterPercent,
    aboveWaterSaturation:sigmaPercent>limit+1e-9,
    classification:sigmaPercent>limit+1e-9?
      'above-water-saturation: laboratory or transient only':
      sigmaPercent===0?'ice-saturated reference':'below-or-at-water-saturation corridor',
    // Not a morphology threshold or strict exclusion for a cloud parcel.
    source:'Murphy & Koop 2005 vapor saturation pressure over ice and water'
  };
}

const baseline=Object.freeze({
  temperatureC:-15,sigmaPercent:14,pressureKPa:101,
  anisotropy6:.50,kineticPhaseDeg:30,surfaceEnergyEpsilon:0,
  prismSigma0:.0035,basalSigma0:.01,basalA:.004,
  seedRadiusUm:1.5,seedShape:'circle',seedPhaseDeg:0,
  gridN:129,localHeight:false,sdak:true,capillarity:true,
  topology:'triangular',autoExpand:true,
  electricEnabled:false,electricCouple:false,eFieldMVm:2,
  eAngleDeg:0,epsIce:3.2,bendGradient:0,bendAngleDeg:0,muFlexo:1.65
});
const source='Nakaya/Libbrecht temperature–habit observations; kinetic parameters illustrative / uncalibrated';
const make=(id,{name,category='plate',supported=true,description,
  controls={},stages=null,limit=null})=>{
  const all=Object.freeze({...baseline,...controls});
  const timeline=stages?Object.freeze(stages.map(s=>Object.freeze({...s}))):null;
  return Object.freeze({id,name,category,supported,description,
    controls:all,stages:timeline,limitations:limit||null,
    provenance:source,kind:stages?'growth-history':'constant-environment',
    calibration:'experiment: uncalibrated attachment and SDAK values'});
};
export const EXPERIMENT_PRESETS=Object.freeze({
  plate:make('plate',{
    name:'Compact hexagonal plate · candidate',
    description:'Around −15 °C, low ice supersaturation typically favors simple plates.',
    controls:{...GROWTH_HABITS.plate,kineticPhaseDeg:0,surfaceEnergyEpsilon:0,
      seedShape:'hex',gridN:97,sdak:false}
  }),
  sector:make('sector',{
    name:'Sectored plate · candidate',
    description:'Intermediate ice supersaturation near −15 °C can favor sectored plates.',
    controls:{...GROWTH_HABITS.sector,kineticPhaseDeg:0,surfaceEnergyEpsilon:0,sdak:false}
  }),
  dendrite:make('dendrite',{
    name:'Stellar dendrite · candidate',
    description:'At −15 °C and high ice supersaturation, six-arm dendritic growth is observed. Not guaranteed by this solver.',
    controls:{...GROWTH_HABITS.dendrite}
  }),
  history:make('history',{
    name:'Plate → arms → facets → branches',
    description:'One continuous ice crystal sees four consecutive humidity stages at −15 °C; shapes are not forced.',
    controls:{sigmaPercent:3,gridN:129,localHeight:true,sdak:true,seedShape:'hex'},
    stages:[
      {t:0,label:'Compact core',temperatureC:-15,sigmaPercent:3},
      {t:15,label:'Arm growth',temperatureC:-15,sigmaPercent:13},
      {t:45,label:'Facet episode',temperatureC:-15,sigmaPercent:5},
      {t:65,label:'Renewed branching',temperatureC:-15,sigmaPercent:14}
    ]
  }),
  warmPlate:make('warmPlate',{
    name:'Warm plate · limited model experiment',
    description:'Near −2 °C plates can occur, but facet kinetic temperature dependence is not calibrated here.',
    controls:{temperatureC:-2,sigmaPercent:1.1,sdak:false,
      anisotropy6:.22,basalA:.009,prismSigma0:.012,gridN:97}
  }),
  needle:make('needle',{
    name:'−5 °C needle / column · observation only',category:'column',supported:false,
    description:'Near −5 °C columns/needles are observed; 2.5-D basal mesh cannot resolve their longitudinal geometry.',
    controls:{temperatureC:-5,sigmaPercent:4,sdak:false},
    limit:'Full 3-D tip, vertical vapor and facet-resolved column solver required'
  }),
  cappedColumn:make('cappedColumn',{
    name:'Capped column · observation only',category:'column',supported:false,
    description:'Vertical habit changes and cap plates require growth of a 3-D column, not a 2.5-D plate.',
    controls:{temperatureC:-5,sigmaPercent:4},
    limit:'Cannot create true capped column without 3-D growth history'
  })
});

export function validatePreset(preset){
  if(!preset||!preset.controls||!preset.id)throw new TypeError('Missing versioned experiment preset');
  const c=preset.controls;
  if(!(c.temperatureC>=-35&&c.temperatureC<=-1&&c.sigmaPercent>=0&&
    c.sigmaPercent<=35&&c.pressureKPa>=60&&c.pressureKPa<=110))
    throw new RangeError('Preset atmospheric conditions outside configured controls');
  if(preset.stages){
    if(c.temperatureC!==preset.stages[0].temperatureC||c.sigmaPercent!==preset.stages[0].sigmaPercent)
      throw new Error('Initial timeline stage differs from preset atmosphere');
    for(let i=0;i<preset.stages.length;i++){
      const s=preset.stages[i];
      if(s.t<0||(i&&s.t<=preset.stages[i-1].t))throw new Error('Timeline stages out of order');
      if(i&&s.t-preset.stages[i-1].t<RAMP_SECONDS)throw new Error('Overlapping transition ramps');
      if(s.temperatureC<-35||s.temperatureC>-1||s.sigmaPercent<0||s.sigmaPercent>35)
        throw new Error('Invalid stage atmosphere');
    }
  }
  return true;
}

/** Single true stage at any time; linear 3-s ramp into each changed regime. */
export function sampleGrowthHistory(stages,time){
  if(!Array.isArray(stages)||!stages.length||time<0||!Number.isFinite(time))
    throw new RangeError('Invalid growth-history query');
  let index=0;
  for(let i=1;i<stages.length;i++)if(time>=stages[i].t)index=i;
  const s=stages[index];
  let temperatureC=s.temperatureC,sigmaPercent=s.sigmaPercent;
  let fraction=1;
  if(index>0&&time<s.t+RAMP_SECONDS){
    fraction=Math.max(0,Math.min(1,(time-s.t)/RAMP_SECONDS));
    const prev=stages[index-1];
    temperatureC=prev.temperatureC+(s.temperatureC-prev.temperatureC)*fraction;
    sigmaPercent=prev.sigmaPercent+(s.sigmaPercent-prev.sigmaPercent)*fraction;
  }
  return {index,label:s.label,temperatureC,sigmaPercent,rampFraction:fraction,
    ...cloudCondition(temperatureC,sigmaPercent)};
}
export function makeExperimentRun(presetId='dendrite'){
  const preset=EXPERIMENT_PRESETS[presetId];
  if(!preset?.supported)throw new Error('Unsupported 2.5-D habit '+presetId);
  validatePreset(preset);
  const stages=preset.stages?
    preset.stages.map(s=>({...s})):
    [{t:0,label:preset.name,temperatureC:preset.controls.temperatureC,
      sigmaPercent:preset.controls.sigmaPercent}];
  return {schema:EXPERIMENT_SCHEMA,version:EXPERIMENT_VERSION,
    id:preset.id,stages,stageIndex:0,
    lastTemperatureC:preset.controls.temperatureC,
    lastSigmaPercent:preset.controls.sigmaPercent,
    transitions:[],source:preset.provenance,
    calibrated:false,verticalVaporTransportSolved:false};
}
/**
 * Rebase c_v into the new physical reference concentration, rather than
 * reinterpreting the old dimensionless supersaturation at a new temperature.
 *
 * On an open Dirichlet domain, the outer shell is overwritten with the NEW
 * prescribed vapor concentration. Such flux at the boundary is an external
 * source/sink and must NOT be presented as a closed gas mass account.
 */
export function makeCustomExperimentRun(temperatureC,sigmaPercent){
  cloudCondition(temperatureC,sigmaPercent);
  return {schema:EXPERIMENT_SCHEMA,version:EXPERIMENT_VERSION,
    id:'custom',stages:[{t:0,label:'Custom constant conditions',
      temperatureC,sigmaPercent}],stageIndex:0,
    lastTemperatureC:temperatureC,lastSigmaPercent:sigmaPercent,
    transitions:[],source:'User-defined atmosphere; uncalibrated kinetics',
    calibrated:false,verticalVaporTransportSolved:false};
}
export function applyGrowthEnvironment(state,run,time){
  if(!run||run.schema!==EXPERIMENT_SCHEMA||run.version!==1)
    throw new Error('No versioned environment history');
  if(state.topology!=='triangular')throw new Error('Growth history requires triangular solver');
  const s=sampleGrowthHistory(run.stages,time);
  if(run.stageIndex!==s.index){
    run.transitions.push({time,stageIndex:s.index,label:s.label,
      temperatureC:s.temperatureC,sigmaPercent:s.sigmaPercent});
    run.stageIndex=s.index;
  }
  if(Math.abs(s.temperatureC-run.lastTemperatureC)<1e-10&&
    Math.abs(s.sigmaPercent-run.lastSigmaPercent)<1e-10){
    run.stageIndex=s.index;return {changed:false,...s};
  }
  const oldT=state.Tinf,oldC=saturationMassDensityIce(oldT),
    nextT=s.temperatureC+273.15,newC=saturationMassDensityIce(nextT);
  // Ice Ih density varies with temperature: preserve the existing ice mass
  // when switching thermodynamic reference, via a uniform volume rescaling.
  // This is a density-compatibility closure, not a solved heat expansion PDE.
  if(state.localHeight&&state.topHeight&&Math.abs(nextT-oldT)>1e-11){
    const oldRho=iceIhProperties(oldT,Math.min(state.pressurePa,210e6)).rho;
    const newRho=iceIhProperties(nextT,Math.min(state.pressurePa,210e6)).rho;
    const volumeScale=oldRho/newRho;
    for(let k=0;k<state.topHeight.length;k++){
      state.topHeight[k]*=volumeScale;state.bottomHeight[k]*=volumeScale;
    }
    state.basalHalfThickness*=volumeScale;
  }
  state.Tinf=nextT;
  state.sigmaInf=s.sigmaPercent/100;
  for(let k=0;k<state.sigma.length;k++){
    if(!state.mask[k])continue;
    state.sigma[k]=Math.max(-.9999999,(oldC*(1+state.sigma[k]))/newC-1);
  }
  for(const k of state.edges)state.sigma[k]=state.sigmaInf;
  // Internal dT is an *absolute temperature anomaly* and is unchanged.
  // Its subsequent heat relaxation uses the newly imposed Tinf.
  run.lastTemperatureC=s.temperatureC;
  run.lastSigmaPercent=s.sigmaPercent;
  return {changed:true,...s};
}
export function cloneExperimentState(run){
  if(!run)return null;
  return {schema:run.schema,version:run.version,id:run.id,
    stageIndex:run.stageIndex,stageLabel:run.stages[run.stageIndex]?.label||'',
    stages:run.stages.map(s=>({...s})),
    lastTemperatureC:run.lastTemperatureC,
    lastSigmaPercent:run.lastSigmaPercent,
    transitions:run.transitions.map(t=>({...t})),
    calibrated:false};
}
for(const p of Object.values(EXPERIMENT_PRESETS))validatePreset(p);
