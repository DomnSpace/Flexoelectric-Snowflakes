/**
 * Nakaya/Libbrecht-informed *scenarios*, not a calibrated ice-habit model.
 * Temperature/humidity tendencies are observational. Their 2.5-D kinetic
 * coefficients are deliberately exposed as uncalibrated exploration values.
 *
 * A basal-plane phase fraction plus one shared thickness cannot reproduce
 * needle/column growth at -5 °C. Never fake it with a tall extrusion.
 */
export const GROWTH_HABITS=Object.freeze({
  plate:Object.freeze({
    name:'Simple hexagonal plate',category:'plate',supported:true,
    temperatureC:-15,sigmaPercent:3,pressureKPa:101,
    anisotropy6:.22,prismSigma0:.012,basalSigma0:.01,basalA:.009,
    seedRadiusUm:1.5,gridN:97,sdak:false,capillarity:true,
    expectation:'Near −15 °C, low supersaturation favors compact plates. A hexagonal outline is not guaranteed.',
    evidence:'Libbrecht 2005, Rep. Prog. Phys. 68 855–895; Nakaya habit diagram',
    calibration:'illustrative facet kinetics; not measured at this temperature'
  }),
  sector:Object.freeze({
    name:'Sectored plate experiment',category:'plate',supported:true,
    temperatureC:-15,sigmaPercent:9,pressureKPa:101,
    anisotropy6:.32,prismSigma0:.007,basalSigma0:.01,basalA:.006,
    seedRadiusUm:1.5,gridN:129,sdak:false,capillarity:true,
    expectation:'Intermediate supersaturation near −15 °C can favor more structured and sectored plates.',
    evidence:'Libbrecht 2005; Nakaya habit diagram',
    calibration:'illustrative kinetic values; sector boundaries not imposed'
  }),
  dendrite:Object.freeze({
    name:'Stellar dendrite candidate',category:'plate',supported:true,
    temperatureC:-15,sigmaPercent:18,pressureKPa:101,
    anisotropy6:.50,prismSigma0:.0035,basalSigma0:.01,basalA:.004,
    seedRadiusUm:1.5,gridN:129,sdak:true,capillarity:true,
    expectation:'Near −15 °C and high supersaturation, stellar dendrites occur in nature. This solver must demonstrate branching rather than assume it.',
    evidence:'Libbrecht 2005, Fig. 2; high supersaturation gives greater morphological complexity',
    calibration:'SDAK and the kinetic values are experimental / not calibrated'
  }),
  needle:Object.freeze({
    name:'Needle / column (reference only)',category:'column',supported:false,
    temperatureC:-5,sigmaPercent:18,pressureKPa:101,
    expectation:'Near −5 °C, high-humidity crystals commonly form needles or columns.',
    evidence:'Libbrecht 2005; Nakaya habit diagram',
    limitation:'This basal-plane 2.5-D solver has no vertical vapor/thermal PDE, column tip, or facet-resolved 3-D growth. No needle preset is simulated.'
  })
});
export const GRID_CHOICES=Object.freeze([73,97,129,161]);
export function modelDomainUm(N,dx){return (N-1)*dx*1e6;}

/**
 * Score the nearest genuinely occupied point to the square Dirichlet boundary.
 * Equivalent-area radius is *not* safe: six long arms may touch the boundary
 * while the same total area has a tiny equivalent radius.
 * Also report plate aspect validity of the shared-height 2.5-D closure.
 */
export function auditGrowthDomain(state,{warningCells=12,stopCells=6,plateRatioWarning=.50}={}){
  const N=state.N,dx=state.dx,center=(N-1)/2;
  let minHalo=(N-1)/2,maxRadius=0,occupied=0;
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    const k=j*N+i;
    if(!state.ice[k]&&state.phi[k]<.15)continue;
    occupied++;
    minHalo=Math.min(minHalo,i,j,N-1-i,N-1-j);
    maxRadius=Math.max(maxRadius,Math.hypot(i-center,j-center)*dx);
  }
  const thickness=2*state.basalHalfThickness;
  const plateRatio=state.basalHalfThickness/Math.max(dx,maxRadius);
  const haloUm=minHalo*dx*1e6;
  const domainUm=modelDomainUm(N,dx);
  const warn=minHalo<=warningCells;
  const stop=minHalo<=stopCells;
  // This is a warning, NOT a physical crystal exclusion condition:
  // a thick crystal may exist but violates the uniform-thickness plate model.
  const plateWarning=plateRatio>plateRatioWarning;
  return {
    N,dx,domainUm,occupied,haloCells:minHalo,haloUm,
    tipRadiusUm:maxRadius*1e6,thicknessUm:thickness*1e6,
    plateRatio,warning:warn,stop,plateWarning,
    reason:stop?'RESOLVED CRYSTAL TOO CLOSE TO SQUARE DIRICHLET BOUNDARY':
      warn?'SQUARE DIRICHLET BOUNDARY MAY BIAS MORPHOLOGY':
      plateWarning?'SHARED-THICKNESS PLATE APPROXIMATION STRETCHED':
      'WITHIN DEMONSTRATOR VALIDITY BOUNDS'
  };
}
