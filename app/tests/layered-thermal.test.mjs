import assert from 'node:assert/strict';
import {createFiveLayerCrystal,stepFiveLayerCrystal} from '../solver/layered-ice-front.js';
import {layeredGasIndex} from '../solver/layered-gas-fv.js';
import {layeredImplicitTrackedWaterKg} from '../solver/layered-vapor-implicit.js';
import {
  attachLayeredThermal,stepLayeredThermal,layeredThermalEnergyJ,
  layeredThermalCapacityAndK
} from '../solver/layered-thermal.js';

function close(a,b,relative=2e-7,absolute=1e-22){
  assert.ok(Number.isFinite(a)&&Number.isFinite(b)&&
    Math.abs(a-b)<=Math.max(absolute,relative*Math.max(Math.abs(a),Math.abs(b))),
    'expected '+a+' ~ '+b);
}
const input={N:11,dx:.75e-6,dz:.75e-6,seedRadiusCells:0,
  Tinf:258.15,sigmaInf:.14,pressurePa:101325,
  prismAttachmentMPerS:5e-3,basalAttachmentMPerS:6e-3};
const closed=createFiveLayerCrystal({...input,
  boundaries:{lateral:'sealed',top:'sealed',bottom:'sealed'}});
attachLayeredThermal(closed);
assert.equal(closed.heatTransportSolved,true);
const cap=layeredThermalCapacityAndK(closed);
const core=layeredGasIndex(closed,0,0,4);
assert.ok(cap.capacity[core]>0);
assert.ok(cap.conductivity[core]>cap.conductivity[layeredGasIndex(closed,3,0,4)]);
const beforeWater=layeredImplicitTrackedWaterKg(closed);
const dm=1e-20;
const warm=stepLayeredThermal(closed,.02,[{toIceIndex:core,massKg:dm}]);
assert.ok(warm.maxDeltaK>0,'latent heat must cause actual local temperature rise');
assert.ok(warm.relativeResidual<2e-9);
close(layeredThermalEnergyJ(closed),dm*closed.Ls,1e-7,1e-22);
close(closed.thermalLatentInjectedJ,dm*closed.Ls,1e-11,1e-24);
assert.ok(Math.abs(closed.thermalBoundaryHeatOutJ)<1e-26,
  'insulated thermal domain cannot leak heat');
assert.ok(closed.localIceEqKgM3[core]>closed.ceq,
  'warmer ice surface must require a higher equilibrium vapor concentration');
close(layeredImplicitTrackedWaterKg(closed),beforeWater,1e-10,1e-28);
const smooth=stepLayeredThermal(closed,.015,[]);
assert.ok(smooth.relativeResidual<2e-9);
close(layeredThermalEnergyJ(closed),dm*closed.Ls,1e-7,1e-22);

const symmetry=createFiveLayerCrystal({...input,
  boundaries:{lateral:'sealed',top:'sealed',bottom:'sealed'},thermalCoupling:true});
const cm=layeredGasIndex(symmetry,0,0,4),delta=2e-20;
stepLayeredThermal(symmetry,.01,[{toIceIndex:cm,massKg:delta}]);
close(symmetry.temperatureDeltaK[layeredGasIndex(symmetry,0,0,3)],
  symmetry.temperatureDeltaK[layeredGasIndex(symmetry,0,0,5)],3e-6,1e-9);
close(symmetry.temperatureDeltaK[layeredGasIndex(symmetry,1,0,4)],
  symmetry.temperatureDeltaK[layeredGasIndex(symmetry,-1,0,4)],3e-6,1e-9);

// Dirichlet ambient Tinf boundaries: heat is transported out and explicitly
// accounted, rather than vanished by an arbitrary cooling multiplier.
const open=createFiveLayerCrystal({...input,thermalCoupling:true});
const oh=stepLayeredThermal(open,.025,[{toIceIndex:core,massKg:dm}]);
assert.ok(oh.heatBoundaryOutJ>0);
assert.ok(oh.energyStoredJ>=0);
close(oh.energyStoredJ+oh.heatBoundaryOutJ,oh.depositedLatentJ,1e-7,1e-22);

// Coupled growth: vapor-phase conservation still holds, heat injection equals
// actual ice mass * Ls, local ice vapor equilibrium is altered on next step.
// Unlike a palette effect, this is a causal change to the actual vapor sink.
const coupled=createFiveLayerCrystal({...input,thermalCoupling:true});
const uncoupled=createFiveLayerCrystal({...input,thermalCoupling:false});
const start=layeredImplicitTrackedWaterKg(coupled);
for(let i=0;i<4;i++){
  const withHeat=stepFiveLayerCrystal(coupled,.015);
  stepFiveLayerCrystal(uncoupled,.015);
  assert.ok(withHeat.thermal);
  assert.equal(withHeat.heatTransportSolved,true);
  assert.ok(withHeat.thermal.relativeResidual<2e-9);
  assert.ok(withHeat.thermal.energyBalanceResidualJ<1e-18);
  assert.ok(withHeat.postGrowthBudgetResidualKg<1e-20);
}
assert.ok(coupled.thermalLatentInjectedJ>0);
assert.ok(coupled.thermalBoundaryHeatOutJ>0);
assert.ok(coupled.temperatureDeltaK.some(t=>t>0));
assert.ok(coupled.localIceEqKgM3.some(v=>v>coupled.ceq));
close(coupled.thermalLatentInjectedJ,coupled.latentReleasedJ,1e-7,1e-23);
close(layeredImplicitTrackedWaterKg(coupled)-start,
  coupled.externalMassInKg-coupled.externalMassOutKg,2e-7,1e-25);
assert.equal(uncoupled.heatTransportSolved,false);
assert.ok(coupled.totalDepositedKg<=uncoupled.totalDepositedKg*(1+2e-5),
  'latent warming should not spuriously accelerate otherwise identical uptake');

// A failed (unconverged) heat solve must not alter its thermal state.
const fail=createFiveLayerCrystal({...input,thermalCoupling:true});
const oldEnergy=layeredThermalEnergyJ(fail),oldTemp=fail.temperatureDeltaK.slice();
assert.throws(()=>stepLayeredThermal(fail,1,[{toIceIndex:core,massKg:dm}],{
  maxIterations:1,tolerance:1e-15
}),/did not converge/);
assert.deepEqual(fail.temperatureDeltaK,oldTemp);
assert.equal(layeredThermalEnergyJ(fail),oldEnergy);
console.log('layered-thermal.test.mjs PASS',{
  latentJ:warm.depositedLatentJ,
  closedHeatResidualJ:warm.energyBalanceResidualJ,
  openHeatExportJ:oh.heatBoundaryOutJ,
  coupledHeatExportJ:coupled.thermalBoundaryHeatOutJ,
  maxDeltaK:warm.maxDeltaK,
  coupledIceKg:coupled.totalDepositedKg,
  uncoupledIceKg:uncoupled.totalDepositedKg
});
