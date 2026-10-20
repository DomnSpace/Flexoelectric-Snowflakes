const tests = [
  './thermo.test.mjs',
  './activation.test.mjs',
  './trajectory.test.mjs',
  './supercooled-eos.test.mjs',
  './composition-feedback.test.mjs',
  './supercool.test.mjs',
  './ice-ih.test.mjs',
  './embryo.test.mjs',
  './nucleation.test.mjs',
  './postseed.test.mjs',
  './growth25d.test.mjs',
  './habit-domain.test.mjs',
  './experiment-regimes.test.mjs',
  './hex25d.test.mjs',
  './anisotropy-phase.test.mjs',
  './pass7-snapshot-layers.test.mjs',
  './pass8-thickness-budget.test.mjs',
  './pass9-actual-height.test.mjs',
  './five-layer-baseline.test.mjs',
  './layered-gas-fv.test.mjs',
  './layered-implicit.test.mjs',
  './five-layer-ice-front.test.mjs',
  './layered-thermal.test.mjs',
  './electric25d.test.mjs',
  './fieldlines.test.mjs',
];

for (const test of tests) {
  console.log(`\n=== ${test} ===`);
  await import(test);
}
console.log(`\nALL ${tests.length} PHYSICS REGRESSION FILES COMPLETED`);
