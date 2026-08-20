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
];

for (const test of tests) {
  console.log(`\n=== ${test} ===`);
  await import(test);
}
console.log(`\nALL ${tests.length} PHYSICS REGRESSION FILES COMPLETED`);
