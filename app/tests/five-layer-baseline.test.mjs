/**
 * Checkpoint: legacy fields MUST remain unchanged while constructing the new
 * opt-in multi-z-layer model. Deliberately compares against externally pinned
 * numerical results from commit 77f118c..., not a freshly regenerated
 * fixture from the changed solver under test.
 */
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {captureLegacyBaseline} from './legacy-baseline-capture.mjs';

const reference=JSON.parse(readFileSync(
  new URL('./fixtures/five-layer-legacy-v1.json',import.meta.url),'utf8'));
const observed=captureLegacyBaseline();
assert.equal(reference.schema,'snowflake/legacy-baseline-v1');
assert.equal(observed.legacyCommit,'77f118c091a2fd1ae278c899e7d5f6a1715ebc54');
assert.equal(reference.legacyCommit,observed.legacyCommit);
assert.equal(observed.cases.length,3);
assert.deepEqual(observed.cases.map(c=>c.case),['hex0','hex30','cartesian']);

let checked=0;
for(let i=0;i<observed.cases.length;i++){
  const got=observed.cases[i],want=reference.cases[i];
  assert.equal(got.case,want.case);
  assert.equal(got.topology,want.topology);
  assert.equal(got.gridN,want.gridN);
  assert.equal(got.gridSpacingM,want.gridSpacingM);
  assert.equal(got.dtSeconds,want.dtSeconds);
  assert.equal(got.steps,want.steps);
  for(let t=0;t<got.snapshots.length;t++){
    assert.equal(got.snapshots[t].step,want.snapshots[t].step);
    assert.equal(got.snapshots[t].fieldsSHA256,want.snapshots[t].fieldsSHA256,
      'Legacy '+got.case+' grid changed at physical step '+want.snapshots[t].step+
      ' (ice/phase/vapor/temperature/source SHA256 mismatch)');
    checked++;
  }
  for(let t=0;t<got.measurements.length;t++){
    const a=got.measurements[t],b=want.measurements[t];
    for(const name of ['t','depositedMassKg','totalLatentHeatJ','halfThicknessM',
      'm6','interfaceCells']){
      const av=a[name],bv=b[name];
      if(av===null||bv===null)assert.equal(av,bv,name);
      else if(name==='interfaceCells'||name==='t')assert.equal(av,bv,name);
      else{
        const absoluteTol=name==='m6'?1e-14:1e-24;
        assert.ok(Math.abs(av-bv)<=Math.max(absoluteTol,1e-11*Math.max(Math.abs(av),Math.abs(bv))),
          got.case+' step '+a.step+' '+name+' changed from '+bv+' to '+av);
      }
    }
  }
}
console.log('five-layer-baseline.test.mjs PASS · 3 reference solvers, '+checked+
  ' pinned SHA256 numerical fields, SI mass/latent heat/thickness and 12 physical steps each');
