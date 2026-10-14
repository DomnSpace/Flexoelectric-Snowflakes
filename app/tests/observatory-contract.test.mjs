import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, resolve } from 'node:path';
const dir=resolve(dirname(fileURLToPath(import.meta.url)),'../..');
const html=readFileSync(resolve(dir,'index.html'),'utf8');
const script=readFileSync(resolve(dir,'app/observatory.js'),'utf8');
const css=readFileSync(resolve(dir,'app/observatory.css'),'utf8');
const ids=new Set([...html.matchAll(/\bid="([^"]+)"/g)].map(x=>x[1]));
const accesses=[...script.matchAll(/\$\('([^']+)'\)/g)].map(x=>x[1]);
for(const id of accesses)assert.ok(ids.has(id),'observatory JS references nonexistent #'+id);
for(const id of ['viewport','play','step','reset','scrub','fieldCanvas','radialGraph','growthGraph','modeGraph','csv','png']){
  assert.ok(ids.has(id),'missing critical UI element #'+id);
}
assert.match(html,/type="module" src="\.\/app\/observatory\.js"/,'root HTML must run module');
assert.match(html,/three\/addons\//,'pinned Three.js add-on import map must be present');
assert.match(script,/from '\.\/solver\/growth25d\.js'/,'must import physical growth solver');
assert.match(script,/resolvedGrowthStep\(state,STEPTIME/,'must step physical kernel');
assert.match(script,/u:new Float32Array\(state\.sigma\)/,'must capture vapor grid');
assert.match(script,/dT:new Float32Array\(state\.dT\)/,'must capture thermal grid');
assert.match(script,/phi:new Float32Array\(state\.phi\)/,'must capture ice field');
assert.ok(!/Math\.random\(/.test(script),'no synthetic scientific data');
assert.ok(css.includes('.graph-grid')&&css.includes('.viewport'));
console.log('observatory-contract.test.mjs PASS', {uiIds:ids.size,refs:accesses.length});
