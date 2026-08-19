import assert from 'node:assert/strict';
import { iceIhProperties } from '../solver/ice-ih-eos.js';

// IAPWS R10-06(2009) verification-style sanity points.
const a=iceIhProperties(273.16,611.657);
assert.ok(a.rho>916 && a.rho<918);
assert.ok(a.cp>1900 && a.cp<2300);
assert.ok(Number.isFinite(a.g));
assert.ok(Number.isFinite(a.h));
assert.ok(a.kappaT>0);

const b=iceIhProperties(250,101325);
assert.ok(b.rho>915 && b.rho<925);
assert.ok(b.cp>1500 && b.cp<2300);
assert.ok(b.h<a.h);

const c=iceIhProperties(200,100e6);
assert.ok(c.rho>920);
assert.ok(c.kappaT>0);
assert.ok(Number.isFinite(c.s));

console.log('Ice Ih EOS tests passed');
