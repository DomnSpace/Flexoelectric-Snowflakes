import assert from 'node:assert/strict';
import { chromium } from 'playwright-core';
import { mkdir } from 'node:fs/promises';
import { resolve } from 'node:path';

const browser=await chromium.launch({
  headless:true,
  executablePath:process.env.CHROME_PATH||'/usr/bin/google-chrome',
  args:['--no-sandbox','--enable-webgl','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']
});
try{
  const page=await browser.newPage({viewport:{width:1550,height:1000},deviceScaleFactor:1,acceptDownloads:true});
  const errors=[];
  page.on('pageerror',err=>errors.push(err.message));
  page.on('console',msg=>{if(msg.type()==='error')errors.push(msg.text());});
  await page.goto('http://127.0.0.1:8088/',{waitUntil:'domcontentloaded'});
  await page.waitForFunction(()=>document.querySelector('#status')?.textContent!=='INITIALIZING PHYSICS',{timeout:45000});
  assert.equal(await page.locator('#errorOverlay').isVisible(),false,'3D/solver error overlay must be hidden');
  assert.match(await page.locator('#status').innerText(),/SOLVER RUNNING/);
  assert.equal(await page.locator('#viewport canvas').count(),1,'WebGL canvas mounted');
  assert.equal(await page.locator('#fieldCanvas').evaluate(el=>el.width),73,'numeric field raster');
  for(const id of ['radialGraph','growthGraph','modeGraph']){
    const dims=await page.locator('#'+id).evaluate(el=>[el.width,el.height]);
    assert.ok(dims[0]>20&&dims[1]>20,'chart '+id+' has dimensions');
  }
  await page.locator('#play').click();
  const t0=parseFloat(await page.locator('#kTime').innerText());
  await page.locator('#step').click();
  const t1=parseFloat(await page.locator('#kTime').innerText());
  assert.ok(t1>t0,'actual physics step must advance time');
  await page.locator('#sigma').evaluate(el=>{el.value='19';el.dispatchEvent(new Event('input',{bubbles:true}));});
  assert.equal(await page.locator('#pending').isVisible(),true,'parameter changes must be pending');
  await page.locator('#apply').click();
  assert.equal(await page.locator('#pending').isVisible(),false,'apply resets model');
  // Exercise the optional electrostatic experiment with a physically
  // computed potential, rather than trusting only the new DOM controls.
  await page.locator('#electricEnabled').check();
  await page.locator('#electricCouple').check();
  await page.locator('#eField').evaluate(el=>{el.value='10';el.dispatchEvent(new Event('input',{bubbles:true}));});
  await page.locator('#bendGradient').evaluate(el=>{el.value='1500';el.dispatchEvent(new Event('input',{bubbles:true}));});
  await page.locator('#apply').click();
  assert.match(await page.locator('#electricReadout').innerText(),/POISSON \+ DIELECTRIC/);
  await page.locator('#field').selectOption('electric');
  assert.match(await page.locator('#fieldMax').innerText(),/MV\/m/);
  await page.locator('#step').click();
  const maxE=parseFloat(await page.locator('#bEField').innerText());
  assert.ok(maxE>0,'physical electric field must be nonzero');
  assert.equal(await page.locator('#errorOverlay').isVisible(),false,'electric solver must not throw');
  await page.locator('#field').selectOption('potential');
  assert.match(await page.locator('#fieldMin').innerText(),/ V/);
  await page.locator('#field').selectOption('temperature');
  await page.locator('#scrub').evaluate(el=>{el.value='0';el.dispatchEvent(new Event('input',{bubbles:true}));});
  const pngValid=await page.locator('#viewport canvas').evaluate(el=>el.toDataURL('image/png').startsWith('data:image/png;base64,'));
  assert.ok(pngValid,'WebGL canvas PNG export works');
  const downloadPromise=page.waitForEvent('download');
  await page.locator('#csv').click();
  const download=await downloadPromise;
  assert.match(download.suggestedFilename(),/\.csv$/);
  // Preset switches are observational hypotheses, not geometry shortcuts.
  await page.locator('#habit').selectOption('dendrite');
  await page.locator('#preset').click();
  assert.equal(await page.locator('#fieldCanvas').evaluate(el=>el.width),129,
    'dendritic experiment uses a larger solver field, not zoom');
  assert.match(await page.locator('#scenario').innerText(),/STELLAR DENDRITE CANDIDATE/);
  assert.match(await page.locator('#waterLimit').innerText(),/above ice/);
  assert.match(await page.locator('#domainDetails').innerText(),/HALO/);
  assert.equal(await page.locator('#electricEnabled').isChecked(),false,
    'atmospheric habit presets must not silently energize the electric model');
  await page.locator('#habit').selectOption('needle');
  assert.equal(await page.locator('#preset').isDisabled(),true,
    'unsupported needle must not manufacture an extruded column');
  assert.match(await page.locator('#habitNote').innerText(),/not simulated/i);
  await page.locator('#habit').selectOption('plate');
  await page.locator('#preset').click();
  assert.equal(await page.locator('#fieldCanvas').evaluate(el=>el.width),97);
  assert.equal(Number(await page.locator('#basalA').inputValue()),.009);
  assert.match(await page.locator('#domainStatus').innerText(),/VALIDITY|CAUTION/);
  assert.equal(await page.locator('#errorOverlay').isVisible(),false);
  await mkdir(resolve('artifacts'),{recursive:true});
  await page.screenshot({path:'artifacts/snowflake-observatory.png',fullPage:true,animations:'disabled'});
  assert.deepEqual(errors,[],'no browser script/console errors');
  console.log('observatory-browser.test.mjs PASS',{t0,t1,download:download.suggestedFilename()});
}finally{
  await browser.close();
}
