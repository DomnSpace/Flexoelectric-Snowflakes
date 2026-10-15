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
  assert.equal(await page.locator('#errorOverlay').isVisible(),false,
    '3D/solver error overlay must be hidden: '+await page.locator('#errorOverlay').innerText()+
    ' / console: '+errors.join(' | '));
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
  assert.equal(await page.locator('#topology').inputValue(),'triangular');
  assert.equal(await page.locator('#gridSize').inputValue(),'129');
  assert.equal(await page.locator('#fieldCanvas').evaluate(el=>el.width),1024,
    'triangular field uses display-only 1024-pixel interpolation, not solver-cell raster');
  assert.match(await page.locator('#fieldGeometry').innerText(),/TRIANGULAR TRANSPORT/);
  assert.match(await page.locator('#gridDim').innerText(),/HEX CELLS/);
  assert.equal(await page.locator('#electricEnabled').isDisabled(),true,
    'Poisson solver is not falsely claimed to support hex mesh');
  await page.locator('#play').click();
  const hexT0=parseFloat(await page.locator('#kTime').innerText());
  await page.locator('#step').click();
  assert.ok(parseFloat(await page.locator('#kTime').innerText())>hexT0,'actual hex physics advances');
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
  assert.equal(await page.locator('#gridSize').inputValue(),'97');
  assert.equal(await page.locator('#fieldCanvas').evaluate(el=>el.width),1024);
  assert.equal(await page.locator('#topology').inputValue(),'triangular');
  assert.equal(Number(await page.locator('#basalA').inputValue()),.009);
  assert.match(await page.locator('#domainStatus').innerText(),/VALIDITY|CAUTION/);
  assert.equal(await page.locator('#errorOverlay').isVisible(),false);
  // Controlled 30° test must reach the HEX PHYSICS, not just its labels.
  await page.locator('#kineticPhase').selectOption('30');
  await page.locator('#surfaceEpsilon').selectOption('.02');
  await page.locator('#seedShape').selectOption('hex');
  await page.locator('#seedPhase').selectOption('30');
  assert.equal(await page.locator('#pending').isVisible(),true,
    'changing a physical orientation requires a solver restart');
  await page.locator('#apply').click();
  assert.equal(await page.locator('#errorOverlay').isVisible(),false);
  assert.equal((await page.locator('#bPhase').innerText()).trim(),'30°');
  assert.equal((await page.locator('#bStiffness').innerText()).trim(),'0.020');
  await page.locator('#step').click();
  assert.equal(await page.locator('#errorOverlay').isVisible(),false,
    'anisotropic surface stiffness and hexagonal seed must not break growth');
  const phaseCsv=page.waitForEvent('download');
  await page.locator('#csv').click();
  const phaseDownload=await phaseCsv;
  assert.match(phaseDownload.suggestedFilename(),/\.csv$/);
    await page.locator('#topology').selectOption('square');
  await page.locator('#apply').click();
  assert.equal(await page.locator('#fieldCanvas').evaluate(el=>el.width),97,
    'Cartesian baseline remains independently functional for A/B tests');
  assert.equal(await page.locator('#electricEnabled').isDisabled(),false);
  await page.locator('#topology').selectOption('triangular');
  await page.locator('#apply').click();
  assert.equal(await page.locator('#errorOverlay').isVisible(),false);
  // Independent high-definition raster: changes display pixels, never
  // solver grid or elapsed model time.
  await page.locator('#fieldResolution').selectOption('1536');
  assert.equal(await page.locator('#fieldCanvas').evaluate(el=>el.width),1536);
  assert.match(await page.locator('#fieldFidelity').innerText(),/PIECEWISE-LINEAR/);
  await page.locator('#fieldResolution').selectOption('1024');
  // Full physical 300 µm hex-domain run at the existing 0.75 µm spacing.
  await page.locator('#gridSize').selectOption('401');
  await page.locator('#apply').click();
  assert.equal(await page.locator('#fieldCanvas').evaluate(el=>el.width),1024);
  assert.match(await page.locator('#gridDim').innerText(),/120601 HEX CELLS/);
  assert.match(await page.locator('#sizeLabel').innerText(),/300\.0 µm DOMAIN/);
  await page.locator('#play').click();
  const largeBefore=parseFloat(await page.locator('#kTime').innerText());
  const largeWallStart=Date.now();
  await page.locator('#step').click({timeout:60000});
  const largeStepMs=Date.now()-largeWallStart;
  assert.ok(parseFloat(await page.locator('#kTime').innerText())>largeBefore,
    '300 µm simulation must advance physics rather than simply display a larger frame');
  assert.ok(largeStepMs<50000,'401-node physical step stays interactive within smoke budget');
  assert.equal(await page.locator('#errorOverlay').isVisible(),false,
    '401-node hex simulation must not fail');
  const maxReplay=Number(await page.locator('#scrub').getAttribute('max'));
  assert.ok(maxReplay<=20,'large-grid replay remains memory bounded');
  await mkdir(resolve('artifacts'),{recursive:true});
  await page.screenshot({path:'artifacts/snowflake-observatory.png',fullPage:true,animations:'disabled'});
  assert.deepEqual(errors,[],'no browser script/console errors');
  console.log('observatory-browser.test.mjs PASS',{t0,t1,largeStepMs,download:download.suggestedFilename()});
}finally{
  await browser.close();
}
