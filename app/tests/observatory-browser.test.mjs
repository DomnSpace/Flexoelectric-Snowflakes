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
  await page.locator('#field').selectOption('temperature');
  await page.locator('#scrub').fill('0');
  const pngValid=await page.locator('#viewport canvas').evaluate(el=>el.toDataURL('image/png').startsWith('data:image/png;base64,'));
  assert.ok(pngValid,'WebGL canvas PNG export works');
  const downloadPromise=page.waitForEvent('download');
  await page.locator('#csv').click();
  const download=await downloadPromise;
  assert.match(download.suggestedFilename(),/\.csv$/);
  await mkdir(resolve('artifacts'),{recursive:true});
  await page.screenshot({path:'artifacts/snowflake-observatory.png',fullPage:true,animations:'disabled'});
  assert.deepEqual(errors,[],'no browser script/console errors');
  console.log('observatory-browser.test.mjs PASS',{t0,t1,download:download.suggestedFilename()});
}finally{
  await browser.close();
}
