import {
  CONSTANTS,
  activationBudget,
  findKohlerCriticalPoint,
  kohlerSaturationRatio,
  kohlerTerms,
  maxwellMasonRadiusRate,
  saturationVaporPressureWaterMK,
} from './thermo.js';

function bisect(fn, a, b, iterations = 70) {
  let fa = fn(a), fb = fn(b);
  if (fa === 0) return a;
  if (fb === 0) return b;
  if (fa * fb > 0) return null;
  for (let i = 0; i < iterations; i++) {
    const m = 0.5 * (a + b), fm = fn(m);
    if (Math.abs(fm) < 1e-13) return m;
    if (fa * fm <= 0) { b = m; fb = fm; } else { a = m; fa = fm; }
  }
  return 0.5 * (a + b);
}

export function findKohlerEquilibria(Dd, kappaHyg, T, Sinf, options = {}) {
  const critical = findKohlerCriticalPoint(Dd, kappaHyg, T, options);
  const minD = Dd * 1.0001;
  const maxD = Math.max(Dd * 2e4, Number.isFinite(critical.Dcrit) ? critical.Dcrit * 40 : 100e-6);
  const n = options.samples ?? 3000;
  const logMin = Math.log(minD), logMax = Math.log(maxD);
  const f = D => kohlerSaturationRatio(D, Dd, kappaHyg, T, options) - Sinf;
  const roots = [];
  let prevD = minD, prevF = f(prevD);
  for (let i = 1; i < n; i++) {
    const D = Math.exp(logMin + (logMax - logMin) * i / (n - 1));
    const curF = f(D);
    if (Number.isFinite(prevF) && Number.isFinite(curF) && prevF * curF < 0) {
      const root = bisect(f, prevD, D);
      if (root && !roots.some(r => Math.abs(r.D - root) / root < 1e-6)) {
        const h = root * 1e-4;
        const slope = (f(root + h) - f(Math.max(minD, root - h))) / (2 * h);
        roots.push({ D: root, slope, stable: slope > 0, type: slope > 0 ? 'stable' : 'unstable' });
      }
    }
    prevD = D; prevF = curF;
  }
  return { roots, critical, barrierGone: Sinf >= critical.Scrit };
}

export function initialEquilibriumDiameter(Dd, kappaHyg, T, Sinf, options = {}) {
  const eq = findKohlerEquilibria(Dd, kappaHyg, T, Sinf, options);
  return eq.roots.find(r => r.stable)?.D ?? Dd * 1.01;
}

export function makeActivationLandscape({ Dd, kappaHyg, T, Sinf, points = 420, maxFactor = 80, ...options }) {
  const eq = findKohlerEquilibria(Dd, kappaHyg, T, Sinf, options);
  const minD = Dd * 1.001;
  const maxD = Math.max(Dd * maxFactor, Number.isFinite(eq.critical.Dcrit) ? eq.critical.Dcrit * 5 : Dd * maxFactor);
  const logMin = Math.log(minD), logMax = Math.log(maxD);
  const rows = [];
  let W = 0;
  let prev = null;
  for (let i = 0; i < points; i++) {
    const D = Math.exp(logMin + (logMax - logMin) * i / (points - 1));
    const terms = kohlerTerms(D, Dd, kappaHyg, T, options);
    const lnRatio = Math.log(Math.max(terms.Seq, 1e-300) / Math.max(Sinf, 1e-300));
    const deltaMu = CONSTANTS.R * T * lnRatio;
    const dNdD = CONSTANTS.RHO_W * Math.PI * D * D / (2 * CONSTANTS.MW);
    if (prev) {
      const integrand = deltaMu * dNdD;
      W += 0.5 * (integrand + prev.integrand) * (D - prev.D);
    }
    rows.push({ D, Dnm: D * 1e9, deltaMu, W, lnRatio });
    prev = { D, integrand: deltaMu * dNdD };
  }
  const minW = Math.min(...rows.map(r => r.W));
  const kBT = 1.380649e-23 * T;
  for (const row of rows) {
    row.Wrel = row.W - minW;
    row.WkBT = row.Wrel / kBT;
  }

  let barrier = null;
  const stable = eq.roots.find(r => r.stable);
  const unstable = eq.roots.find(r => !r.stable);
  if (stable && unstable) {
    const nearest = D => rows.reduce((a, b) => Math.abs(b.D - D) < Math.abs(a.D - D) ? b : a);
    const lo = nearest(stable.D), hi = nearest(unstable.D);
    barrier = { stableD: stable.D, unstableD: unstable.D, heightJ: hi.W - lo.W, heightKBT: (hi.W - lo.W) / kBT };
  }
  return { rows, equilibria: eq, barrier };
}

export function createAmbientDriver(mode, ambient, config = {}) {
  const Sw0 = 1 + ambient.supersaturationPct / 100;
  return {
    mode,
    T0: ambient.T,
    Sw0,
    e0: Sw0 * saturationVaporPressureWaterMK(ambient.T),
    startSupersaturationPct: config.startSupersaturationPct ?? ambient.supersaturationPct,
    targetSupersaturationPct: config.targetSupersaturationPct ?? 0.45,
    rampDuration: config.rampDuration ?? 12,
    coolingRateKPerMin: config.coolingRateKPerMin ?? 4,
    targetTemperatureK: config.targetTemperatureK ?? ambient.T - 20,
  };
}

export function ambientFromDriver(driver, t) {
  if (driver.mode === 'supersat_ramp') {
    const f = Math.min(1, Math.max(0, t / Math.max(driver.rampDuration, 1e-6)));
    const sPct = driver.startSupersaturationPct + (driver.targetSupersaturationPct - driver.startSupersaturationPct) * f;
    return { T: driver.T0, supersaturationPct: sPct, Sw: 1 + sPct / 100 };
  }
  if (driver.mode === 'cool_fixed_e') {
    const T = Math.max(driver.targetTemperatureK, driver.T0 - driver.coolingRateKPerMin * t / 60);
    const Sw = driver.e0 / saturationVaporPressureWaterMK(T);
    return { T, supersaturationPct: (Sw - 1) * 100, Sw };
  }
  return { T: driver.T0, supersaturationPct: (driver.Sw0 - 1) * 100, Sw: driver.Sw0 };
}

export function integrateWetDiameter({ D, Dd, kappaHyg, T, Sinf, pressurePa, dt, maxFraction = 0.015 }) {
  let radius = Math.max(D / 2, Dd * 0.500001);
  let remaining = dt;
  let substeps = 0;
  let lastTransport = null;
  while (remaining > 1e-12 && substeps < 20000) {
    const terms = kohlerTerms(2 * radius, Dd, kappaHyg, T);
    const transport = maxwellMasonRadiusRate(radius, T, Sinf, terms.Seq, pressurePa);
    lastTransport = transport;
    const rate = transport.rate;
    const maxDr = Math.max(radius * maxFraction, 1e-12);
    const h = Math.min(remaining, Math.abs(rate) > 0 ? maxDr / Math.abs(rate) : remaining);
    radius = Math.max(Dd * 0.500001, radius + rate * h);
    remaining -= h;
    substeps++;
  }
  return { D: 2 * radius, substeps, transport: lastTransport };
}

export function activationSnapshot({ t, T, supersaturationPct, pressurePa, Dd, kappaHyg, D }) {
  const Sinf = 1 + supersaturationPct / 100;
  const critical = findKohlerCriticalPoint(Dd, kappaHyg, T);
  const equilibria = findKohlerEquilibria(Dd, kappaHyg, T, Sinf);
  const terms = kohlerTerms(D, Dd, kappaHyg, T);
  const budget = activationBudget(D, Dd, kappaHyg, T, Sinf);
  const transport = maxwellMasonRadiusRate(D / 2, T, Sinf, terms.Seq, pressurePa);
  return { t, T, supersaturationPct, pressurePa, Dd, kappaHyg, D, Sinf, critical, equilibria, terms, budget, transport };
}
