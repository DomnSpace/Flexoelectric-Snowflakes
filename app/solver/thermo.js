/**
 * Milestone-1 thermodynamics for the activation laboratory.
 *
 * Primary references:
 * - Petters & Kreidenweis (2007), ACP 7, 1961–1971, doi:10.5194/acp-7-1961-2007
 * - Murphy & Koop (2005), QJRMS 131, 1539–1565, doi:10.1256/qj.04.94
 * - IAPWS R1-76(2014), surface tension of ordinary water.
 *
 * Scientific status:
 * - saturation vapor pressure: measured-fit / reviewed thermodynamics
 * - kappa-Kohler: semi-empirical single-parameter water-activity representation
 * - surface tension: IAPWS equilibrium formulation; metastable extrapolation only to its stated range
 */

export const CONSTANTS = Object.freeze({
  R: 8.31446261815324,
  RV: 461.52280830495,
  MW: 0.01801528,
  RHO_W: 997.0,
  TC: 647.096,
  P0: 101325,
  D_V0: 2.2e-5,
  K_AIR0: 0.0241,
});

export function clamp(x, lo, hi) { return Math.max(lo, Math.min(hi, x)); }

export function saturationVaporPressureWaterMK(T) {
  if (!(T > 0)) throw new RangeError('Temperature must be in kelvin.');
  const lnT = Math.log(T);
  const lnP = 54.842763 - 6763.22 / T - 4.210 * lnT + 0.000367 * T
    + Math.tanh(0.0415 * (T - 218.8))
      * (53.878 - 1331.22 / T - 9.44523 * lnT + 0.014025 * T);
  return Math.exp(lnP);
}

export function saturationVaporPressureIceMK(T) {
  if (!(T > 110)) throw new RangeError('Murphy–Koop ice fit requires T > 110 K.');
  return Math.exp(9.550426 - 5723.265 / T + 3.53068 * Math.log(T) - 0.00728332 * T);
}

export function surfaceTensionWaterIAPWS(T) {
  const tau = 1 - T / CONSTANTS.TC;
  if (tau <= 0) return 0;
  return 0.2358 * Math.pow(tau, 1.256) * (1 - 0.625 * tau);
}

export function kelvinLength(T, sigma = surfaceTensionWaterIAPWS(T), rhoW = CONSTANTS.RHO_W) {
  return 4 * sigma * CONSTANTS.MW / (CONSTANTS.R * T * rhoW);
}

export function waterActivityKappa(D, Dd, kappaHyg) {
  if (!(D > Dd) || !(Dd > 0)) return 0;
  const D3 = D ** 3, Dd3 = Dd ** 3;
  return clamp((D3 - Dd3) / (D3 - Dd3 * (1 - kappaHyg)), 0, 1);
}

export function kohlerTerms(D, Dd, kappaHyg, T, options = {}) {
  const sigma = options.sigma ?? surfaceTensionWaterIAPWS(T);
  const rhoW = options.rhoW ?? CONSTANTS.RHO_W;
  const A = kelvinLength(T, sigma, rhoW);
  const aw = waterActivityKappa(D, Dd, kappaHyg);
  const lnAw = aw > 0 ? Math.log(aw) : -Infinity;
  const lnKelvin = A / D;
  const lnSeq = lnAw + lnKelvin;
  return {
    aw, A, sigma, lnAw, lnKelvin, lnSeq,
    Seq: Number.isFinite(lnSeq) ? Math.exp(lnSeq) : 0,
  };
}

export function kohlerSaturationRatio(D, Dd, kappaHyg, T, options = {}) {
  return kohlerTerms(D, Dd, kappaHyg, T, options).Seq;
}

function goldenSectionMax(fn, a, b, iterations = 64) {
  const gr = (Math.sqrt(5) - 1) / 2;
  let c = b - gr * (b - a), d = a + gr * (b - a);
  let fc = fn(c), fd = fn(d);
  for (let i = 0; i < iterations; i++) {
    if (fc > fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = fn(c); }
    else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = fn(d); }
  }
  const x = (a + b) / 2;
  return { x, y: fn(x) };
}

export function findKohlerCriticalPoint(Dd, kappaHyg, T, options = {}) {
  if (!(Dd > 0) || !(kappaHyg >= 0)) throw new RangeError('Invalid dry diameter or kappa.');
  if (kappaHyg === 0) return { Dcrit: Infinity, Scrit: 1, sCritPct: 0, numerical: true };
  const minD = Dd * (1 + 1e-7);
  const maxD = Math.max(Dd * 1e4, 50e-6);
  const logMin = Math.log(minD), logMax = Math.log(maxD);
  const n = options.samples ?? 2048;
  let bestI = 1, bestS = -Infinity;
  const Ds = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const D = Math.exp(logMin + (logMax - logMin) * i / (n - 1));
    Ds[i] = D;
    const S = kohlerSaturationRatio(D, Dd, kappaHyg, T, options);
    if (S > bestS) { bestS = S; bestI = i; }
  }
  const lo = Ds[Math.max(0, bestI - 1)], hi = Ds[Math.min(n - 1, bestI + 1)];
  const refined = goldenSectionMax(D => kohlerSaturationRatio(D, Dd, kappaHyg, T, options), lo, hi);
  return { Dcrit: refined.x, Scrit: refined.y, sCritPct: (refined.y - 1) * 100, numerical: true };
}

export function kohlerCriticalApproximation(Dd, kappaHyg, T, options = {}) {
  const A = kelvinLength(T, options.sigma ?? surfaceTensionWaterIAPWS(T), options.rhoW ?? CONSTANTS.RHO_W);
  if (kappaHyg <= 0) return { Dcrit: Infinity, Scrit: 1, sCritPct: 0 };
  const Dcrit = Math.sqrt(3 * kappaHyg * Dd ** 3 / A);
  const s = Math.sqrt(4 * A ** 3 / (27 * kappaHyg * Dd ** 3));
  return { Dcrit, Scrit: 1 + s, sCritPct: 100 * s };
}

export function saturationState(T, ambientSw) {
  const esw = saturationVaporPressureWaterMK(T);
  const esi = saturationVaporPressureIceMK(Math.max(110.01, T));
  const e = ambientSw * esw;
  return { esw, esi, e, Sw: ambientSw, Si: e / esi };
}

export function sphericalVaporField(r, radius, Ssurface, Sinf) {
  if (r <= radius) return Ssurface;
  return Sinf + (Ssurface - Sinf) * radius / r;
}

export function sphericalSaturationGradient(r, radius, Ssurface, Sinf) {
  const rr = Math.max(radius, r);
  return Math.abs((Ssurface - Sinf) * radius / (rr * rr));
}

export function diffusivityWaterVaporAir(T, pressurePa = CONSTANTS.P0) {
  return CONSTANTS.D_V0 * Math.pow(T / 273.15, 1.94) * (CONSTANTS.P0 / pressurePa);
}

export function thermalConductivityAir(T) {
  return CONSTANTS.K_AIR0 * Math.pow(T / 273.15, 0.90);
}

export function latentHeatVaporization(T) {
  return 2.501e6 - 2.37e3 * (T - 273.15);
}

export function maxwellMasonRadiusRate(radius, T, Sinf, Ssurface, pressurePa = CONSTANTS.P0, options = {}) {
  const rhoW = options.rhoW ?? CONSTANTS.RHO_W;
  const Dv = options.Dv ?? diffusivityWaterVaporAir(T, pressurePa);
  const Kair = options.Kair ?? thermalConductivityAir(T);
  const Lv = options.Lv ?? latentHeatVaporization(T);
  const esw = saturationVaporPressureWaterMK(T);
  const Fd = rhoW * CONSTANTS.RV * T / (Dv * esw);
  const Fk = rhoW * Lv / (Kair * T) * (Lv / (CONSTANTS.RV * T) - 1);
  const denominator = Math.max(radius, 1e-12) * (Fd + Fk);
  return { rate: (Sinf - Ssurface) / denominator, Fd, Fk, Dv, Kair, Lv };
}

export function dropletRadiusRate(radius, T, Sinf, Ssurface, options = {}) {
  return maxwellMasonRadiusRate(radius, T, Sinf, Ssurface, options.pressurePa ?? CONSTANTS.P0, options).rate;
}

export function activationBudget(D, Dd, kappaHyg, T, Sinf, options = {}) {
  const terms = kohlerTerms(D, Dd, kappaHyg, T, options);
  const ambient = Math.log(Math.max(Sinf, 1e-12));
  const hygroscopic = -terms.lnAw;
  const kelvin = -terms.lnKelvin;
  const net = ambient + hygroscopic + kelvin;
  const deltaMuTransfer = -CONSTANTS.R * T * net;
  return { ambient, hygroscopic, kelvin, net, deltaMuTransfer, terms };
}

export function makeKohlerCurve({ Dd, kappaHyg, T, points = 420, maxFactor = 60, ...options }) {
  const critical = findKohlerCriticalPoint(Dd, kappaHyg, T, options);
  const minD = Dd * 1.001;
  const maxD = Math.max(Dd * maxFactor, Number.isFinite(critical.Dcrit) ? critical.Dcrit * 4 : Dd * maxFactor);
  const logMin = Math.log(minD), logMax = Math.log(maxD);
  const data = [];
  for (let i = 0; i < points; i++) {
    const D = Math.exp(logMin + (logMax - logMin) * i / (points - 1));
    const t = kohlerTerms(D, Dd, kappaHyg, T, options);
    data.push({
      D, Dnm: D * 1e9,
      totalPct: (t.Seq - 1) * 100,
      activityPct: (t.aw - 1) * 100,
      kelvinPct: (Math.exp(t.lnKelvin) - 1) * 100,
      lnTotal: t.lnSeq,
    });
  }
  return { data, critical };
}
