/**
 * Milestone-1 thermodynamics and spherical condensational growth.
 *
 * References:
 * - Petters & Kreidenweis (2007), ACP 7, 1961–1971, doi:10.5194/acp-7-1961-2007
 * - Murphy & Koop (2005), QJRMS 131, 1539–1565, doi:10.1256/qj.04.94
 * - IAPWS R1-76(2014), surface tension of ordinary water.
 *
 * Scientific status:
 * - saturation vapor pressures: measured-fit / reviewed thermodynamics
 * - kappa-Kohler: semi-empirical single-parameter activity representation
 * - droplet growth below: transport approximation (quasi-steady vapor diffusion only;
 *   latent-heat resistance is deliberately deferred to the coupled field solver)
 */

export const CONSTANTS = Object.freeze({
  R: 8.31446261815324,      // J mol^-1 K^-1
  MW: 0.01801528,           // kg mol^-1
  RHO_W: 997.0,             // kg m^-3, fixed Milestone-1 approximation
  TC: 647.096,              // K, water critical temperature
  D_V: 2.2e-5,              // m^2 s^-1, nominal water-vapor diffusivity in air near ambient
});

export function clamp(x, lo, hi) {
  return Math.max(lo, Math.min(hi, x));
}

export function saturationVaporPressureWaterMK(T) {
  if (!(T > 0)) throw new RangeError('Temperature must be in kelvin.');
  const lnT = Math.log(T);
  const lnP = 54.842763 - 6763.22 / T - 4.210 * lnT + 0.000367 * T
    + Math.tanh(0.0415 * (T - 218.8))
      * (53.878 - 1331.22 / T - 9.44523 * lnT + 0.014025 * T);
  return Math.exp(lnP); // Pa
}

export function saturationVaporPressureIceMK(T) {
  if (!(T > 110)) throw new RangeError('Murphy–Koop ice fit requires T > 110 K.');
  return Math.exp(9.550426 - 5723.265 / T + 3.53068 * Math.log(T) - 0.00728332 * T); // Pa
}

/** IAPWS surface tension, N m^-1. Recommended on saturation curve; reasonable metastable extrapolation to ~248 K. */
export function surfaceTensionWaterIAPWS(T) {
  const tau = 1 - T / CONSTANTS.TC;
  if (tau <= 0) return 0;
  return 0.2358 * Math.pow(tau, 1.256) * (1 - 0.625 * tau);
}

export function kelvinLength(T, sigma = surfaceTensionWaterIAPWS(T), rhoW = CONSTANTS.RHO_W) {
  return 4 * sigma * CONSTANTS.MW / (CONSTANTS.R * T * rhoW); // m; diameter-form A
}

export function waterActivityKappa(D, Dd, kappaHyg) {
  if (!(D > Dd) || !(Dd > 0)) return 0;
  const D3 = D ** 3;
  const Dd3 = Dd ** 3;
  const numerator = D3 - Dd3;
  const denominator = D3 - Dd3 * (1 - kappaHyg);
  return clamp(numerator / denominator, 0, 1);
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
    aw,
    A,
    sigma,
    lnAw,
    lnKelvin,
    lnSeq,
    Seq: Number.isFinite(lnSeq) ? Math.exp(lnSeq) : 0,
  };
}

export function kohlerSaturationRatio(D, Dd, kappaHyg, T, options = {}) {
  return kohlerTerms(D, Dd, kappaHyg, T, options).Seq;
}

function goldenSectionMax(fn, a, b, iterations = 60) {
  const gr = (Math.sqrt(5) - 1) / 2;
  let c = b - gr * (b - a);
  let d = a + gr * (b - a);
  let fc = fn(c), fd = fn(d);
  for (let i = 0; i < iterations; i++) {
    if (fc > fd) {
      b = d; d = c; fd = fc; c = b - gr * (b - a); fc = fn(c);
    } else {
      a = c; c = d; fc = fd; d = a + gr * (b - a); fd = fn(d);
    }
  }
  const x = (a + b) / 2;
  return { x, y: fn(x) };
}

/** Numerically locate the maximum of the full kappa-Kohler curve; no large-particle approximation in the solver. */
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
  const lo = Ds[Math.max(0, bestI - 1)];
  const hi = Ds[Math.min(n - 1, bestI + 1)];
  const refined = goldenSectionMax(
    D => kohlerSaturationRatio(D, Dd, kappaHyg, T, options), lo, hi
  );
  return {
    Dcrit: refined.x,
    Scrit: refined.y,
    sCritPct: (refined.y - 1) * 100,
    numerical: true,
  };
}

/** Standard large-particle approximation retained only as a validation diagnostic. */
export function kohlerCriticalApproximation(Dd, kappaHyg, T, options = {}) {
  const A = kelvinLength(T, options.sigma ?? surfaceTensionWaterIAPWS(T), options.rhoW ?? CONSTANTS.RHO_W);
  if (kappaHyg <= 0) return { Dcrit: Infinity, Scrit: 1, sCritPct: 0 };
  const Dcrit = Math.sqrt(3 * kappaHyg * Dd ** 3 / A);
  const s = Math.sqrt(4 * A ** 3 / (27 * kappaHyg * Dd ** 3)); // approx S-1
  return { Dcrit, Scrit: 1 + s, sCritPct: 100 * s };
}

export function saturationState(T, ambientSw) {
  const esw = saturationVaporPressureWaterMK(T);
  const esi = saturationVaporPressureIceMK(Math.max(110.01, T));
  const e = ambientSw * esw;
  return { esw, esi, e, Sw: ambientSw, Si: e / esi };
}

/**
 * Quasi-steady spherical vapor diffusion field outside the wet particle.
 * S(r)=S_inf+(S_surface-S_inf) R/r.
 */
export function sphericalVaporField(r, radius, Ssurface, Sinf) {
  if (r <= radius) return Ssurface;
  return Sinf + (Ssurface - Sinf) * radius / r;
}

/** Magnitude proportional to |dS/dr|; returns saturation-ratio gradient, m^-1. */
export function sphericalSaturationGradient(r, radius, Ssurface, Sinf) {
  const rr = Math.max(radius, r);
  return Math.abs((Ssurface - Sinf) * radius / (rr * rr));
}

/**
 * Vapor-diffusion-only condensational growth of droplet radius.
 * Derived from steady Fick flux to a sphere at fixed T.
 */
export function dropletRadiusRate(radius, T, Sinf, Ssurface, options = {}) {
  const Dv = options.Dv ?? CONSTANTS.D_V;
  const rhoW = options.rhoW ?? CONSTANTS.RHO_W;
  const esw = saturationVaporPressureWaterMK(T);
  const coeff = Dv * CONSTANTS.MW * esw / (rhoW * CONSTANTS.R * T);
  return coeff * (Sinf - Ssurface) / Math.max(radius, 1e-12); // m s^-1
}

export function activationBudget(D, Dd, kappaHyg, T, Sinf, options = {}) {
  const terms = kohlerTerms(D, Dd, kappaHyg, T, options);
  const ambient = Math.log(Math.max(Sinf, 1e-12));
  const hygroscopic = -terms.lnAw;
  const kelvin = -terms.lnKelvin;
  const net = ambient + hygroscopic + kelvin; // ln(Sinf/Seq)
  return { ambient, hygroscopic, kelvin, net, terms };
}

export function makeKohlerCurve({ Dd, kappaHyg, T, points = 360, maxFactor = 60, ...options }) {
  const critical = findKohlerCriticalPoint(Dd, kappaHyg, T, options);
  const minD = Dd * 1.001;
  const maxD = Math.max(Dd * maxFactor, Number.isFinite(critical.Dcrit) ? critical.Dcrit * 4 : Dd * maxFactor);
  const logMin = Math.log(minD), logMax = Math.log(maxD);
  const data = [];
  for (let i = 0; i < points; i++) {
    const D = Math.exp(logMin + (logMax - logMin) * i / (points - 1));
    const t = kohlerTerms(D, Dd, kappaHyg, T, options);
    data.push({
      D,
      Dnm: D * 1e9,
      totalPct: (t.Seq - 1) * 100,
      activityPct: (t.aw - 1) * 100,
      kelvinPct: (Math.exp(t.lnKelvin) - 1) * 100,
      lnTotal: t.lnSeq,
    });
  }
  return { data, critical };
}
