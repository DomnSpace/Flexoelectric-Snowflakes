/**
 * IAPWS G12-15 supercooled-water Gibbs-energy equation of state.
 * Reference: IAPWS, Guideline on Thermodynamic Properties of Supercooled Water (2015).
 *
 * Density is evaluated from the analytic pressure derivative in Eq. (13).
 * cp, alpha_p, and kappa_T are evaluated by centered derivatives of the same
 * fundamental Gibbs/volume functions. This keeps every returned property tied
 * to one thermodynamic potential while avoiding duplicated derivative algebra.
 */

const EOS = {
  omega0: 0.5212269,
  L0: 0.76317954,
  k0: 0.072158686,
  k1: -0.31569232,
  k2: 5.2992608,
  TLL: 228.2,
  rho0: 1081.6482,
  R: 461.523087,
};
EOS.pi0 = 300e6 / (EOS.rho0 * EOS.R * EOS.TLL);

const C = [
  [-8.1570681381655, 0, 0, 0],
  [1.2875032, 0, 1, 0],
  [7.0901673598012, 1, 0, 0],
  [-3.2779161e-2, -0.2555, 2.1051, -0.0016],
  [7.3703949e-1, 1.5762, 1.1422, 0.6894],
  [-2.1628622e-1, 1.6400, 0.9510, 0.0130],
  [-5.1782479, 3.6385, 0, 0.0002],
  [4.2293517e-4, -0.3828, 3.6402, 0.0435],
  [2.3592109e-2, 1.6219, 2.0760, 0.0500],
  [4.3773754, 4.3287, -0.0016, 0.0004],
  [-2.9967770e-3, 3.4763, 2.2769, 0.0528],
  [-9.6558018e-1, 5.1556, 0.0008, 0.0147],
  [3.7595286, -0.3593, 0.3706, 0.8584],
  [1.2632441, 5.0361, -0.3975, 0.9924],
  [2.8542697e-1, 2.9786, 2.9730, 1.0041],
  [-8.5994947e-1, 6.2373, -0.3180, 1.0961],
  [-3.2916153e-1, 4.0460, 2.9805, 1.0228],
  [9.0019616e-2, 5.3558, 2.9265, 1.0303],
  [8.1149726e-2, 9.0157, 0.4456, 1.6180],
  [-3.2788213, 1.2194, 0.1298, 0.5213],
];

function clamp(x, lo, hi) { return Math.max(lo, Math.min(hi, x)); }

function solveX(L, omega) {
  const c1 = (10 / 9) * (Math.log(19) - L);
  const c2 = (50 / 49) * (Math.log(99) - L);
  let lo, hi;
  if (omega < c1) { lo = 0.049; hi = 0.5; }
  else if (omega < c2) { lo = 0.0099; hi = 0.051; }
  else {
    lo = 0.99 * Math.exp(-(50 / 49) * L - omega);
    hi = Math.min(1.1 * Math.exp(-L - omega), 0.0101);
  }
  lo = clamp(lo, 1e-14, 1 - 1e-14);
  hi = clamp(hi, lo + 1e-14, 1 - 1e-14);
  const f = x => L + Math.log(x / (1 - x)) + omega * (1 - 2 * x);
  let flo = f(lo), fhi = f(hi);
  if (!(flo * fhi <= 0)) {
    lo = 1e-12; hi = 0.5; flo = f(lo); fhi = f(hi);
  }
  for (let i = 0; i < 100; i++) {
    const mid = 0.5 * (lo + hi), fm = f(mid);
    if (flo * fm <= 0) { hi = mid; fhi = fm; }
    else { lo = mid; flo = fm; }
  }
  return 0.5 * (lo + hi);
}

function core(T, p) {
  const { omega0, L0, k0, k1, k2, TLL, rho0, R } = EOS;
  const tau = T / TLL - 1;
  const pi = p / (rho0 * R * TLL);
  const tauBar = tau + 1;
  const piBar = pi + EOS.pi0;
  const K2 = Math.sqrt(1 + k2 * k2);
  const q = pi - k2 * tau;
  const A = 1 + k0 * k2 + k1 * q;
  const K1 = Math.sqrt(A * A - 4 * k0 * k1 * k2 * q);
  const L = L0 * K2 / (2 * k1 * k2) * (1 + k0 * k2 + k1 * (pi + k2 * tau) - K1);
  const omega = 2 + omega0 * pi;
  const x = solveX(L, omega);
  const phi = 2 * x - 1;

  const Lpi = L0 * K2 * (K1 + k0 * k2 - k1 * pi + k1 * k2 * tau - 1) / (2 * k2 * K1);
  let psiR = 0, psiPi = 0;
  for (const [c, a, b, d] of C) {
    const e = Math.exp(-d * piBar);
    psiR += c * Math.pow(tauBar, a) * Math.pow(piBar, b) * e;
    psiPi += c * Math.pow(tauBar, a) * Math.pow(piBar, b - 1) * (b - d * piBar) * e;
  }

  const mix = x * L + x * Math.log(x) + (1 - x) * Math.log(1 - x) + omega * x * (1 - x);
  const psi = psiR + tauBar * mix;
  const B = omega0 / 2 * (1 - phi * phi) + Lpi * (phi + 1);
  const v = (tauBar / 2 * B + psiPi) / rho0;
  return { T, p, tau, pi, L, omega, x, phi, psi, g: psi * R * TLL, v, rho: 1 / v };
}

export function supercooledWaterCore(T, pPa = 101325) {
  return core(T, pPa);
}

export function homogeneousNucleationPressureLimit(T) {
  const theta = T / 235.15;
  return (0.1 + 228.27 * (1 - Math.pow(theta, 6.243)) + 15.724 * (1 - Math.pow(theta, 79.81))) * 1e6;
}

export function homogeneousNucleationTemperatureLimit(pPa = 101325) {
  const pMPa = pPa / 1e6;
  if (pMPa > 198.9) {
    return 172.82 + 0.03718 * pMPa + 3.403e-5 * pMPa ** 2 - 1.573e-8 * pMPa ** 3;
  }
  let lo = 175, hi = 240;
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    const pm = homogeneousNucleationPressureLimit(mid) / 1e6;
    if (pm > pMPa) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

export function supercooledWaterProperties(T, pPa = 101325) {
  const c = core(T, pPa);
  const hT = 0.01;
  const gp = core(T + hT, pPa).g;
  const gm = core(T - hT, pPa).g;
  const cp = -T * (gp - 2 * c.g + gm) / (hT * hT);

  const vpT = core(T + hT, pPa).v;
  const vmT = core(T - hT, pPa).v;
  const alphaP = (vpT - vmT) / (2 * hT * c.v);

  const hp = Math.max(100, Math.abs(pPa) * 1e-6);
  const vpP = core(T, pPa + hp).v;
  const vmP = core(T, Math.max(1, pPa - hp)).v;
  const kappaT = -(vpP - vmP) / (2 * hp * c.v);

  const cv = cp - T * alphaP * alphaP / (c.rho * kappaT);
  const speedSound = Math.pow(c.rho * kappaT * cv / cp, -0.5);
  const TH = homogeneousNucleationTemperatureLimit(pPa);
  return {
    ...c,
    cp,
    cv,
    alphaP,
    kappaT,
    speedSound,
    homogeneousNucleationLimitK: TH,
    insideRecommendedMetastableRange: T >= TH && T <= 300 && pPa >= 0 && pPa <= 400e6,
  };
}

export const IAPWS_G12_15 = Object.freeze({ ...EOS, coefficients: C.map(row => [...row]) });
