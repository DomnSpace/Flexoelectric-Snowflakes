/**
 * IAPWS R10-06(2009) — Revised Release on the Equation of State 2006 for H2O Ice Ih.
 * Fundamental Gibbs potential g(T,p), with the IAPWS-95 entropy reference for phase equilibria.
 *
 * The implementation follows Eq. (1), Tables 2–4 of the Release and exposes the
 * computer-verification properties in Table 6.
 */

const Tt = 273.16;
const pt = 611.657;
const p0 = 101325;
const pi0 = p0 / pt;
const s0 = -0.332733756492168e4; // J kg^-1 K^-1, IAPWS-95 reference

const g0c = [
  -0.632020233335886e6,
   0.655022213658955,
  -0.189369929326131e-7,
   0.339746123271053e-14,
  -0.556464869058991e-21,
];

const t1 = { re: 0.368017112855051e-1, im: 0.510878114959572e-1 };
const r1 = { re: 0.447050716285388e2, im: 0.656876847463481e2 };
const t2 = { re: 0.337315741065416, im: 0.335449415919309 };
const r20 = { re: -0.725974574329220e2, im: -0.781008427112870e2 };
const r21 = { re: -0.557107698030123e-4, im: 0.464578634580806e-4 };
const r22 = { re: 0.234801409215913e-10, im: -0.285651142904972e-10 };

const C = (re, im=0) => ({re,im});
const add = (a,b) => C(a.re+b.re,a.im+b.im);
const sub = (a,b) => C(a.re-b.re,a.im-b.im);
const mul = (a,b) => C(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re);
const scale = (a,s) => C(a.re*s,a.im*s);
const div = (a,b) => { const d=b.re*b.re+b.im*b.im; return C((a.re*b.re+a.im*b.im)/d,(a.im*b.re-a.re*b.im)/d); };
const logc = a => C(Math.log(Math.hypot(a.re,a.im)), Math.atan2(a.im,a.re));

function F(t, tau) {
  const tm=sub(t,C(tau)), tp=add(t,C(tau));
  return sub(
    add(mul(tm,logc(tm)),mul(tp,logc(tp))),
    add(scale(mul(t,logc(t)),2),div(C(tau*tau),t))
  );
}
function FT(t,tau) {
  return add(sub(logc(add(t,C(tau))),logc(sub(t,C(tau)))),scale(div(C(tau),t),-2));
}
function FTT(t,tau) {
  return add(add(div(C(1),sub(t,C(tau))),div(C(1),add(t,C(tau)))),scale(div(C(1),t),-2));
}

function poly(coeff, x) { let y=0,p=1; for(const c of coeff){y+=c*p;p*=x;} return y; }
function poly1(coeff,x,den){ let y=0,p=1; for(let k=1;k<coeff.length;k++){y+=k*coeff[k]*p/den;p*=x;} return y; }
function poly2(coeff,x,den){ let y=0,p=1; for(let k=2;k<coeff.length;k++){y+=k*(k-1)*coeff[k]*p/(den*den);p*=x;} return y; }
function cpoly(c0,c1,c2,x){return add(add(c0,scale(c1,x)),scale(c2,x*x));}

export function iceIhGibbsDerivatives(T, pPa = p0) {
  if (!(T>0 && T<=Tt)) throw new RangeError('Ice Ih EOS requires 0 < T <= 273.16 K.');
  if (!(pPa>=0 && pPa<=210e6)) throw new RangeError('Ice Ih EOS release range is 0 <= p <= 210 MPa.');
  const tau=T/Tt, pi=pPa/pt, dpi=pi-pi0;
  const g0=poly(g0c,dpi), g0p=poly1(g0c,dpi,pt), g0pp=poly2(g0c,dpi,pt);
  const r2=cpoly(r20,r21,r22,dpi);
  const r2p=scale(add(r21,scale(r22,2*dpi)),1/pt);
  const r2pp=scale(r22,2/(pt*pt));
  const f1=F(t1,tau),f2=F(t2,tau);
  const ft1=FT(t1,tau),ft2=FT(t2,tau);
  const ftt1=FTT(t1,tau),ftt2=FTT(t2,tau);
  const g = g0 - s0*Tt*tau + Tt*(mul(r1,f1).re + mul(r2,f2).re);
  const gT = -s0 + mul(r1,ft1).re + mul(r2,ft2).re;
  const gp = g0p + Tt*mul(r2p,f2).re;
  const gTT = (mul(r1,ftt1).re + mul(r2,ftt2).re)/Tt;
  const gTp = mul(r2p,ft2).re;
  const gpp = g0pp + Tt*mul(r2pp,f2).re;
  return {T,p:pPa,g,gT,gp,gTT,gTp,gpp};
}

export function iceIhProperties(T,pPa=p0){
  const d=iceIhGibbsDerivatives(T,pPa);
  const v=d.gp, rho=1/v, s=-d.gT, cp=-T*d.gTT;
  const h=d.g-T*d.gT;
  const f=d.g-pPa*d.gp;
  const u=d.g-T*d.gT-pPa*d.gp;
  const alpha=d.gTp/d.gp;
  const beta=-d.gTp/d.gpp;
  const kappaT=-d.gpp/d.gp;
  const kappaS=(d.gTp*d.gTp-d.gTT*d.gpp)/(d.gp*d.gTT);
  return {...d,v,rho,s,cp,h,f,u,alpha,beta,kappaT,kappaS};
}

export const IAPWS_ICE_IH = Object.freeze({
  Tt,pt,p0,pi0,s0,g0c:[...g0c],t1:{...t1},r1:{...r1},t2:{...t2},r20:{...r20},r21:{...r21},r22:{...r22},
  reference:'IAPWS R10-06(2009), IAPWS-95 entropy reference',
  validity:'0 < T <= 273.16 K, 0 <= p <= 210 MPa within the Ice-Ih existence range',
});
