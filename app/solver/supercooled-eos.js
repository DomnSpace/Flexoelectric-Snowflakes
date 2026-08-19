/**
 * IAPWS G12-15 supercooled-water Gibbs-energy equation of state.
 * Reference: IAPWS, Guideline on Thermodynamic Properties of Supercooled Water (2015).
 *
 * Density is evaluated from the analytic pressure derivative. Temperature derivatives
 * are centered derivatives of the same fundamental Gibbs function, so h, s, cp,
 * alpha_p and kappa_T remain thermodynamically tied to one potential.
 */

const EOS = {omega0:0.5212269,L0:0.76317954,k0:0.072158686,k1:-0.31569232,k2:5.2992608,TLL:228.2,rho0:1081.6482,R:461.523087};
EOS.pi0=300e6/(EOS.rho0*EOS.R*EOS.TLL);
const C=[[-8.1570681381655,0,0,0],[1.2875032,0,1,0],[7.0901673598012,1,0,0],[-3.2779161e-2,-.2555,2.1051,-.0016],[7.3703949e-1,1.5762,1.1422,.6894],[-2.1628622e-1,1.64,.951,.013],[-5.1782479,3.6385,0,.0002],[4.2293517e-4,-.3828,3.6402,.0435],[2.3592109e-2,1.6219,2.076,.05],[4.3773754,4.3287,-.0016,.0004],[-2.996777e-3,3.4763,2.2769,.0528],[-9.6558018e-1,5.1556,.0008,.0147],[3.7595286,-.3593,.3706,.8584],[1.2632441,5.0361,-.3975,.9924],[2.8542697e-1,2.9786,2.973,1.0041],[-8.5994947e-1,6.2373,-.318,1.0961],[-3.2916153e-1,4.046,2.9805,1.0228],[9.0019616e-2,5.3558,2.9265,1.0303],[8.1149726e-2,9.0157,.4456,1.618],[-3.2788213,1.2194,.1298,.5213]];
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
function solveX(L,o){const c1=10/9*(Math.log(19)-L),c2=50/49*(Math.log(99)-L);let lo,hi;if(o<c1){lo=.049;hi=.5}else if(o<c2){lo=.0099;hi=.051}else{lo=.99*Math.exp(-50/49*L-o);hi=Math.min(1.1*Math.exp(-L-o),.0101)}lo=clamp(lo,1e-14,1-1e-14);hi=clamp(hi,lo+1e-14,1-1e-14);const f=x=>L+Math.log(x/(1-x))+o*(1-2*x);let fl=f(lo);if(fl*f(hi)>0){lo=1e-12;hi=.5;fl=f(lo)}for(let i=0;i<100;i++){const m=.5*(lo+hi),fm=f(m);if(fl*fm<=0)hi=m;else{lo=m;fl=fm}}return .5*(lo+hi)}
function core(T,p){const {omega0,L0,k0,k1,k2,TLL,rho0,R}=EOS,tau=T/TLL-1,pi=p/(rho0*R*TLL),tb=tau+1,pb=pi+EOS.pi0,K2=Math.sqrt(1+k2*k2),q=pi-k2*tau,A=1+k0*k2+k1*q,K1=Math.sqrt(A*A-4*k0*k1*k2*q),L=L0*K2/(2*k1*k2)*(1+k0*k2+k1*(pi+k2*tau)-K1),o=2+omega0*pi,x=solveX(L,o),phi=2*x-1,Lpi=L0*K2*(K1+k0*k2-k1*pi+k1*k2*tau-1)/(2*k2*K1);let pr=0,pp=0;for(const[c,a,b,d]of C){const e=Math.exp(-d*pb);pr+=c*tb**a*pb**b*e;pp+=c*tb**a*pb**(b-1)*(b-d*pb)*e}const mix=x*L+x*Math.log(x)+(1-x)*Math.log(1-x)+o*x*(1-x),psi=pr+tb*mix,B=omega0/2*(1-phi*phi)+Lpi*(phi+1),v=(tb/2*B+pp)/rho0;return{T,p,tau,pi,L,omega:o,x,phi,psi,g:psi*R*TLL,v,rho:1/v}}
export function supercooledWaterCore(T,pPa=101325){return core(T,pPa)}
export function homogeneousNucleationPressureLimit(T){const th=T/235.15;return(.1+228.27*(1-th**6.243)+15.724*(1-th**79.81))*1e6}
export function homogeneousNucleationTemperatureLimit(pPa=101325){const p=pPa/1e6;if(p>198.9)return 172.82+.03718*p+3.403e-5*p*p-1.573e-8*p**3;let lo=175,hi=240;for(let i=0;i<80;i++){const m=.5*(lo+hi);if(homogeneousNucleationPressureLimit(m)/1e6>p)lo=m;else hi=m}return .5*(lo+hi)}
export function supercooledWaterProperties(T,pPa=101325){const c=core(T,pPa),hT=.01,gpT=core(T+hT,pPa).g,gmT=core(T-hT,pPa).g,gT=(gpT-gmT)/(2*hT),gTT=(gpT-2*c.g+gmT)/(hT*hT),s=-gT,h=c.g-T*gT,cp=-T*gTT,vpT=core(T+hT,pPa).v,vmT=core(T-hT,pPa).v,alphaP=(vpT-vmT)/(2*hT*c.v),hp=Math.max(100,Math.abs(pPa)*1e-6),vpP=core(T,pPa+hp).v,vmP=core(T,Math.max(1,pPa-hp)).v,kappaT=-(vpP-vmP)/(2*hp*c.v),cv=cp-T*alphaP*alphaP/(c.rho*kappaT),speedSound=Math.pow(c.rho*kappaT*cv/cp,-.5),TH=homogeneousNucleationTemperatureLimit(pPa);return{...c,s,h,cp,cv,alphaP,kappaT,speedSound,homogeneousNucleationLimitK:TH,insideRecommendedMetastableRange:T>=TH&&T<=300&&pPa>=0&&pPa<=400e6}}
export const IAPWS_G12_15=Object.freeze({...EOS,coefficients:C.map(r=>[...r])});
