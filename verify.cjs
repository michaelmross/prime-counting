'use strict';
const assert=require('node:assert/strict');
const oldShort=require('../research_cycle_1029/results.json').results;
const oldLong=require('../research_cycle_1031/results.json').rows;
const roots=[9239,9767,100151,300149,325019,350087,375017,400031];
function powFloor(P,a,b){
 const x=BigInt(P)**BigInt(a),k=BigInt(b);
 let u=BigInt(Math.floor(P**(a/b)));
 while((u+1n)**k<=x)++u;
 while(u**k>x)--u;
 return Number(u);
}
function sieve(n){
 const f=new Uint8Array(n+1).fill(1);f[0]=f[1]=0;
 for(let p=2;p*p<=n;++p)if(f[p])for(let i=p*p;i<=n;i+=p)f[i]=0;
 const ps=[];for(let i=2;i<=n;++i)if(f[i])ps.push(i);
 return {f,ps};
}
const cap=Math.max(...roots.map(P=>Math.floor(((P+2)**2+powFloor(P,5,3))/(powFloor(P,19,25)+1))));
const {f,ps}=sieve(cap);
const sums=new Float64Array(ps.length+1);
for(let i=0;i<ps.length;++i)sums[i+1]=sums[i]+ps[i];
const round=(x,k=8)=>+x.toFixed(k);
function upper(x){
 let lo=0,hi=ps.length;
 while(lo<hi){const m=(lo+hi)>>>1;if(ps[m]<=x)lo=m+1;else hi=m;}
 return lo;
}
function symbol(a,p){
 let r=1n,b=BigInt(a),m=BigInt(p),k=(m-1n)/2n;
 while(k){if(k&1n)r=r*b%m;b=b*b%m;k>>=1n;}
 return r===m-1n?-1:Number(r);
}
let checked=0;
function calc(P,a,b,variance,checkSymbols){
 const X=P*P,B=(P+2)**2,L=B-X,T=powFloor(P,a,b),D=powFloor(P,19,25);
 assert(T>L && f[P] && f[P+2]);
 const events=variance?new Int8Array(T):null;
 let endpoint=0,numer=0,diagNumerator=0n;
 for(let i=upper(D);i<ps.length&&4*ps[i]<3*P;++i){
  const r=ps[i];
  const j0=upper(Math.floor(X/r)),j1=upper(Math.floor(B/r));
  const j2=upper(Math.floor((X+T)/r)),j3=upper(Math.floor((B+T)/r));
  const left=r*(sums[j1]-sums[j0])-X*(j1-j0);
  const middle=L*(j2-j1);
  const right=(B+T)*(j3-j2)-r*(sums[j3]-sums[j2]);
  const cur=left+middle+right;
  assert([left,middle,right,cur].every(Number.isSafeInteger));
  assert(cur>=0 && (j0===j3 || ps[j0]>P));
  numer+=cur;
  endpoint+=j1-j0;
  if(checkSymbols)for(let j=j0;j<j1;++j){
   const s=ps[j],h=r*s-X;
   assert(h>=1&&h<=L);
   const a1=(r-h%r)%r,a2=(s-h%s)%s;
   assert(a1===X%r && a2===X%s && a1!==0 && a2!==0);
   assert.equal(symbol(a1,r),1);assert.equal(symbol(a2,s),1);
   checked+=2;
  }
  if(events){
   let within=0,kernelSum=0;
   for(let j=j0;j<j3;++j){
    const u=r*ps[j]-X,lo=Math.max(0,u-L),hi=Math.min(T,u);
    if(hi<=lo)continue;
    kernelSum+=hi-lo;within+=hi-lo;
    if(lo>0)++events[lo];
    if(hi<T)--events[hi];
    for(let k=j+1;k<j3;++k){
     const v=r*ps[k]-X,start=Math.max(0,v-L);
     if(start>=hi)break;
     within+=2*(hi-start);
    }
   }
   assert.equal(kernelSum,cur);
   assert(Number.isSafeInteger(within));
   diagNumerator+=BigInt(T)*BigInt(within)-BigInt(cur)**2n;
  }
 }
 assert(Number.isSafeInteger(numer));
 const mean=numer/T;
 const result={T,endpoint,meanNumerator:numer,mean:round(mean),residual:round(endpoint-mean)};
 if(events){
  const center=Math.floor(mean);
  let y=endpoint,sumsq=0,windowSum=0;
  for(let t=0;t<T;++t){y+=events[t];assert(y>=0);windowSum+=y;sumsq+=(y-center)**2;}
  assert.equal(windowSum,numer);assert(Number.isSafeInteger(sumsq));
  const varianceNumerator=BigInt(T)*BigInt(sumsq)-BigInt(numer-T*center)**2n;
  const crossNumerator=varianceNumerator-diagNumerator;
  const varTotal=Number(varianceNumerator)/T**2;
  Object.assign(result,{variance:round(varTotal),standardDeviation:round(Math.sqrt(varTotal)),diagonalVariance:round(Number(diagNumerator)/T**2),crossCovarianceSum:round(Number(crossNumerator)/T**2),crossCovarianceNumerator:String(crossNumerator),center,centeredSquareNumerator:sumsq});
 }
 return result;
}
const t0=Date.now();
const rows=roots.map((P,i)=>{
 const long=calc(P,5,3,i<2,true),short=calc(P,41,30,true,false);
 assert.equal(long.endpoint,short.endpoint);
 const a=oldShort.find(x=>x.P===P),b=oldLong.find(x=>x.P===P);
 if(a){assert.equal(short.endpoint,a.semiprimeEndpoint[1]);assert.equal(short.meanNumerator,a.semiprimeShiftNumerators[1]);}
 if(b)assert(Math.abs(long.mean-b.averageRemoval[1])<0.000001);
 return {P,long,short};
});
const witness={P:71,long:calc(71,5,3,true,true)};
assert(BigInt(witness.long.crossCovarianceNumerator)>0n);
console.log(JSON.stringify({checked,rows,positiveCrossCovarianceWitness:witness},null,2));
console.error(`PASS: ${rows.length} roots, ${checked} factor-symbol checks, ${(Date.now()-t0)/1000}s`);
