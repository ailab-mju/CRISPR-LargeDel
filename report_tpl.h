/* report_tpl.h - kmer_align HTML 리포트 템플릿.
   CSS 는 dataviz 검증 팔레트(라이트/다크 모두 통과)를 그대로 씁니다.
   JS 는 R (report.cpp 가 내보내는 JSON) 하나로 페이지 전체를 그리므로
   단일 실행 리포트와 여러 실행을 합친 리포트가 같은 템플릿을 씁니다.
   외부 라이브러리/CDN 없음 - 파일 하나로 어디서나 열립니다. */
#ifndef KA_REPORT_TPL_H
#define KA_REPORT_TPL_H
static const char *KA_CSS = R"KACSS(
*{box-sizing:border-box}
:root{
 color-scheme:light;
 --bg:#f7f7f5; --surface:#fcfcfb; --rule:#e6e5e0; --rule-2:#d2d1ca;
 --ink:#0b0b0b; --ink-2:#52514e; --ink-3:#85847e;
 --ins:#2a78d6; --del:#eb6834; --clip:#1baf7a; --sub:#eda100; --ldel:#e87ba4; --wt:#b9b8b1;
 --accent:#2a78d6;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
 color-scheme:dark;
 --bg:#111110; --surface:#1a1a19; --rule:#2e2e2b; --rule-2:#44443f;
 --ink:#fff; --ink-2:#c3c2b7; --ink-3:#8f8e86;
 --ins:#3987e5; --del:#d95926; --clip:#199e70; --sub:#c98500; --ldel:#d55181; --wt:#5a5a55;
 --accent:#3987e5;
}}
:root[data-theme="dark"]{
 color-scheme:dark;
 --bg:#111110; --surface:#1a1a19; --rule:#2e2e2b; --rule-2:#44443f;
 --ink:#fff; --ink-2:#c3c2b7; --ink-3:#8f8e86;
 --ins:#3987e5; --del:#d95926; --clip:#199e70; --sub:#c98500; --ldel:#d55181; --wt:#5a5a55;
 --accent:#3987e5;
}
html,body{margin:0;background:var(--bg);color:var(--ink)}
body{font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans KR",sans-serif;
 -webkit-font-smoothing:antialiased}
.wrap{max-width:1180px;margin:0 auto;padding:40px 16px 80px}
header{border-bottom:1px solid var(--rule-2);padding-bottom:22px;margin-bottom:8px}
h1{font-size:27px;line-height:1.25;margin:0 0 6px;letter-spacing:-.01em}
.sub{color:var(--ink-2);margin:0;font-size:14.5px}
.meta{display:flex;flex-wrap:wrap;gap:6px 18px;margin-top:14px;font:12.5px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--ink-3)}
.meta b{color:var(--ink-2);font-weight:500}
section{background:var(--surface);border:1px solid var(--rule);border-radius:14px;
 padding:24px 22px;margin-top:22px}
h2{font-size:18.5px;margin:0 0 4px;letter-spacing:-.005em}
.note{color:var(--ink-2);font-size:13.5px;margin:6px 0 0}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-top:16px}
.tile{border:1px solid var(--rule);border-radius:10px;padding:13px 14px;background:var(--bg)}
.tile .k{display:block;font-size:11.5px;color:var(--ink-3);text-transform:uppercase;letter-spacing:.04em}
.tile .v{display:block;font:600 22px/1.2 ui-monospace,SFMono-Regular,Menlo,monospace;margin-top:5px}
.tile .n{display:block;font-size:12px;color:var(--ink-3);margin-top:3px}
.legend{display:flex;flex-wrap:wrap;gap:8px 16px;margin:14px 0 4px;font-size:13px;color:var(--ink-2)}
.legend span{display:inline-flex;align-items:center;gap:6px}
.sw{width:11px;height:11px;border-radius:3px;display:inline-block;flex:none}
.chart{position:relative;margin-top:14px}
.chart svg{display:block;width:100%;height:auto;overflow:visible}
.tip{position:absolute;pointer-events:none;opacity:0;transition:opacity .09s;
 background:var(--surface);border:1px solid var(--rule-2);border-radius:9px;padding:9px 11px;
 font-size:12.5px;box-shadow:0 6px 22px rgba(0,0,0,.14);min-width:130px;z-index:5;color:var(--ink)}
.tip b{font:600 13px ui-monospace,SFMono-Regular,Menlo,monospace}
.tip .r{display:flex;justify-content:space-between;gap:14px;margin-top:4px;color:var(--ink-2)}
.tabs{display:flex;gap:6px;margin-top:16px;flex-wrap:wrap}
.tabs button{font:inherit;font-size:13px;padding:5px 12px;border-radius:99px;cursor:pointer;
 border:1px solid var(--rule-2);background:transparent;color:var(--ink-2)}
.tabs button[aria-selected="true"]{background:var(--accent);border-color:var(--accent);color:#fff}
table{width:100%;border-collapse:collapse;font-size:13.5px;margin-top:14px}
th,td{text-align:left;padding:7px 9px;border-bottom:1px solid var(--rule)}
th{font-size:11.5px;text-transform:uppercase;letter-spacing:.04em;color:var(--ink-3);font-weight:600}
td.n,th.n{text-align:right;font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.seq{font:12px ui-monospace,SFMono-Regular,Menlo,monospace;word-break:break-all}
.flag{display:inline-block;font-size:11px;padding:1px 7px;border-radius:99px;
 background:color-mix(in srgb,var(--del) 16%,transparent);color:var(--del);font-weight:600}
.ok{background:color-mix(in srgb,var(--clip) 16%,transparent);color:var(--clip)}
.thin{background:color-mix(in srgb,var(--ink-3) 14%,transparent);color:var(--ink-2)}
.tblwrap{overflow-x:auto}
details{margin-top:14px}
summary{cursor:pointer;color:var(--ink-2);font-size:13.5px}
.warn{border-left:3px solid var(--sub);padding:10px 0 10px 14px;margin-top:16px;
 color:var(--ink-2);font-size:13.5px;background:color-mix(in srgb,var(--sub) 7%,transparent);
 border-radius:0 8px 8px 0}
footer{margin-top:34px;color:var(--ink-3);font-size:12.5px;text-align:center}
@media(max-width:640px){.wrap{padding:26px 16px 60px}h1{font-size:22px}}
)KACSS";
static const char *KA_JS = R"KAJS("use strict";
const F=n=>Number(n||0).toLocaleString();
const css=v=>getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const esc=x=>String(x).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const CAT=[["wt","야생형","--wt"],["small_insertion","작은 삽입","--ins"],["small_deletion","작은 결실","--del"],
           ["terminal_insertion_clipped","말단 삽입","--clip"],["substitution","치환","--sub"],["large_deletion","대형 결실","--ldel"]];
const STACK=["small_insertion","small_deletion","terminal_insertion_clipped","substitution","large_deletion","wt"];
function sumMaps(list,key){const o={};for(const s of list){const m=s[key]||{};for(const k in m){o[k]=(o[k]||0)+m[k];}}return o;}
function pct(s){const t=(s.treated||{}).total_analyzed||0,o={};for(const [k] of CAT){o[k]=t?100*((s.treated||{})[k]||0)/t:0;}return o;}
const baseOf=l=>String(l).replace(/[_-](ON|OFF|on|off)$/,"");
function svg(w,h,inner,label){return `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="${label}">${inner}</svg>`;}
function hookTip(box,tip,hit){
  const s=box.querySelector("svg");
  s.addEventListener("pointermove",e=>{
    const r=s.getBoundingClientRect();
    const res=hit((e.clientX-r.left)/r.width,(e.clientY-r.top)/r.height);
    if(!res){tip.style.opacity=0;return;}
    tip.innerHTML=res.html;tip.style.opacity=1;
    const bx=box.getBoundingClientRect();
    let x=(e.clientX-bx.left)+14,y=(e.clientY-bx.top)-10;
    x=Math.min(Math.max(x,4),bx.width-tip.offsetWidth-4);
    y=Math.min(Math.max(y,4),bx.height-tip.offsetHeight-4);
    tip.style.left=x+"px";tip.style.top=y+"px";
  });
  s.addEventListener("pointerleave",()=>tip.style.opacity=0);
}
function yAxis(max,logY,py,L,W,R){
  let g="";
  for(let i=0;i<=4;i++){
    const v=logY?Math.round(Math.pow(10,i*Math.log10(max+1)/4)-1):max*i/4;
    g+=`<line x1="${L}" y1="${py(v)}" x2="${W-R}" y2="${py(v)}" stroke="${css('--rule')}"/>`;
    g+=`<text x="${L-8}" y="${py(v)+4}" text-anchor="end" font-size="11.5" fill="${css('--ink-3')}" font-family="ui-monospace,monospace">${F(Math.round(v))}</text>`;
  }
  return g;
}

/* ---------- 인델 크기 스펙트럼 ---------- */
function spectrum(box,tip,data,lo,hi,logY){
  box.querySelectorAll("svg").forEach(e=>e.remove());
  const W=1000,H=300,L=58,R=16,T=18,B=52;
  const bins=[];for(let k=lo;k<=hi;k++){if(k!==0)bins.push(k);}
  const val=k=>k<0?(data.del[-k]||0):(data.ins[k]||0);
  const max=Math.max(1,...bins.map(val));
  const sc=v=>logY?Math.log10(v+1)/Math.log10(max+1):v/max;
  const bw=(W-L-R)/bins.length,py=v=>T+(1-sc(v))*(H-T-B);
  let g=yAxis(max,logY,py,L,W,R),bars="";
  bins.forEach((k,i)=>{const v=val(k);if(v<=0)return;const x=L+i*bw,y=py(v);
    bars+=`<rect x="${(x+0.6).toFixed(1)}" y="${y.toFixed(1)}" width="${Math.max(bw-1.2,0.8).toFixed(1)}" height="${Math.max(py(0)-y,1).toFixed(1)}" rx="${Math.min(2,bw/3).toFixed(1)}" fill="${css(k<0?'--del':'--ins')}"/>`;});
  const step=bins.length>40?10:(bins.length>16?5:1);
  bins.forEach((k,i)=>{if(Math.abs(k)%step)return;
    g+=`<text x="${(L+i*bw+bw/2).toFixed(1)}" y="${H-B+17}" text-anchor="middle" font-size="11" fill="${css('--ink-3')}" font-family="ui-monospace,monospace">${k>0?'+':''}${k}</text>`;});
  const zi=bins.findIndex(k=>k>0),zx=L+(zi<0?bins.length:zi)*bw;
  g+=`<line x1="${zx}" y1="${T}" x2="${zx}" y2="${H-B}" stroke="${css('--rule-2')}" stroke-dasharray="3 3"/>`;
  g+=`<text x="${(L+zx)/2}" y="${H-B+36}" text-anchor="middle" font-size="12" fill="${css('--del')}">← 결실 (bp)</text>`;
  g+=`<text x="${(zx+W-R)/2}" y="${H-B+36}" text-anchor="middle" font-size="12" fill="${css('--ins')}">삽입 (bp) →</text>`;
  box.insertAdjacentHTML("afterbegin",svg(W,H,g+bars,"인델 크기 스펙트럼"));
  hookTip(box,tip,fx=>{const i=Math.floor((fx*W-L)/bw);if(i<0||i>=bins.length)return null;
    const k=bins[i],v=val(k);if(!v)return null;
    return {html:`<b>${k>0?'+':''}${k} bp ${k<0?'결실':'삽입'}</b><div class="r"><span>read</span><b>${F(v)}</b></div>`};});
}

/* ---------- 분류 누적 막대 ---------- */
function stacks(box,tip,rows){
  const W=1000,rh=34,L=110,R=56,T=8,H=T+rows.length*rh+18;
  const tok=Object.fromEntries(CAT.map(([k,,t])=>[k,t])),lab=Object.fromEntries(CAT.map(([k,l])=>[k,l]));
  let g="";
  rows.forEach((r,i)=>{const y=T+i*rh+5,bh=rh-14;let acc=0;
    STACK.forEach(k=>{const p=r.pct[k]||0;if(p<=0)return;const x=L+(W-L-R)*acc/100,w=(W-L-R)*p/100;
      g+=`<rect x="${x.toFixed(1)}" y="${y}" width="${Math.max(w-2,0.8).toFixed(1)}" height="${bh}" rx="3" fill="${css(tok[k])}"/>`;acc+=p;});
    g+=`<text x="${L-10}" y="${y+bh/2+4}" text-anchor="end" font-size="12.5" fill="${css('--ink-2')}" font-family="ui-monospace,monospace">${esc(r.label)}</text>`;
    g+=`<text x="${W-R+8}" y="${y+bh/2+4}" font-size="12" fill="${css('--ink-3')}" font-family="ui-monospace,monospace">${F(r.total)}</text>`;});
  box.insertAdjacentHTML("afterbegin",svg(W,H,g,"샘플별 분류 구성"));
  hookTip(box,tip,(fx,fy)=>{const i=Math.floor((fy*H-T)/rh);if(i<0||i>=rows.length)return null;const r=rows[i];
    return {html:`<b>${esc(r.label)}</b>`+STACK.map(k=>`<div class="r"><span><i class="sw" style="background:${css(tok[k])}"></i> ${lab[k]}</span><b>${(r.pct[k]||0).toFixed(2)}%</b></div>`).join("")};});
}

/* ---------- 구간 막대 (대형 결실 / 삽입 길이) ---------- */
function binned(box,tip,series,BINS,LAB,opt){
  const W=1000,H=290,L=58,R=16,T=18,B=54;
  const gw=(W-L-R)/BINS.length,n=series.length,bw=(gw-12)/n;
  const agg=series.map(s=>BINS.map(([a,b])=>{let t=0;for(const k in s.data){const v=+k;if(v>=a&&v<=b)t+=s.data[k];}return t;}));
  const max=Math.max(1,...agg.flat()),logY=!!opt.logY;
  const sc=v=>logY?Math.log10(v+1)/Math.log10(max+1):v/max,py=v=>T+(1-sc(v))*(H-T-B);
  let g=yAxis(max,logY,py,L,W,R),bars="";
  BINS.forEach((b,bi)=>{series.forEach((s,si)=>{const v=agg[si][bi];if(!v)return;const x=L+bi*gw+6+si*bw,y=py(v);
      bars+=`<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${Math.max(bw-2,1).toFixed(1)}" height="${(py(0)-y).toFixed(1)}" rx="3" fill="${css(s.tok)}"/>`;});
    g+=`<text x="${(L+bi*gw+gw/2).toFixed(1)}" y="${H-B+17}" text-anchor="middle" font-size="11.5" fill="${css('--ink-3')}" font-family="ui-monospace,monospace">${LAB[bi]}</text>`;});
  if(opt.mark!=null){const x=L+opt.mark*gw;
    g+=`<line x1="${x}" y1="${T}" x2="${x}" y2="${H-B}" stroke="${css(opt.markTok||'--rule-2')}" stroke-dasharray="4 3"/>`;
    g+=`<text x="${x+6}" y="${T+13}" font-size="11.5" fill="${css(opt.markTok||'--ink-2')}">${opt.markText}</text>`;}
  g+=`<text x="${(L+W-R)/2}" y="${H-B+38}" text-anchor="middle" font-size="12" fill="${css('--ink-2')}">${opt.xlabel}</text>`;
  box.insertAdjacentHTML("afterbegin",svg(W,H,g+bars,opt.xlabel));
  hookTip(box,tip,fx=>{const bi=Math.floor((fx*W-L)/gw);if(bi<0||bi>=BINS.length)return null;
    return {html:`<b>${LAB[bi]} bp</b>`+series.map((s,si)=>`<div class="r"><span><i class="sw" style="background:${css(s.tok)}"></i> ${s.label}</span><b>${F(agg[si][bi])}</b></div>`).join("")};});
}

/* ---------- 페이지 ---------- */
function tile(k,v,n){return `<div class="tile"><span class="k">${k}</span><span class="v">${v}</span><span class="n">${n}</span></div>`;}
function legend(items){return `<div class="legend">`+items.map(([t,l])=>`<span><i class="sw" style="background:var(${t})"></i> ${l}</span>`).join("")+`</div>`;}

function build(){
  const S=R.samples,C=R.compare||[],p0=S[0].params||{};
  const ins=sumMaps(S,"mut_ins"),del=sumMaps(S,"mut_del"),ld=sumMaps(S,"ldel"),ldc=sumMaps(S,"ldel_ctr");
  const rows=S.map(s=>({label:s.label,total:(s.treated||{}).total_analyzed||0,pct:pct(s)}));
  const totA=S.reduce((a,s)=>a+((s.treated||{}).aligned_pairs_all_positions||0),0);
  const totN=rows.reduce((a,r)=>a+r.total,0);
  const edited=rows.reduce((a,r)=>a+r.total*(100-r.pct.wt)/100,0);
  const ik=Object.keys(ins).map(Number),maxI=ik.length?Math.max(...ik):0;
  const big=ik.filter(k=>k>=8).reduce((a,k)=>a+ins[k],0);

  /* 8 bp 이상 삽입 서열 — 원본 read 스캔 결과를 샘플에 걸쳐 합산 */
  const calls={};S.forEach(s=>(s.seqs||[]).forEach(([q,n])=>calls[q]=(calls[q]||0)+n));
  const top=Object.entries(calls).sort((a,b)=>b[1]-a[1]||(a[0]<b[0]?-1:1)).slice(0,12);
  const anyScan=S.some(s=>s.scan);
  let seqRows="";
  for(const [q,n] of top){
    let mh=0,mr=0,ch=0,cr=0,ns=0;
    for(const s of S){const sc=s.scan;if(!sc||!(q in sc.mut))continue;ns++;
      mh+=sc.mut[q];mr+=sc.mut_reads;ch+=sc.ctr[q]||0;cr+=sc.ctr_reads;}
    if(ns){const mp=1e6*mh/Math.max(mr,1),cp=1e6*ch/Math.max(cr,1),art=cp>0&&cp>=0.5*mp,thin=!art&&mh<5;
      /* 원본에서 5번도 안 나오면 어느 쪽으로도 단정하지 않습니다 */
      const vt=art?'대조군에도 있음 · 아티팩트':(thin?`근거 부족 · 원본 ${mh}회`:'처리군 특이'),vc=art?'':(thin?'thin':'ok');
      seqRows+=`<tr><td class="n">${q.length}</td><td class="n">${F(n)}</td><td class="n">${F(Math.round(mp))}</td><td class="n">${F(Math.round(cp))}</td>`+
        `<td><span class="flag ${vc}">${vt}</span>${S.length>1?` <span style="color:var(--ink-3);font-size:11.5px">${ns}/${S.length}</span>`:''}</td><td class="seq">${esc(q)}</td></tr>`;}
    else seqRows+=`<tr><td class="n">${q.length}</td><td class="n">${F(n)}</td><td class="n">–</td><td class="n">–</td><td><span class="flag">미확인</span></td><td class="seq">${esc(q)}</td></tr>`;
  }
  if(!seqRows)seqRows=`<tr><td colspan="6" style="color:var(--ink-3)">8 bp 이상 삽입 없음</td></tr>`;
  let scanNote="";
  if(anyScan){const mr=S.reduce((a,s)=>a+(s.scan?s.scan.mut_reads:0),0),cr=S.reduce((a,s)=>a+(s.scan?s.scan.ctr_reads:0),0);
    scanNote=`처리군 원본 ${F(mr)} read, 대조군 원본 ${F(cr)} read 를 직접 세어 백만 read 당 빈도로 비교했습니다.`;}
  else scanNote=`원본 FASTQ 를 읽지 못해 판정하지 않았습니다.`;

  const meta=[["gene",p0.gene],["targetst",p0.targetst],["cle",p0.cle],["절단",p0.cut],["K",p0.K],["dist_cut",p0.dist_cut],
              ["anchor_min",p0.anchor_min],["cut_between",p0.cut_between],["ins_fix",p0.ins_fix],["ins_max",p0.ins_max],["ref",p0.ref_len?p0.ref_len+" bp":undefined]]
    .filter(([,v])=>v!==undefined&&v!=="").map(([k,v])=>`<span><b>${k}</b> ${esc(v)}</span>`).join(" ");

  const clsRows=rows.map(r=>`<tr><td>${esc(r.label)}</td>`+CAT.map(([k])=>`<td class="n">${r.pct[k].toFixed(2)}%</td>`).join("")+`<td class="n">${F(r.total)}</td></tr>`).join("");

  let compHTML="";
  if(C.length){
    const cins=sumMaps(C,"mut_ins"),ck=Object.keys(cins).map(Number);
    const bigB=ck.filter(k=>k>=8).reduce((a,k)=>a+cins[k],0),mxB=ck.length?Math.max(...ck):0;
    const byBase={};C.forEach((s,i)=>byBase[baseOf(s.label)]=s);
    const trs=S.map((s,i)=>{const o=byBase[baseOf(s.label)]||(C.length===S.length?C[i]:null);if(!o)return "";
      const a=pct(o),b=pct(s);
      return `<tr><td>${esc(baseOf(s.label))}</td>`+CAT.map(([k])=>`<td class="n">${a[k].toFixed(2)} → ${b[k].toFixed(2)}</td>`).join("")+`</tr>`;}).join("");
    compHTML=`<section><h2>삽입 정렬 수정 전후</h2>
      <p class="note">같은 데이터를 <code>ins_fix=0</code> 과 <code>ins_fix=1</code> 로 돌린 결과입니다.</p>
      <div class="tiles">${tile("최대 삽입 길이",`${mxB} → ${maxI} bp`,"CIGAR 의 I")}${tile("8 bp 이상 삽입",`${F(bigB)} → ${F(big)}`,"read 수")}</div>
      ${legend([["--wt","수정 전"],["--ins","수정 후"]])}
      <div class="chart" id="b-ins"><div class="tip" id="t-ins"></div></div>
      <details open><summary>분류 비율 전 → 후 (%)</summary><div class="tblwrap"><table>
      <thead><tr><th>샘플</th>${CAT.map(([,l])=>`<th class="n">${l}</th>`).join("")}</tr></thead><tbody>${trs}</tbody></table></div></details></section>`;
    R._cins=cins;
  }

  document.getElementById("app").innerHTML=`
  <header><h1>${esc(R.title)}</h1>
    <p class="sub">kmer_align ${esc(R.version)} · 샘플 ${S.length}개 · 생성 ${esc(R.date)}</p>
    <div class="meta">${meta}</div></header>
  <section><h2>요약</h2><div class="tiles">
    ${tile("정렬된 pair",F(totA),"전체 위치")}${tile("분류 대상",F(totN),"절단 부위를 지나는 pair")}
    ${tile("편집된 비율",(100*edited/Math.max(totN,1)).toFixed(1)+"%","야생형이 아닌 것")}
    ${tile("최대 삽입 길이",maxI+" bp","CIGAR 의 I")}${tile("8 bp 이상 삽입",F(big),"read 수")}</div></section>
  <section><h2>인델 크기 스펙트럼</h2>
    <p class="note">절단 부위 ±13 bp 안의 인델을 길이별로 센 것입니다. 왼쪽이 결실, 오른쪽이 삽입.</p>
    ${legend([["--del","결실"],["--ins","삽입"]])}
    <div class="chart" id="b-spec"><div class="tip" id="t-spec"></div></div>
    <div class="tabs" role="tablist">
      <button role="tab" aria-selected="true" data-lo="-30" data-hi="30">±30 bp</button>
      <button role="tab" aria-selected="false" data-lo="-10" data-hi="10">±10 bp</button>
      <button role="tab" aria-selected="false" data-lo="-49" data-hi="100">전체</button>
      <button id="logbtn" aria-pressed="false" style="margin-left:10px">y축 로그</button></div></section>
  <section><h2>${S.length>1?"샘플별 분류 구성":"분류 구성"}</h2>
    ${legend(STACK.map(k=>[CAT.find(c=>c[0]===k)[2],CAT.find(c=>c[0]===k)[1]]))}
    <div class="chart" id="b-stack"><div class="tip" id="t-stack"></div></div>
    <details><summary>표로 보기</summary><div class="tblwrap"><table>
      <thead><tr><th>샘플</th>${CAT.map(([,l])=>`<th class="n">${l}</th>`).join("")}<th class="n">합계</th></tr></thead><tbody>${clsRows}</tbody></table></div></details></section>
  <section><h2>대형 결실 크기 분포</h2>
    <p class="note">split read 의 CIGAR 에 적힌 결실 길이입니다. 대조군을 나란히 두어 아티팩트 수준을 함께 봅니다.</p>
    ${legend([["--ldel","처리군"],["--wt","대조군"]])}
    <div class="chart" id="b-ldel"><div class="tip" id="t-ldel"></div></div>
    <div class="warn"><b>2 kb 이상 구간 주의.</b> only/common 분할은 입력에 민감해 단일 실행으로 재현되지 않습니다.
      입력의 1%만 달라져도 최대 30 pp 움직입니다. 이 구간을 보고할 때는 반복 표본으로 변동폭을 함께 제시하십시오.</div></section>
  <section><h2>8 bp 이상 삽입 서열</h2>
    <p class="note">대조군에서 같은 서열이 비슷하게 나오면 CRISPR 산물이 아니라 공동 증폭된 배경 서열입니다.
      파이프라인의 only/common 대조군 필터는 <b>대형 결실에만</b> 적용되고 삽입에는 적용되지 않으므로 여기서 확인합니다.</p>
    <p class="note">${scanNote}</p>
    <div class="tblwrap"><table><thead><tr><th class="n">길이</th><th class="n">삽입 호출</th><th class="n">처리군 원본 /M</th><th class="n">대조군 원본 /M</th><th>판정</th><th>서열</th></tr></thead>
      <tbody>${seqRows}</tbody></table></div></section>
  ${compHTML}
  <footer>kmer_align · CRISPR-LargeDel</footer>`;

  stacks(document.getElementById("b-stack"),document.getElementById("t-stack"),rows);
  const LDB=[[50,99],[100,249],[250,499],[500,999],[1000,1999],[2000,3999],[4000,1e12]];
  binned(document.getElementById("b-ldel"),document.getElementById("t-ldel"),
    [{label:"처리군",tok:"--ldel",data:ld},{label:"대조군",tok:"--wt",data:ldc}],LDB,
    ["50–99","100–249","250–499","500–999","1k–2k","2k–4k","4k+"],
    {mark:5,markTok:"--sub",markText:"2 kb 이상: 재현성 주의",xlabel:"결실 길이 (bp)"});
  if(C.length){
    binned(document.getElementById("b-ins"),document.getElementById("t-ins"),
      [{label:"수정 전",tok:"--wt",data:R._cins},{label:"수정 후",tok:"--ins",data:ins}],
      [[1,1],[2,3],[4,7],[8,15],[16,30],[31,60],[61,1e12]],["1","2–3","4–7","8–15","16–30","31–60","61+"],
      {logY:true,mark:3,markText:"← 수정 전에는 크기를 붙이지 못하던 구간",xlabel:"삽입 길이 (bp) · y축 로그"});
  }
  const sb=document.getElementById("b-spec"),st=document.getElementById("t-spec"),data={ins:ins,del:del};
  let RANGE=[-30,30],LOGY=false;const draw=()=>spectrum(sb,st,data,RANGE[0],RANGE[1],LOGY);draw();
  document.querySelectorAll('.tabs button[data-lo]').forEach(b=>b.addEventListener('click',()=>{
    document.querySelectorAll('.tabs button[data-lo]').forEach(x=>x.setAttribute('aria-selected','false'));
    b.setAttribute('aria-selected','true');RANGE=[+b.dataset.lo,+b.dataset.hi];draw();}));
  const lb=document.getElementById('logbtn');
  lb.addEventListener('click',()=>{LOGY=!LOGY;lb.setAttribute('aria-pressed',LOGY);lb.setAttribute('aria-selected',LOGY);draw();});
}
build();
)KAJS";
#endif
