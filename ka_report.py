#!/usr/bin/env python3
"""
ka_report.py - kmer_align 결과 폴더들을 하나의 자기완결 HTML 리포트로 묶습니다.

  python3 ka_report.py -o report.html result/sg1_1_ON result/sg1_2_ON ...
  python3 ka_report.py -o report.html --pair result/sg1_1_OFF result/sg1_1_ON

외부 라이브러리/CDN 없이 동작하며, 결과 HTML 은 파일 하나로 어디서나 열립니다.
"""
import argparse, csv, glob, json, os, re, sys, collections, html, datetime

CAT = [("wt","야생형","wt"),
       ("small_insertion","작은 삽입","ins"),
       ("small_deletion","작은 결실","del"),
       ("terminal_insertion_clipped","말단 삽입","clip"),
       ("substitution","치환","sub"),
       ("large_deletion","대형 결실","ldel")]

def find_prefix(d):
    m = glob.glob(os.path.join(d, "*_mutation_summary.tsv"))
    if not m: return None
    return os.path.basename(m[0])[:-len("_mutation_summary.tsv")]

def read_summary(d, pre):
    out = {"control": {}, "treated": {}}
    with open(os.path.join(d, pre + "_mutation_summary.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try: v = int(r["read_pairs"])
            except ValueError: continue
            out[r["sample"]][r["pattern"]] = v
    return out

def cigar_ops(c):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", c)]

def scan_sam(path, want_seq=False):
    """CIGAR 의 I / D 길이와, 원하면 삽입 서열까지."""
    ins, dele, seqs = collections.Counter(), collections.Counter(), collections.Counter()
    if not os.path.exists(path): return ins, dele, seqs
    with open(path) as f:
        for line in f:
            if line.startswith("@"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10: continue
            ops = cigar_ops(p[5])
            i_tot = sum(n for n, o in ops if o == "I")
            d_max = max([n for n, o in ops if o == "D"] or [0])
            if i_tot: ins[i_tot] += 1
            if d_max: dele[d_max] += 1
            if want_seq and i_tot >= 8:
                q = 0
                for n, o in ops:
                    if o == "I":
                        seqs[p[9][q:q+n]] += 1
                        break
                    if o in "MIS=X": q += n
    return ins, dele, seqs

def collect(d):
    pre = find_prefix(d)
    if pre is None: return None
    s = read_summary(d, pre)
    mi, md, mseq = collections.Counter(), collections.Counter(), collections.Counter()
    for suf in ("site_mut_sm", "site_mut_clip"):
        a, b, c = scan_sam(os.path.join(d, f"{pre}_{suf}.sam"), want_seq=True)
        mi += a; md += b; mseq += c
    ci, cd, cseq = collections.Counter(), collections.Counter(), collections.Counter()
    for suf in ("site_ctr_sm", "site_ctr_clip"):
        a, b, c = scan_sam(os.path.join(d, f"{pre}_{suf}.sam"), want_seq=True)
        ci += a; cd += b; cseq += c
    ld = collections.Counter()
    for suf in ("mut_only", "mut_common"):
        _, b, _ = scan_sam(os.path.join(d, f"{pre}_{suf}.sam"))
        ld += b
    ldc = collections.Counter()
    for suf in ("ctr_only", "ctr_common"):
        _, b, _ = scan_sam(os.path.join(d, f"{pre}_{suf}.sam"))
        ldc += b
    log = os.path.join(os.path.dirname(d.rstrip("/")), "log_%s.txt" % os.path.basename(d.rstrip("/")))
    params = {}
    if os.path.exists(log):
        vals = [l.strip()[3:] for l in open(log, errors="replace") if l.startswith("---")]
        keys = ["ref","ctr_1","ctr_2","mut_1","mut_2","targetst","targeted","cle","gene",
                "sample","K","outdir","dist_cut","full_sam","out_fmt","min_gap","anchor_min",
                "cut_between","ins_fix","ins_max"]
        params = dict(zip(keys, vals))
    return {"dir": d, "name": pre, "summary": s, "params": params,
            "mut_ins": dict(mi), "mut_del": dict(md), "ctr_ins": dict(ci), "ctr_del": dict(cd),
            "mut_ins_seq": mseq.most_common(40), "ctr_ins_seq": dict(cseq),
            "ldel": dict(ld), "ldel_ctr": dict(ldc)}

# ------------------------------------------------------------------ HTML
CSS = """
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
.tblwrap{overflow-x:auto}
details{margin-top:14px}
summary{cursor:pointer;color:var(--ink-2);font-size:13.5px}
.warn{border-left:3px solid var(--sub);padding:10px 0 10px 14px;margin-top:16px;
 color:var(--ink-2);font-size:13.5px;background:color-mix(in srgb,var(--sub) 7%,transparent);
 border-radius:0 8px 8px 0}
footer{margin-top:34px;color:var(--ink-3);font-size:12.5px;text-align:center}
@media(max-width:640px){.wrap{padding:26px 16px 60px}h1{font-size:22px}}
"""

JS = r"""
const F=n=>n.toLocaleString();
const css=v=>getComputedStyle(document.documentElement).getPropertyValue(v).trim();
function svg(w,h,inner,label){
  return `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="${label}">${inner}</svg>`;
}
function hookTip(box,tip,hit){
  const s=box.querySelector("svg");
  s.addEventListener("pointermove",e=>{
    const r=s.getBoundingClientRect();
    const res=hit((e.clientX-r.left)/r.width,(e.clientY-r.top)/r.height);
    if(!res){tip.style.opacity=0;return;}
    tip.innerHTML=res.html; tip.style.opacity=1;
    const bx=box.getBoundingClientRect();
    let x=(e.clientX-bx.left)+14, y=(e.clientY-bx.top)-10;
    x=Math.min(Math.max(x,4),bx.width-tip.offsetWidth-4);
    y=Math.min(Math.max(y,4),bx.height-tip.offsetHeight-4);
    tip.style.left=x+"px"; tip.style.top=y+"px";
  });
  s.addEventListener("pointerleave",()=>tip.style.opacity=0);
}

/* ---------- 1. 인델 크기 스펙트럼 (결실 - / 삽입 +) ---------- */
function spectrum(boxId,tipId,data,lo,hi,logY){
  const box=document.getElementById(boxId),tip=document.getElementById(tipId);
  const W=1000,H=300,L=54,R=16,T=18,B=52;
  const bins=[]; for(let k=lo;k<=hi;k++){ if(k!==0) bins.push(k); }
  const val=k=> k<0 ? (data.del[-k]||0) : (data.ins[k]||0);
  const max=Math.max(1,...bins.map(val));
  const sc=v=> logY ? Math.log10(v+1)/Math.log10(max+1) : v/max;
  const bw=(W-L-R)/bins.length, py=v=>T+(1-sc(v))*(H-T-B);
  let g="",bars="";
  const ticks=4;
  const tickVals = logY
    ? Array.from({length:ticks+1},(_,i)=>Math.round(Math.pow(10,i*Math.log10(max+1)/ticks)-1))
    : Array.from({length:ticks+1},(_,i)=>max*i/ticks);
  for(let i=0;i<=ticks;i++){
    const v=tickVals[i];
    g+=`<line x1="${L}" y1="${py(v)}" x2="${W-R}" y2="${py(v)}" stroke="${css('--rule')}"/>`;
    g+=`<text x="${L-8}" y="${py(v)+4}" text-anchor="end" font-size="11.5" fill="${css('--ink-3')}"
         font-family="ui-monospace,monospace">${F(Math.round(v))}</text>`;
  }
  bins.forEach((k,i)=>{
    const v=val(k); if(v<=0) return;
    const x=L+i*bw, y=py(v), h=py(0)-y;
    const c=k<0?css('--del'):css('--ins');
    bars+=`<rect x="${(x+0.6).toFixed(1)}" y="${y.toFixed(1)}" width="${Math.max(bw-1.2,0.8).toFixed(1)}"
            height="${Math.max(h,1).toFixed(1)}" rx="${Math.min(2,bw/3).toFixed(1)}" fill="${c}"/>`;
  });
  const step = bins.length>40 ? 10 : (bins.length>16 ? 5 : 1);
  bins.forEach((k,i)=>{
    if(Math.abs(k)%step!==0) return;
    g+=`<text x="${(L+i*bw+bw/2).toFixed(1)}" y="${H-B+17}" text-anchor="middle" font-size="11"
        fill="${css('--ink-3')}" font-family="ui-monospace,monospace">${k>0?'+':''}${k}</text>`;
  });
  const zx=L+bins.findIndex(k=>k>0)*bw;
  g+=`<line x1="${zx}" y1="${T}" x2="${zx}" y2="${H-B}" stroke="${css('--rule-2')}" stroke-dasharray="3 3"/>`;
  g+=`<text x="${(L+zx)/2}" y="${H-B+36}" text-anchor="middle" font-size="12" fill="${css('--del')}">← 결실 (bp)</text>`;
  g+=`<text x="${(zx+W-R)/2}" y="${H-B+36}" text-anchor="middle" font-size="12" fill="${css('--ins')}">삽입 (bp) →</text>`;
  box.insertAdjacentHTML("afterbegin",svg(W,H,g+bars,"인델 크기 스펙트럼"));
  hookTip(box,tip,(fx)=>{
    const i=Math.floor((fx*W-L)/bw); if(i<0||i>=bins.length) return null;
    const k=bins[i],v=val(k); if(!v) return null;
    return {html:`<b>${k>0?'+':''}${k} bp ${k<0?'결실':'삽입'}</b>
      <div class="r"><span>read</span><b>${F(v)}</b></div>`};
  });
}

/* ---------- 2. 분류 누적 막대 ---------- */
function stacks(boxId,tipId,rows){
  const box=document.getElementById(boxId),tip=document.getElementById(tipId);
  const keys=[["ins","작은 삽입","--ins"],["del","작은 결실","--del"],["clip","말단 삽입","--clip"],
              ["sub","치환","--sub"],["ldel","대형 결실","--ldel"],["wt","야생형","--wt"]];
  const W=1000,rh=34,L=96,R=52,T=8;
  const H=T+rows.length*rh+18;
  let g="";
  rows.forEach((r,i)=>{
    const y=T+i*rh+5, bh=rh-14;
    let acc=0;
    keys.forEach(([k,,tok])=>{
      const p=r.pct[k]||0; if(p<=0) return;
      const x=L+(W-L-R)*acc/100, w=(W-L-R)*p/100;
      g+=`<rect x="${x.toFixed(1)}" y="${y}" width="${Math.max(w-2,0.8).toFixed(1)}" height="${bh}"
           rx="3" fill="${css(tok)}"/>`;
      acc+=p;
    });
    g+=`<text x="${L-10}" y="${y+bh/2+4}" text-anchor="end" font-size="12.5" fill="${css('--ink-2')}"
        font-family="ui-monospace,monospace">${r.label}</text>`;
    g+=`<text x="${W-R+8}" y="${y+bh/2+4}" font-size="12" fill="${css('--ink-3')}"
        font-family="ui-monospace,monospace">${F(r.total)}</text>`;
  });
  box.insertAdjacentHTML("afterbegin",svg(W,H,g,"샘플별 분류 구성"));
  hookTip(box,tip,(fx,fy)=>{
    const i=Math.floor((fy*H-T)/rh); if(i<0||i>=rows.length) return null;
    const r=rows[i];
    return {html:`<b>${r.label}</b>`+keys.map(([k,lab,tok])=>
      `<div class="r"><span><i class="sw" style="background:${css(tok)}"></i> ${lab}</span>
       <b>${(r.pct[k]||0).toFixed(2)}%</b></div>`).join("")};
  });
}

/* ---------- 2b. 삽입 길이 구간 (수정 전/후) ---------- */
function insBins(boxId,tipId,series){
  const box=document.getElementById(boxId),tip=document.getElementById(tipId);
  const BINS=[[1,1],[2,3],[4,7],[8,15],[16,30],[31,60],[61,1e9]];
  const LAB=["1","2–3","4–7","8–15","16–30","31–60","61+"];
  const W=1000,H=280,L=58,R=16,T=18,B=54;
  const gw=(W-L-R)/BINS.length,n=series.length,bw=(gw-12)/n;
  const agg=series.map(s=>BINS.map(([a,b])=>{
    let t=0; for(const k in s.data){const v=+k; if(v>=a&&v<=b) t+=s.data[k];} return t;}));
  const max=Math.max(1,...agg.flat());
  const sc=v=>Math.log10(v+1)/Math.log10(max+1);
  const py=v=>T+(1-sc(v))*(H-T-B);
  let g="",bars="";
  for(let i=0;i<=4;i++){
    const v=Math.round(Math.pow(10,i*Math.log10(max+1)/4)-1);
    g+=`<line x1="${L}" y1="${py(v)}" x2="${W-R}" y2="${py(v)}" stroke="${css('--rule')}"/>`;
    g+=`<text x="${L-8}" y="${py(v)+4}" text-anchor="end" font-size="11.5" fill="${css('--ink-3')}"
        font-family="ui-monospace,monospace">${F(v)}</text>`;}
  BINS.forEach((b,bi)=>{
    series.forEach((s,si)=>{
      const v=agg[si][bi]; if(!v) return;
      const x=L+bi*gw+6+si*bw,y=py(v);
      bars+=`<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${Math.max(bw-2,1).toFixed(1)}"
              height="${(py(0)-y).toFixed(1)}" rx="3" fill="${css(s.tok)}"/>`;
    });
    g+=`<text x="${(L+bi*gw+gw/2).toFixed(1)}" y="${H-B+17}" text-anchor="middle" font-size="11.5"
        fill="${css('--ink-3')}" font-family="ui-monospace,monospace">${LAB[bi]}</text>`;
  });
  const x8=L+3*gw;
  g+=`<line x1="${x8}" y1="${T}" x2="${x8}" y2="${H-B}" stroke="${css('--rule-2')}" stroke-dasharray="4 3"/>`;
  g+=`<text x="${x8+6}" y="${T+13}" font-size="11.5" fill="${css('--ink-2')}">← 수정 전에는 크기를 붙이지 못하던 구간</text>`;
  g+=`<text x="${(L+W-R)/2}" y="${H-B+38}" text-anchor="middle" font-size="12" fill="${css('--ink-2')}">삽입 길이 (bp) · y축 로그</text>`;
  box.insertAdjacentHTML("afterbegin",svg(W,H,g+bars,"삽입 길이 구간별 read 수"));
  hookTip(box,tip,(fx)=>{
    const bi=Math.floor((fx*W-L)/gw); if(bi<0||bi>=BINS.length) return null;
    return {html:`<b>${LAB[bi]} bp 삽입</b>`+series.map((s,si)=>
      `<div class="r"><span><i class="sw" style="background:${css(s.tok)}"></i> ${s.label}</span>
       <b>${F(agg[si][bi])}</b></div>`).join("")};
  });
}

/* ---------- 3. 대형 결실 크기 분포 ---------- */
function ldelChart(boxId,tipId,series){
  const box=document.getElementById(boxId),tip=document.getElementById(tipId);
  const BINS=[[50,99],[100,249],[250,499],[500,999],[1000,1999],[2000,3999],[4000,1e9]];
  const LAB=["50–99","100–249","250–499","500–999","1k–2k","2k–4k","4k+"];
  const W=1000,H=290,L=54,R=16,T=18,B=54;
  const gw=(W-L-R)/BINS.length, n=series.length, bw=(gw-10)/n;
  const agg=series.map(s=>BINS.map(([a,b])=>{
    let t=0; for(const k in s.data){const v=+k; if(v>=a&&v<=b) t+=s.data[k];} return t;}));
  const max=Math.max(1,...agg.flat());
  const py=v=>T+(1-v/max)*(H-T-B);
  let g="",bars="";
  for(let i=0;i<=4;i++){const v=max*i/4;
    g+=`<line x1="${L}" y1="${py(v)}" x2="${W-R}" y2="${py(v)}" stroke="${css('--rule')}"/>`;
    g+=`<text x="${L-8}" y="${py(v)+4}" text-anchor="end" font-size="11.5" fill="${css('--ink-3')}"
        font-family="ui-monospace,monospace">${F(Math.round(v))}</text>`;}
  BINS.forEach((b,bi)=>{
    series.forEach((s,si)=>{
      const v=agg[si][bi]; if(!v) return;
      const x=L+bi*gw+5+si*bw, y=py(v);
      bars+=`<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${Math.max(bw-2,1).toFixed(1)}"
              height="${(py(0)-y).toFixed(1)}" rx="3" fill="${css(s.tok)}"/>`;
    });
    g+=`<text x="${(L+bi*gw+gw/2).toFixed(1)}" y="${H-B+17}" text-anchor="middle" font-size="11.5"
        fill="${css('--ink-3')}" font-family="ui-monospace,monospace">${LAB[bi]}</text>`;
  });
  const x2k=L+5*gw;
  g+=`<line x1="${x2k}" y1="${T}" x2="${x2k}" y2="${H-B}" stroke="${css('--sub')}" stroke-dasharray="4 3"/>`;
  g+=`<text x="${x2k+6}" y="${T+13}" font-size="11.5" fill="${css('--sub')}">2 kb 이상: 재현성 주의</text>`;
  g+=`<text x="${(L+W-R)/2}" y="${H-B+38}" text-anchor="middle" font-size="12" fill="${css('--ink-2')}">결실 길이 (bp)</text>`;
  box.insertAdjacentHTML("afterbegin",svg(W,H,g+bars,"대형 결실 크기 분포"));
  hookTip(box,tip,(fx)=>{
    const bi=Math.floor((fx*W-L)/gw); if(bi<0||bi>=BINS.length) return null;
    return {html:`<b>${LAB[bi]} bp</b>`+series.map((s,si)=>
      `<div class="r"><span><i class="sw" style="background:${css(s.tok)}"></i> ${s.label}</span>
       <b>${F(agg[si][bi])}</b></div>`).join("")};
  });
}
"""

def open_any(path):
    """.gz/.bz2/.xz 도 그대로 읽습니다."""
    if path.endswith(".gz"):
        import gzip; return gzip.open(path, "rt", errors="replace")
    if path.endswith(".bz2"):
        import bz2; return bz2.open(path, "rt", errors="replace")
    if path.endswith(".xz"):
        import lzma; return lzma.open(path, "rt", errors="replace")
    return open(path, errors="replace")

RCT = str.maketrans("ACGTN", "TGCAN")

def resolve(p, hints):
    if p and os.path.exists(p): return p
    for h in hints:
        for cand in (os.path.join(h, p or ""), os.path.join(h, os.path.basename(p or ""))):
            if cand and os.path.exists(cand): return cand
    return None

def scan_fastq(paths, queries, limit=0):
    """원본 read 에서 각 서열의 출현 횟수를 셉니다 (양쪽 가닥)."""
    pat = []
    for q in queries:
        pat.append((q, q, q[::-1].translate(RCT)))
    cnt = collections.Counter(); total = 0
    for p in paths:
        try: fh = open_any(p)
        except OSError: continue
        with fh:
            for i, line in enumerate(fh):
                if limit and i >= limit*4: break
                if i & 3 != 1: continue
                total += 1
                for key, fwd, rev in pat:
                    if fwd in line or rev in line: cnt[key] += 1
    return cnt, total

def pct_map(c):
    tot = c.get("total_analyzed", 0) or 1
    return {k: 100.0*c.get(src, 0)/tot for src, _, k in CAT}

def build(samples, title, outpath, scan=None, comp=None):
    E = html.escape
    rows = [{"label": s["name"].replace("PD1_", ""), "total": s["summary"]["treated"].get("total_analyzed", 0),
             "pct": pct_map(s["summary"]["treated"])} for s in samples]
    agg_ins, agg_del = collections.Counter(), collections.Counter()
    agg_ld, agg_ldc = collections.Counter(), collections.Counter()
    seqs = collections.Counter(); ctr_seqs = collections.Counter()
    for s in samples:
        agg_ins.update({int(k): v for k, v in s["mut_ins"].items()})
        agg_del.update({int(k): v for k, v in s["mut_del"].items()})
        agg_ld.update({int(k): v for k, v in s["ldel"].items()})
        agg_ldc.update({int(k): v for k, v in s["ldel_ctr"].items()})
        for q, n in s["mut_ins_seq"]: seqs[q] += n
        ctr_seqs.update(s["ctr_ins_seq"])
    p0 = samples[0]["params"]
    tot_analyzed = sum(r["total"] for r in rows)
    tot_aligned = sum(s["summary"]["treated"].get("aligned_pairs_all_positions", 0) for s in samples)
    edited = sum(r["total"]*(100-r["pct"]["wt"])/100 for r in rows)
    big_ins = sum(v for k, v in agg_ins.items() if k >= 8)
    maxi = max(agg_ins) if agg_ins else 0

    data = {"ins": {str(k): v for k, v in agg_ins.items()},
            "del": {str(k): v for k, v in agg_del.items()}}
    ld_series = [{"label": "처리군", "tok": "--ldel", "data": {str(k): v for k, v in agg_ld.items()}},
                 {"label": "대조군", "tok": "--wt", "data": {str(k): v for k, v in agg_ldc.items()}}]

    seq_rows = ""; scan_note = ""
    if scan:
        mc, mt = scan["mut"]; cc, ct = scan["ctr"]
        scan_note = (f'처리군 원본 {mt:,} read, 대조군 원본 {ct:,} read 를 직접 세었습니다. '
                     f'백만 read 당 출현 빈도로 비교합니다.')
        for q, n in seqs.most_common(12):
            mppm = 1e6*mc.get(q,0)/max(mt,1); cppm = 1e6*cc.get(q,0)/max(ct,1)
            art = cppm >= 0.5*mppm and cppm > 0
            seq_rows += (f'<tr><td class="n">{len(q)}</td><td class="n">{n}</td>'
                         f'<td class="n">{mppm:,.0f}</td><td class="n">{cppm:,.0f}</td>'
                         f'<td><span class="flag {"" if art else "ok"}">'
                         f'{"대조군에도 있음 · 아티팩트" if art else "처리군 특이"}</span></td>'
                         f'<td class="seq">{E(q)}</td></tr>')
        hdr = ('<tr><th class="n">길이</th><th class="n">삽입 호출</th>'
               '<th class="n">처리군 원본 /M</th><th class="n">대조군 원본 /M</th>'
               '<th>판정</th><th>서열</th></tr>')
        ncol = 6
    else:
        scan_note = ('원본 FASTQ 를 찾지 못해 분류된 read 안에서만 비교했습니다. '
                     '<b>대조군 원본을 직접 세지 않으면 아티팩트를 놓칠 수 있습니다</b> — '
                     '<code>--fastq-dir</code> 로 FASTQ 위치를 알려주십시오.')
        for q, n in seqs.most_common(12):
            c = ctr_seqs.get(q, 0)
            art = c >= max(1, n*0.5)
            seq_rows += (f'<tr><td class="n">{len(q)}</td><td class="n">{n}</td><td class="n">{c}</td>'
                         f'<td><span class="flag {"" if art else "ok"}">'
                         f'{"대조군에도 있음" if art else "미확인"}</span></td>'
                         f'<td class="seq">{E(q)}</td></tr>')
        hdr = ('<tr><th class="n">길이</th><th class="n">처리군</th><th class="n">대조군(분류분)</th>'
               '<th>판정</th><th>서열</th></tr>')
        ncol = 5
    if not seq_rows:
        seq_rows = f'<tr><td colspan="{ncol}" style="color:var(--ink-3)">8 bp 이상 삽입 없음</td></tr>'

    cls_rows = ""
    for r in rows:
        cls_rows += f'<tr><td>{E(r["label"])}</td>' + "".join(
            f'<td class="n">{r["pct"][k]:.2f}%</td>' for _, _, k in CAT) + f'<td class="n">{r["total"]:,}</td></tr>'

    meta = " ".join(f"<span><b>{k}</b> {E(str(v))}</span>" for k, v in
                    [("gene", p0.get("gene", "-")), ("targetst", p0.get("targetst", "-")),
                     ("cle", p0.get("cle", "-")), ("K", p0.get("K", "-")),
                     ("dist_cut", p0.get("dist_cut", "-")), ("anchor_min", p0.get("anchor_min", "-")),
                     ("cut_between", p0.get("cut_between", "-")),
                     ("ins_fix", p0.get("ins_fix", "1")), ("ins_max", p0.get("ins_max", "0"))])

    comp_html = ""
    if comp:
        cins = collections.Counter()
        for s_ in comp: cins.update({int(k): v for k, v in s_["mut_ins"].items()})
        crows = [{"label": s_["name"].replace("PD1_", ""), "pct": pct_map(s_["summary"]["treated"])}
                 for s_ in comp]
        cmap = {r["label"].replace("_OFF", "").replace("_ON", ""): r for r in crows}
        trs = ""
        for r in rows:
            base = r["label"].replace("_ON", "").replace("_OFF", "")
            o = cmap.get(base)
            if not o: continue
            trs += f'<tr><td>{E(base)}</td>' + "".join(
                f'<td class="n">{o["pct"][k]:.2f} → {r["pct"][k]:.2f}</td>' for _, _, k in CAT) + "</tr>"
        big_before = sum(v for k, v in cins.items() if k >= 8)
        mx_before = max(cins) if cins else 0
        comp_html = f"""
<section>
  <h2>삽입 정렬 수정 전후</h2>
  <p class="note">같은 데이터를 <code>ins_fix=0</code> 과 <code>ins_fix=1</code> 로 두 번 돌린 결과입니다.</p>
  <div class="tiles">
    <div class="tile"><span class="k">최대 삽입 길이</span><span class="v">{mx_before} → {maxi} bp</span><span class="n">CIGAR 의 I</span></div>
    <div class="tile"><span class="k">8 bp 이상 삽입</span><span class="v">{big_before:,} → {big_ins:,}</span><span class="n">read 수</span></div>
  </div>
  <div class="legend">
    <span><i class="sw" style="background:var(--wt)"></i> 수정 전</span>
    <span><i class="sw" style="background:var(--ins)"></i> 수정 후</span>
  </div>
  <div class="chart" id="b-ins"><div class="tip" id="t-ins"></div></div>
  <details open><summary>분류 비율 전 → 후 (%)</summary><div class="tblwrap"><table>
    <thead><tr><th>샘플</th>{''.join(f'<th class="n">{ko}</th>' for _, ko, _ in CAT)}</tr></thead>
    <tbody>{trs}</tbody></table></div></details>
</section>"""
        comp_js = (f'insBins("b-ins","t-ins",[{{"label":"수정 전","tok":"--wt","data":'
                   f'{json.dumps({str(k): v for k, v in cins.items()})}}},'
                   f'{{"label":"수정 후","tok":"--ins","data":{json.dumps(data["ins"])}}}]);')
    else:
        comp_js = ""

    doc = f"""<!DOCTYPE html><html lang="ko"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{E(title)}</title><style>{CSS}</style></head><body><div class="wrap">
<header>
  <h1>{E(title)}</h1>
  <p class="sub">kmer_align 실행 결과 · 샘플 {len(samples)}개 · 생성 {datetime.date.today().isoformat()}</p>
  <div class="meta">{meta}</div>
</header>

<section>
  <h2>요약</h2>
  <div class="tiles">
    <div class="tile"><span class="k">정렬된 pair</span><span class="v">{tot_aligned:,}</span><span class="n">전체 위치</span></div>
    <div class="tile"><span class="k">분류 대상</span><span class="v">{tot_analyzed:,}</span><span class="n">절단 부위를 지나는 pair</span></div>
    <div class="tile"><span class="k">편집된 비율</span><span class="v">{100.0*edited/max(tot_analyzed,1):.1f}%</span><span class="n">야생형이 아닌 것</span></div>
    <div class="tile"><span class="k">최대 삽입 길이</span><span class="v">{maxi} bp</span><span class="n">CIGAR 의 I</span></div>
    <div class="tile"><span class="k">8 bp 이상 삽입</span><span class="v">{big_ins:,}</span><span class="n">read 수</span></div>
  </div>
</section>

<section>
  <h2>인델 크기 스펙트럼</h2>
  <p class="note">절단 부위 ±13 bp 안의 인델을 길이별로 센 것입니다. 왼쪽이 결실, 오른쪽이 삽입.</p>
  <div class="legend">
    <span><i class="sw" style="background:var(--del)"></i> 결실</span>
    <span><i class="sw" style="background:var(--ins)"></i> 삽입</span>
  </div>
  <div class="chart" id="b-spec"><div class="tip" id="t-spec"></div></div>
  <div class="tabs" role="tablist">
    <button role="tab" aria-selected="true" data-lo="-30" data-hi="30">±30 bp</button>
    <button role="tab" aria-selected="false" data-lo="-10" data-hi="10">±10 bp</button>
    <button role="tab" aria-selected="false" data-lo="-49" data-hi="100">전체</button>
    <button id="logbtn" aria-pressed="false" style="margin-left:10px">y축 로그</button>
  </div>
</section>

<section>
  <h2>샘플별 분류 구성</h2>
  <div class="legend">
    <span><i class="sw" style="background:var(--ins)"></i> 작은 삽입</span>
    <span><i class="sw" style="background:var(--del)"></i> 작은 결실</span>
    <span><i class="sw" style="background:var(--clip)"></i> 말단 삽입</span>
    <span><i class="sw" style="background:var(--sub)"></i> 치환</span>
    <span><i class="sw" style="background:var(--ldel)"></i> 대형 결실</span>
    <span><i class="sw" style="background:var(--wt)"></i> 야생형</span>
  </div>
  <div class="chart" id="b-stack"><div class="tip" id="t-stack"></div></div>
  <details><summary>표로 보기</summary><div class="tblwrap"><table>
    <thead><tr><th>샘플</th>{''.join(f'<th class="n">{ko}</th>' for _,ko,_ in CAT)}<th class="n">합계</th></tr></thead>
    <tbody>{cls_rows}</tbody></table></div></details>
</section>

<section>
  <h2>대형 결실 크기 분포</h2>
  <p class="note">split read 의 CIGAR 에 적힌 결실 길이입니다. 대조군을 나란히 두어 아티팩트 수준을 함께 봅니다.</p>
  <div class="legend">
    <span><i class="sw" style="background:var(--ldel)"></i> 처리군</span>
    <span><i class="sw" style="background:var(--wt)"></i> 대조군</span>
  </div>
  <div class="chart" id="b-ldel"><div class="tip" id="t-ldel"></div></div>
  <div class="warn"><b>2 kb 이상 구간 주의.</b> only/common 분할은 입력에 민감해 단일 실행으로 재현되지 않습니다.
   입력의 1%만 달라져도 최대 30 pp 움직입니다. 이 구간을 보고할 때는 반복 표본으로 변동폭을 함께 제시하십시오.</div>
</section>

<section>
  <h2>8 bp 이상 삽입 서열</h2>
  <p class="note">대조군에서 같은 서열이 비슷하게 나오면 CRISPR 산물이 아니라 공동 증폭된 배경 서열입니다.
   파이프라인의 only/common 대조군 필터는 <b>대형 결실에만</b> 적용되고 삽입에는 적용되지 않으므로, 여기서 직접 확인해야 합니다.</p>
  <p class="note">{scan_note}</p>
  <div class="tblwrap"><table>
    <thead>{hdr}</thead>
    <tbody>{seq_rows}</tbody></table></div>
</section>

{comp_html}
<footer>kmer_align · CRISPR-LargeDel</footer>
</div>
<script>
const DATA={json.dumps(data)}, ROWS={json.dumps(rows)}, LD={json.dumps(ld_series)};
{JS}
stacks("b-stack","t-stack",ROWS);
ldelChart("b-ldel","t-ldel",LD);
{comp_js}
let RANGE=[-30,30], LOGY=false;
function drawSpec(){{
  const b=document.getElementById("b-spec");
  b.querySelectorAll("svg").forEach(e=>e.remove());
  spectrum("b-spec","t-spec",DATA,RANGE[0],RANGE[1],LOGY);
}}
drawSpec();
document.querySelectorAll('.tabs button[data-lo]').forEach(btn=>btn.addEventListener('click',()=>{{
  document.querySelectorAll('.tabs button[data-lo]').forEach(b=>b.setAttribute('aria-selected','false'));
  btn.setAttribute('aria-selected','true');
  RANGE=[+btn.dataset.lo,+btn.dataset.hi]; drawSpec();
}}));
const lb=document.getElementById('logbtn');
lb.addEventListener('click',()=>{{
  LOGY=!LOGY; lb.setAttribute('aria-pressed',LOGY);
  lb.setAttribute('aria-selected',LOGY); drawSpec();
}});
</script></body></html>"""
    with open(outpath, "w", encoding="utf-8") as f:
        f.write(doc)
    return outpath

def main():
    ap = argparse.ArgumentParser(description="kmer_align 결과를 단일 HTML 리포트로 묶습니다.")
    ap.add_argument("dirs", nargs="+", help="결과 폴더 (여러 개 가능)")
    ap.add_argument("-o", "--out", default="ka_report.html")
    ap.add_argument("-t", "--title", default="kmer_align 분석 리포트")
    ap.add_argument("--fastq-dir", action="append", default=[],
                    help="원본 FASTQ 가 있는 폴더 (여러 번 지정 가능). 8 bp 이상 삽입이 "
                         "대조군에도 있는지 직접 세는 데 씁니다.")
    ap.add_argument("--compare", nargs="+", default=[], metavar="DIR",
                    help="수정 전(ins_fix=0) 결과 폴더들. 주면 전후 비교 섹션이 붙습니다.")
    ap.add_argument("--no-scan", action="store_true", help="원본 FASTQ 스캔 생략")
    ap.add_argument("--scan-limit", type=int, default=0, help="스캔할 read 수 상한 (0=전부)")
    a = ap.parse_args()
    samples = []
    for d in a.dirs:
        s = collect(d)
        if s is None:
            print(f"[skip] {d}: mutation_summary.tsv 없음", file=sys.stderr); continue
        samples.append(s)
    if not samples:
        print("[error] 읽을 결과가 없습니다", file=sys.stderr); sys.exit(1)
    scan = None
    if not a.no_scan:
        hints = list(a.fastq_dir) + [os.path.dirname(d.rstrip("/")) or "." for d in a.dirs] \
                + [os.path.join(os.path.dirname(d.rstrip("/")) or ".", "..", "work") for d in a.dirs]
        q = collections.Counter()
        for s_ in samples:
            for seq, n in s_["mut_ins_seq"]: q[seq] += n
        queries = [k for k, _ in q.most_common(12)]
        mut_p, ctr_p = [], []
        for s_ in samples:
            for key, bucket in (("mut_1", mut_p), ("mut_2", mut_p), ("ctr_1", ctr_p), ("ctr_2", ctr_p)):
                r = resolve(s_["params"].get(key), hints)
                if r and r not in bucket: bucket.append(r)
        if queries and mut_p and ctr_p:
            print(f"[scan] 원본 FASTQ {len(mut_p)+len(ctr_p)}개에서 서열 {len(queries)}종 확인 중...",
                  file=sys.stderr)
            scan = {"mut": scan_fastq(mut_p, queries, a.scan_limit),
                    "ctr": scan_fastq(ctr_p, queries, a.scan_limit)}
    comp = [x for x in (collect(d) for d in a.compare) if x] or None
    p = build(samples, a.title, a.out, scan, comp)
    print(f"[ok] {p}  ({os.path.getsize(p)/1024:.0f} KB, 샘플 {len(samples)}개)")

if __name__ == "__main__":
    main()
