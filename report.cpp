/* report.cpp - kmer_align HTML 리포트 (report.h 참고) */
#include "report.h"
#include "report_tpl.h"
#include "arg_parse.h"
#include "func.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
using namespace std;
extern Args arg;

#define KA_VERSION "2026.09"
#define SEQ_MIN     8     /* 이 길이 이상의 삽입만 서열을 모읍니다 */
#define SEQ_SCAN   24     /* 원본 FASTQ 에서 셀 후보 수 */

namespace {

struct Sample
{
	string name, label;
	vector<pair<string,string> > params;
	map<string,long long> cat[2];                     /* 0 control, 1 treated */
	map<int,long long> mut_ins, mut_del, ctr_ins, ctr_del, ldel, ldel_ctr;
	map<string,long long> seqs;                       /* 첫 I 서열 -> 레코드 수 */
	bool scanned = false;
	long long scan_reads[2] = {0,0};                  /* 0 control, 1 treated */
	map<string,long long> scan_hit[2];
};

Sample g_run;

/* ---------------------------------------------------------------- 유틸 */
string fmt_ll(long long v){ char b[32]; snprintf(b,sizeof b,"%lld",v); return b; }
string fmt_i(int v){ char b[16]; snprintf(b,sizeof b,"%d",v); return b; }

string json_str(const string &s)
{
	string o = "\"";
	for(size_t i=0;i<s.size();i++)
	{
		unsigned char c = (unsigned char)s[i];
		switch(c)
		{
			case '"':  o += "\\\""; break;
			case '\\': o += "\\\\"; break;
			case '\n': o += "\\n";  break;
			case '\r': o += "\\r";  break;
			case '\t': o += "\\t";  break;
			case '<':  o += "\\u003c"; break;   /* </script> 방지 */
			case '>':  o += "\\u003e"; break;
			case '&':  o += "\\u0026"; break;
			default:
				if(c < 0x20){ char b[8]; snprintf(b,sizeof b,"\\u%04x",c); o += b; }
				else o += (char)c;
		}
	}
	return o + "\"";
}
string html_esc(const string &s)
{
	string o;
	for(char c: s){ if(c=='<')o+="&lt;"; else if(c=='>')o+="&gt;"; else if(c=='&')o+="&amp;"; else if(c=='"')o+="&quot;"; else o+=c; }
	return o;
}
string tsv_clean(string s)
{
	for(char &c: s){ if(c=='\t'||c=='\n'||c=='\r') c=' '; }
	return s;
}
string today()
{
	time_t t = time(NULL); struct tm *tm = localtime(&t);
	char b[16]; if(tm) strftime(b,sizeof b,"%Y-%m-%d",tm); else strcpy(b,"-");
	return b;
}
string revcomp(const string &s)
{
	string r(s.rbegin(), s.rend());
	for(char &c: r){ switch(c){case 'A':c='T';break;case 'T':c='A';break;case 'C':c='G';break;case 'G':c='C';break;default:break;} }
	return r;
}

/* ---------------------------------------------------------------- JSON */
string json_hist(const map<int,long long> &m)
{
	string o = "{"; bool first = true;
	for(auto &kv: m){ if(!first) o += ","; first=false; o += "\"" + fmt_i(kv.first) + "\":" + fmt_ll(kv.second); }
	return o + "}";
}
string json_sample(const Sample &s)
{
	string o = "{\"name\":" + json_str(s.name) + ",\"label\":" + json_str(s.label) + ",\"params\":{";
	for(size_t i=0;i<s.params.size();i++)
	{
		if(i) o += ",";
		o += json_str(s.params[i].first) + ":" + json_str(s.params[i].second);
	}
	o += "}";
	const char *side[2] = {"control","treated"};
	for(int k=0;k<2;k++)
	{
		o += string(",\"") + side[k] + "\":{";
		bool first = true;
		for(auto &kv: s.cat[k]){ if(!first) o += ","; first=false; o += json_str(kv.first) + ":" + fmt_ll(kv.second); }
		o += "}";
	}
	o += ",\"mut_ins\":" + json_hist(s.mut_ins) + ",\"mut_del\":" + json_hist(s.mut_del);
	o += ",\"ctr_ins\":" + json_hist(s.ctr_ins) + ",\"ctr_del\":" + json_hist(s.ctr_del);
	o += ",\"ldel\":" + json_hist(s.ldel) + ",\"ldel_ctr\":" + json_hist(s.ldel_ctr);
	/* 서열은 많이 나온 순 상위만 */
	vector<pair<string,long long> > v(s.seqs.begin(), s.seqs.end());
	sort(v.begin(), v.end(), [](const pair<string,long long>&a,const pair<string,long long>&b){
		return a.second!=b.second ? a.second>b.second : a.first<b.first; });
	if(v.size() > 200) v.resize(200);
	o += ",\"seqs\":[";
	for(size_t i=0;i<v.size();i++){ if(i) o += ","; o += "[" + json_str(v[i].first) + "," + fmt_ll(v[i].second) + "]"; }
	o += "]";
	if(s.scanned)
	{
		o += ",\"scan\":{\"ctr_reads\":" + fmt_ll(s.scan_reads[0]) + ",\"mut_reads\":" + fmt_ll(s.scan_reads[1]);
		const char *key[2] = {"ctr","mut"};
		for(int k=0;k<2;k++)
		{
			o += string(",\"") + key[k] + "\":{"; bool first = true;
			for(auto &kv: s.scan_hit[k]){ if(!first) o += ","; first=false; o += json_str(kv.first) + ":" + fmt_ll(kv.second); }
			o += "}";
		}
		o += "}";
	}
	else o += ",\"scan\":null";
	return o + "}";
}

bool write_html(const vector<Sample> &S, const vector<Sample> &C, const string &title, const string &path)
{
	FILE *f = fopen(path.c_str(),"w");
	if(f==NULL){ fprintf(stderr,"[warning] could not write '%s'\n",path.c_str()); return false; }
	string R = "{\"title\":" + json_str(title) + ",\"date\":" + json_str(today()) +
	           ",\"version\":" + json_str(KA_VERSION) + ",\"samples\":[";
	for(size_t i=0;i<S.size();i++){ if(i) R += ","; R += json_sample(S[i]); }
	R += "],\"compare\":[";
	for(size_t i=0;i<C.size();i++){ if(i) R += ","; R += json_sample(C[i]); }
	R += "]}";
	fprintf(f,"<!DOCTYPE html>\n<html lang=\"ko\"><head><meta charset=\"utf-8\">\n"
	          "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\n"
	          "<title>%s</title>\n<style>%s</style></head>\n<body><div class=\"wrap\" id=\"app\"></div>\n"
	          "<script>const R=%s;</script>\n<script>%s</script>\n</body></html>\n",
	        html_esc(title).c_str(), KA_CSS, R.c_str(), KA_JS);
	fclose(f);
	return true;
}

/* ---------------------------------------------------------------- 집계 TSV (합치기용) */
void write_hist(FILE *f, const char *key, const map<int,long long> &m)
{
	for(auto &kv: m) fprintf(f,"h\t%s\t%d\t%lld\n",key,kv.first,kv.second);
}
bool write_stats(const Sample &s, const string &path)
{
	FILE *f = fopen(path.c_str(),"w");
	if(f==NULL){ fprintf(stderr,"[warning] could not write '%s'\n",path.c_str()); return false; }
	fprintf(f,"#kmer_align_report\t1\n");
	fprintf(f,"name\t%s\nlabel\t%s\n",tsv_clean(s.name).c_str(),tsv_clean(s.label).c_str());
	for(auto &p: s.params) fprintf(f,"param\t%s\t%s\n",tsv_clean(p.first).c_str(),tsv_clean(p.second).c_str());
	const char *side[2] = {"control","treated"};
	for(int k=0;k<2;k++) for(auto &kv: s.cat[k]) fprintf(f,"cat\t%s\t%s\t%lld\n",side[k],kv.first.c_str(),kv.second);
	write_hist(f,"mut_ins",s.mut_ins); write_hist(f,"mut_del",s.mut_del);
	write_hist(f,"ctr_ins",s.ctr_ins); write_hist(f,"ctr_del",s.ctr_del);
	write_hist(f,"ldel",s.ldel);       write_hist(f,"ldel_ctr",s.ldel_ctr);
	for(auto &kv: s.seqs) fprintf(f,"seq\t%s\t%lld\n",kv.first.c_str(),kv.second);
	if(s.scanned)
	{
		fprintf(f,"scan\tctr_reads\t%lld\nscan\tmut_reads\t%lld\n",s.scan_reads[0],s.scan_reads[1]);
		const char *key[2] = {"ctr","mut"};
		for(int k=0;k<2;k++) for(auto &kv: s.scan_hit[k]) fprintf(f,"hit\t%s\t%s\t%lld\n",key[k],kv.first.c_str(),kv.second);
	}
	fclose(f);
	return true;
}
vector<string> split_tab(const string &line)
{
	vector<string> v; size_t st=0;
	for(size_t i=0;i<=line.size();i++) if(i==line.size()||line[i]=='\t'){ v.push_back(line.substr(st,i-st)); st=i+1; }
	return v;
}
bool read_stats(const string &path, Sample &s)
{
	FILE *f = fopen(path.c_str(),"r");
	if(f==NULL) return false;
	char buf[1<<14]; bool ok=false;
	while(fgets(buf,sizeof buf,f))
	{
		string line(buf);
		while(!line.empty() && (line.back()=='\n'||line.back()=='\r')) line.pop_back();
		vector<string> v = split_tab(line);
		if(v.empty()) continue;
		const string &t = v[0];
		if(t=="#kmer_align_report"){ ok=true; }
		else if(t=="name" && v.size()>=2) s.name=v[1];
		else if(t=="label" && v.size()>=2) s.label=v[1];
		else if(t=="param" && v.size()>=3) s.params.push_back(make_pair(v[1],v[2]));
		else if(t=="cat" && v.size()>=4) s.cat[v[1]=="treated"?1:0][v[2]] = atoll(v[3].c_str());
		else if(t=="h" && v.size()>=4)
		{
			map<int,long long> *m = NULL;
			if(v[1]=="mut_ins")m=&s.mut_ins; else if(v[1]=="mut_del")m=&s.mut_del;
			else if(v[1]=="ctr_ins")m=&s.ctr_ins; else if(v[1]=="ctr_del")m=&s.ctr_del;
			else if(v[1]=="ldel")m=&s.ldel; else if(v[1]=="ldel_ctr")m=&s.ldel_ctr;
			if(m) (*m)[atoi(v[2].c_str())] += atoll(v[3].c_str());
		}
		else if(t=="seq" && v.size()>=3) s.seqs[v[1]] += atoll(v[2].c_str());
		else if(t=="scan" && v.size()>=3){ s.scanned=true; s.scan_reads[v[1]=="mut_reads"?1:0]=atoll(v[2].c_str()); }
		else if(t=="hit" && v.size()>=4) s.scan_hit[v[1]=="mut"?1:0][v[2]] = atoll(v[3].c_str());
	}
	fclose(f);
	return ok;
}

/* ---------------------------------------------------------------- 원본 FASTQ 스캔
   8 bp 이상 삽입 서열이 대조군 원본에도 있는지 셉니다.  파이프라인의 only/common
   대조군 필터는 대형 결실에만 적용되므로, 삽입은 여기서 확인하지 않으면 공동
   증폭된 배경 서열이 그대로 "삽입 변이" 로 보고됩니다. */
void scan_fastq(const char *path, const vector<string> &q, vector<long long> &hit, long long &reads)
{
	if(path==NULL || path[0]=='\0') return;
	int pipe=0;
	FILE *f = ka_open(path,&pipe);
	if(f==NULL){ fprintf(stderr,"[report] cannot open '%s' for the insertion check\n",path); return; }
	vector<string> rc(q.size());
	for(size_t k=0;k<q.size();k++) rc[k]=revcomp(q[k]);
	char line[4096]; long long ln=0;
	while(fgets(line,sizeof line,f))
	{
		if((ln++ & 3) != 1) continue;
		reads++;
		for(size_t k=0;k<q.size();k++)
			if(strstr(line,q[k].c_str()) || strstr(line,rc[k].c_str())) hit[k]++;
	}
	ka_close(f,pipe);
}

} /* namespace */

/* ================================================================ 실행 중 집계 */
int report_role_for(const string &b)
{
	if(b.empty() || !arg.report) return RR_NONE;
	if(b==arg.mutsite_sm || b==arg.mutsite_clip) return RR_MUT_SITE;
	if(b==arg.ctrsite_sm || b==arg.ctrsite_clip) return RR_CTR_SITE;
	if(b==arg.resmutonly || b==arg.resmutcom)    return RR_MUT_LDEL;
	if(b==arg.resctronly || b==arg.resctrcom)    return RR_CTR_LDEL;
	return RR_NONE;
}

void report_observe(int role, const string &line)
{
	if(role==RR_NONE) return;
	/* SAM 필드 6(CIGAR), 10(SEQ) */
	size_t cs=0,ce=0,ss=0,se=0,st=0; int fld=0;
	for(size_t i=0;i<=line.size();i++)
	{
		if(i==line.size() || line[i]=='\t')
		{
			if(fld==5){cs=st;ce=i;}
			else if(fld==9){ss=st;se=i;break;}
			fld++; st=i+1;
		}
	}
	if(ce<=cs) return;
	/* itot: 모든 I 의 합 (스펙트럼용).  서열 표에는 가장 긴 I 하나를 씁니다 -
	   합으로 8 bp 를 넘겨도 첫 I 가 1 bp 면 "C" 같은 한 글자가 표에 올라와
	   모든 read 에 들어 있는 것으로 셉니다. */
	long long itot=0, n=0, q=0, best=0; int dmax=0; string ins;
	const long long L=(long long)(se-ss);
	for(size_t i=cs;i<ce;i++)
	{
		char c=line[i];
		if(c>='0'&&c<='9'){ n=n*10+(c-'0'); continue; }
		if(c=='I')
		{
			itot+=n;
			if(n>best && L>0 && q<L){ best=n; ins=line.substr(ss+q,(size_t)min(n,L-q)); }
		}
		else if(c=='D'){ if(n>dmax) dmax=(int)n; }
		if(c=='M'||c=='I'||c=='S'||c=='='||c=='X') q+=n;   /* query 를 소비하는 연산 */
		n=0;
	}
	switch(role)
	{
		case RR_MUT_SITE:
			if(itot>0) g_run.mut_ins[(int)itot]++;
			if(dmax>0) g_run.mut_del[dmax]++;
			if(best>=SEQ_MIN && (long long)ins.size()>=SEQ_MIN) g_run.seqs[ins]++;
			break;
		case RR_CTR_SITE:
			if(itot>0) g_run.ctr_ins[(int)itot]++;
			if(dmax>0) g_run.ctr_del[dmax]++;
			break;
		case RR_MUT_LDEL: if(dmax>0) g_run.ldel[dmax]++;     break;
		case RR_CTR_LDEL: if(dmax>0) g_run.ldel_ctr[dmax]++; break;
		default: break;
	}
}

static void put_pat(map<string,long long> &m, const ka_pat &p)
{
	m["wt"]=p.wt; m["small_insertion"]=p.ins; m["small_deletion"]=p.del;
	m["terminal_insertion_clipped"]=p.clip; m["substitution"]=p.sub; m["large_deletion"]=p.ldel;
	m["total_analyzed"]=p.total; m["not_covered_excluded"]=p.nocov;
	m["large_deletion_common_excluded"]=p.lcommon; m["aligned_pairs_all_positions"]=p.aligned;
}

void report_finish(const ka_pat &ctr, const ka_pat &mut)
{
	if(!arg.report) return;
	Sample &s = g_run;
	s.name  = arg.genename + "_" + arg.samplename;
	s.label = arg.samplename;
	put_pat(s.cat[0],ctr); put_pat(s.cat[1],mut);
	int cut = arg.targetst + arg.cle;
	s.params = {
		{"gene",arg.genename},{"sample",arg.samplename},{"ref",arg.refname?arg.refname:""},
		{"ref_len",fmt_i(arg.ref_len)},{"targetst",fmt_i(arg.targetst)},{"targeted",fmt_i(arg.targeted)},
		{"cle",fmt_i(arg.cle)},{"cut",fmt_i(cut)},{"K",fmt_i(arg.K)},{"dist_cut",fmt_i(arg.dist_cut)},
		{"min_gap",fmt_i(arg.min_gap)},{"anchor_min",fmt_i(arg.anchor_min)},{"cut_between",fmt_i(arg.cut_between)},
		{"ins_fix",fmt_i(arg.ins_fix)},{"ins_max",fmt_i(arg.ins_max)},{"large_del_len",fmt_i(arg.large_del_len)},
		{"ctr_1",arg.ctrfastq_l?arg.ctrfastq_l:""},{"ctr_2",arg.ctrfastq_r?arg.ctrfastq_r:""},
		{"mut_1",arg.mutfastq_l?arg.mutfastq_l:""},{"mut_2",arg.mutfastq_r?arg.mutfastq_r:""}};

	/* 원본 FASTQ 스캔 - 후보가 있을 때만 */
	vector<pair<string,long long> > v(s.seqs.begin(), s.seqs.end());
	sort(v.begin(), v.end(), [](const pair<string,long long>&a,const pair<string,long long>&b){
		return a.second!=b.second ? a.second>b.second : a.first<b.first; });
	if(v.size() > SEQ_SCAN) v.resize(SEQ_SCAN);
	if(!v.empty())
	{
		vector<string> q; for(auto &p: v) q.push_back(p.first);
		vector<long long> hc(q.size(),0), hm(q.size(),0);
		clock_t t0 = clock();
		scan_fastq(arg.ctrfastq_l,q,hc,s.scan_reads[0]); scan_fastq(arg.ctrfastq_r,q,hc,s.scan_reads[0]);
		scan_fastq(arg.mutfastq_l,q,hm,s.scan_reads[1]); scan_fastq(arg.mutfastq_r,q,hm,s.scan_reads[1]);
		if(s.scan_reads[0]>0 && s.scan_reads[1]>0)
		{
			s.scanned = true;
			for(size_t k=0;k<q.size();k++){ s.scan_hit[0][q[k]]=hc[k]; s.scan_hit[1][q[k]]=hm[k]; }
		}
		printf("report: insertion check over %lld + %lld raw reads, %.1lf s\n",
		       s.scan_reads[0],s.scan_reads[1],(double)(clock()-t0)/CLOCKS_PER_SEC);
	}

	string base = arg.outdir + "/" + s.name + "_report";
	vector<Sample> S(1,s), C;
	string title = s.name + " — kmer_align";
	bool a = write_stats(s, base+".tsv");
	bool b = write_html(S, C, title, base+".html");
	if(a && b) printf("report written to %s.html\n", base.c_str());
}

/* ================================================================ `kmer_align report` */
static void report_usage(const char *prog)
{
	fprintf(stderr,
	"usage: %s report [-o OUT.html] [-t TITLE] [--compare DIR...] [--] DIR...\n"
	"\n"
	"  Each DIR is a kmer_align output directory (or its *_report.tsv file).\n"
	"  All DIRs are combined into one HTML report.\n"
	"\n"
	"  -o FILE          output path (default kmer_align_report.html)\n"
	"  -t TITLE         report title\n"
	"  --compare DIR... runs to compare against, e.g. the same samples with\n"
	"                   ins_fix=0; ends at '--' or the next option\n"
	"\n"
	"  example:\n"
	"    %s report -o PD1.html --compare result/*_OFF -- result/*_ON\n", prog, prog);
}
static bool ends_with(const string &s, const string &x){ return s.size()>=x.size() && s.compare(s.size()-x.size(),x.size(),x)==0; }
static bool load_dir(const string &d, vector<Sample> &out)
{
	struct stat stt;
	if(stat(d.c_str(),&stt)==0 && S_ISREG(stt.st_mode))
	{
		Sample s; if(read_stats(d,s)){ out.push_back(s); return true; }
		fprintf(stderr,"[report] '%s' is not a kmer_align report file\n",d.c_str()); return false;
	}
	DIR *dp = opendir(d.c_str());
	if(dp==NULL){ fprintf(stderr,"[report] cannot open '%s'\n",d.c_str()); return false; }
	vector<string> hits; struct dirent *e;
	while((e=readdir(dp))!=NULL){ string n=e->d_name; if(ends_with(n,"_report.tsv")) hits.push_back(n); }
	closedir(dp);
	sort(hits.begin(),hits.end());
	if(hits.empty()){ fprintf(stderr,"[report] no *_report.tsv in '%s' (run with report=1)\n",d.c_str()); return false; }
	for(auto &h: hits){ Sample s; if(read_stats(d+"/"+h,s)) out.push_back(s); }
	return true;
}
int report_main(int argc, char **argv)
{
	/* argv[0] == "report" */
	string out="kmer_align_report.html", title;
	vector<string> dirs, cmp;
	int mode=0;   /* 0 dirs, 1 compare */
	for(int i=1;i<argc;i++)
	{
		string a=argv[i];
		if(a=="-h"||a=="--help"){ report_usage("kmer_align"); return 0; }
		else if(a=="-o" && i+1<argc){ out=argv[++i]; mode=0; }
		else if(a=="-t" && i+1<argc){ title=argv[++i]; mode=0; }
		else if(a=="--compare"){ mode=1; }
		else if(a=="--"){ mode=0; }
		else if(!a.empty() && a[0]=='-'){ fprintf(stderr,"[report] unknown option '%s'\n\n",a.c_str()); report_usage("kmer_align"); return 1; }
		else (mode==1?cmp:dirs).push_back(a);
	}
	if(dirs.empty()){ report_usage("kmer_align"); return 1; }
	vector<Sample> S, C;
	for(auto &d: dirs) load_dir(d,S);
	for(auto &d: cmp)  load_dir(d,C);
	if(S.empty()){ fprintf(stderr,"[report] nothing to report\n"); return 1; }
	if(title.empty())
	{
		string g; for(auto &p: S[0].params) if(p.first=="gene") g=p.second;
		title = (g.empty()?string("kmer_align"):g) + " — " + fmt_i((int)S.size()) + " samples";
	}
	if(!write_html(S,C,title,out)) return 1;
	printf("report: %zu sample(s)%s -> %s\n", S.size(), C.empty()?"":(" vs "+fmt_i((int)C.size())+" compared").c_str(), out.c_str());
	return 0;
}
