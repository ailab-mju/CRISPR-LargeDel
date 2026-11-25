#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>
#include "func.h"
#include "arg_parse.h"
#include "samout.h"
#include "report.h"

#ifdef USE_HTSLIB
#include <htslib/sam.h>
#include <htslib/kstring.h>
#endif

using namespace std;

extern Args arg;

static int g_fmt = OUT_BAM;
void samout_set_format(int fmt){ g_fmt = fmt; }
int  samout_format(void){ return g_fmt; }
bool samout_have_bam(void)
{
#ifdef USE_HTSLIB
	return true;
#else
	return false;
#endif
}

/* 256 buckets over the amplicon.  Reads tile the whole reference in the
   whole-alignment output, so the buckets come out roughly even; the small
   per-site outputs pile into one or two buckets, which is harmless because
   they are small to begin with. */
#define NBUCKET 256
#define BUCKET_FLUSH (64 * 1024)

struct SamOut
{
	string base;
	int role;       /* report.h: which report tally this file feeds */
	string buf[NBUCKET];
	bool spilled[NBUCKET];
	long long n;
};

static string header_text(const char *so)
{
	string h = "@HD\tVN:1.6\tSO:";
	h += so; h += "\n";
	h += "@SQ\tSN:" + arg.genename + "\tLN:" + to_string(arg.ref_len) + "\n";
	h += "@PG\tID:kmer_align\tPN:CRISPR-LargeDel\tCL:" + g_cmdline + "\n";
	return h;
}

static string tmp_name(const string &base, int i)
{
	return base + ".ka_tmp." + to_string(i);
}

SamOut *samout_open(const string &base)
{
	if(base.empty()){return NULL;}
	SamOut *o = new SamOut();
	o->base = base;
	o->role = report_role_for(base);
	o->n = 0;
	for(int i=0;i<NBUCKET;i++){o->spilled[i]=false;}
	/* a stale temp file from a crashed run would be appended to, so clear */
	for(int i=0;i<NBUCKET;i++){remove(tmp_name(base,i).c_str());}
	return o;
}

static void spill(SamOut *o, int i)
{
	if(o->buf[i].empty()){return;}
	FILE *f = fopen(tmp_name(o->base,i).c_str(), o->spilled[i] ? "a" : "w");
	if(f==NULL)
	{
		fprintf(stderr,"[error] cannot write %s\n", tmp_name(o->base,i).c_str());
		o->buf[i].clear();
		return;
	}
	fwrite(o->buf[i].data(),1,o->buf[i].size(),f);
	fclose(f);
	o->spilled[i]=true;
	o->buf[i].clear();
}

void samout_write(SamOut *o, read &r)
{
	if(o==NULL){return;}
	string line = sam_record(r);
	report_observe(o->role, line);
	line += "\n";
	int pos = r.pos; if(pos<0){pos=0;}
	long long span = arg.ref_len>0 ? arg.ref_len : 1;
	int i = (int)((long long)pos * NBUCKET / (span + 1));
	if(i<0){i=0;} if(i>=NBUCKET){i=NBUCKET-1;}
	o->buf[i] += line;
	o->n++;
	if(o->buf[i].size() >= BUCKET_FLUSH){spill(o,i);}
}

static int line_pos(const char *s)
{
	/* field 4 of a SAM record */
	int tab = 0;
	for(const char *p=s; *p; p++)
	{
		if(*p=='\t'){ if(++tab==3){ return atoi(p+1); } }
	}
	return 0;
}

static void load_bucket(SamOut *o, int i, vector<string> &out)
{
	string all;
	if(o->spilled[i])
	{
		FILE *f = fopen(tmp_name(o->base,i).c_str(),"rb");
		if(f)
		{
			char chunk[1<<16]; size_t k;
			while((k=fread(chunk,1,sizeof(chunk),f))>0){all.append(chunk,k);}
			fclose(f);
		}
		remove(tmp_name(o->base,i).c_str());
	}
	all += o->buf[i];
	o->buf[i].clear();
	size_t st = 0;
	while(st < all.size())
	{
		size_t e = all.find('\n', st);
		if(e==string::npos){e=all.size();}
		if(e>st){out.push_back(all.substr(st,e-st));}
		st = e+1;
	}
	stable_sort(out.begin(), out.end(), [](const string &a, const string &b){
		return line_pos(a.c_str()) < line_pos(b.c_str());
	});
}

void samout_close(SamOut *o)
{
	if(o==NULL){return;}
	int fmt = g_fmt;
#ifndef USE_HTSLIB
	fmt = OUT_SAM;                     /* no htslib -> text only */
#endif
	FILE *fsam = NULL;
	if(fmt==OUT_SAM || fmt==OUT_BOTH)
	{
		fsam = fopen((o->base+".sam").c_str(),"w");
		if(fsam){string h=header_text("coordinate"); fwrite(h.data(),1,h.size(),fsam);}
	}
#ifdef USE_HTSLIB
	samFile *fbam = NULL; sam_hdr_t *hdr = NULL; bam1_t *brec = NULL;
	bool idx_on = false;
	/* htslib keeps the pointer it is handed, so this must outlive the call -
	   a temporary's c_str() dangles the moment sam_idx_init() returns. */
	string idxpath = o->base + ".bam.bai";
	long long werr = 0;
	if(fmt==OUT_BAM || fmt==OUT_BOTH)
	{
		string hdrtxt = header_text("coordinate");
		hdr = sam_hdr_parse(hdrtxt.size(), hdrtxt.c_str());
		fbam = sam_open((o->base+".bam").c_str(),"wb");
		if(fbam==NULL || hdr==NULL)
		{
			fprintf(stderr,"[error] cannot open %s.bam\n",o->base.c_str());
		}
		else
		{
			if(sam_hdr_write(fbam,hdr)<0){fprintf(stderr,"[error] header write failed\n");}
			/* min_shift 0 -> classic .bai (CSI would be min_shift 14) */
			int ir = sam_idx_init(fbam,hdr,0,idxpath.c_str());
			if(ir==0){idx_on=true;}
			else{fprintf(stderr,"[warn] sam_idx_init failed (%d) for %s.bam\n",ir,o->base.c_str());}
			brec = bam_init1();
		}
	}
#endif
	for(int i=0;i<NBUCKET;i++)
	{
		vector<string> lines;
		load_bucket(o,i,lines);
		for(size_t j=0;j<lines.size();j++)
		{
			if(fsam){fwrite(lines[j].data(),1,lines[j].size(),fsam); fputc('\n',fsam);}
#ifdef USE_HTSLIB
			if(fbam && brec)
			{
				kstring_t ks = {0,0,NULL};
				kputsn(lines[j].data(),lines[j].size(),&ks);
				if(sam_parse1(&ks,hdr,brec)<0)
				{
					fprintf(stderr,"[warn] unparseable record skipped in %s\n",o->base.c_str());
				}
				else if(sam_write1(fbam,hdr,brec)<0)
				{
					/* a failed write here means the BAM is truncated - disk
					   full, or a bad record - and must not pass silently */
					fprintf(stderr,"[error] BAM write failed for %s (disk full?)\n",o->base.c_str());
					werr++;
				}
				free(ks.s);
			}
#endif
		}
	}
	if(fsam){fclose(fsam);}
#ifdef USE_HTSLIB
	if(fbam)
	{
		if(idx_on && sam_idx_save(fbam)<0){fprintf(stderr,"[warn] index write failed for %s.bam\n",o->base.c_str());}
		if(sam_close(fbam)<0){fprintf(stderr,"[error] closing %s.bam failed\n",o->base.c_str());}
		else{fbam=NULL;}
		if(werr){fprintf(stderr,"[error] %lld record(s) lost writing %s.bam\n",werr,o->base.c_str());}
	}
	if(fbam){sam_close(fbam);}
	if(brec){bam_destroy1(brec);}
	if(hdr){sam_hdr_destroy(hdr);}
#endif
	delete o;
}
