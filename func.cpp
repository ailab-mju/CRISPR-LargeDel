#include "func.h"
#include "arg_parse.h"
#include "samout.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <list>
#include <unordered_map>
#include <vector>
#include <string.h>
using namespace std;
extern char g_ref_seq[100000];
extern SamOut *f_out,*f_out2,*f_out3;
extern Args arg;
/* ------------------------------------------------------------------ input
   kmer_align reads plain text with fgets().  Sequencing data almost always
   arrives gzip compressed, so rather than making the caller decompress
   everything first, a compressed path is opened through the corresponding
   decompressor.  Nothing else in the code has to change. */
FILE *ka_open(const char *path,int *is_pipe)
{
	const char *dec=NULL;
	size_t n=strlen(path);
	if(n>3 && strcmp(path+n-3,".gz")==0)       {dec="gzip -dc";}
	else if(n>4 && strcmp(path+n-4,".bz2")==0) {dec="bzip2 -dc";}
	else if(n>3 && strcmp(path+n-3,".xz")==0)  {dec="xz -dc";}
	if(dec==NULL){*is_pipe=0;return fopen(path,"r");}
#ifdef __EMSCRIPTEN__
	/* no subprocesses in the browser build; the web app hands over plain text */
	fprintf(stderr,"[error] '%s' is compressed; the browser build reads plain FASTQ only\n",path);
	*is_pipe=0; return NULL;
#else
	/* single-quote the path and escape any quote inside it */
	string cmd=string(dec)+" '";
	for(size_t i=0;i<n;i++)
	{
		if(path[i]=='\''){cmd+="'\\''";}
		else{cmd+=path[i];}
	}
	cmd+="'";
	FILE *f=popen(cmd.c_str(),"r");
	*is_pipe = (f!=NULL);
	return f;
#endif
}
void ka_close(FILE *f,int is_pipe)
{
	if(f==NULL){return;}
#ifdef __EMSCRIPTEN__
	(void)is_pipe; fclose(f);
#else
	if(is_pipe){pclose(f);}else{fclose(f);}
#endif
}

void parse_read(read &r)
{
	r.query_len = strlen(r.query);
	char_to_num(r.query,r.query_len);
}



int make_ref(string refname,int K,unordered_map<long long int,vector<int>> &kmer_ref)
{
	char t[200];
	t['A'] = 0,t['C'] = 1,t['G'] = 2,t['T'] = 3;
	FILE *f1;
	int i,l;
	int ref_pipe=0;
	long long int h=0,last_bit=(1LL<<(2*(K-1)));	
	f1 = ka_open(refname.c_str(),&ref_pipe);
	if(f1==NULL)
	{
		fprintf(stderr,"[error] cannot open reference '%s'\n",refname.c_str());
		exit(1);
	}
	/* Read the header, then concatenate every following sequence line.
	   The old code used two fscanf("%s") calls, so a FASTA wrapped at 60 or
	   80 columns - which is what almost every tool writes - silently became a
	   60 bp reference and the run produced a 0% alignment rate with no error.
	   Wrapped and single-line FASTA now both work; a second record is
	   ignored, as before. */
	{
		char line[8192];
		size_t used=0;
		const size_t cap=sizeof(g_ref_seq)-1;
		int seen_header=0,truncated=0;
		while(fgets(line,sizeof(line),f1)!=NULL)
		{
			size_t ll=strlen(line);
			while(ll>0 && (line[ll-1]=='\n'||line[ll-1]=='\r')){line[--ll]='\0';}
			if(ll==0){continue;}
			if(line[0]=='>')
			{
				if(seen_header){break;}   /* stop at the second record */
				seen_header=1; continue;
			}
			if(!seen_header){continue;}
			for(size_t k=0;k<ll;k++)
			{
				char c=line[k];
				if(c>='a'&&c<='z'){c=(char)(c-'a'+'A');}
				if(c!='A'&&c!='C'&&c!='G'&&c!='T'&&c!='N'){continue;}
				if(used>=cap){truncated=1;break;}
				g_ref_seq[used++]=c;
			}
			if(truncated){break;}
		}
		g_ref_seq[used]='\0';
		if(truncated)
		{
			fprintf(stderr,"[error] reference '%s' is longer than %lu bp\n",
				refname.c_str(),(unsigned long)cap);
			ka_close(f1,ref_pipe);
			exit(1);
		}
		if(!seen_header)
		{
			fprintf(stderr,"[error] '%s' is not a FASTA (no '>' header line)\n",refname.c_str());
			ka_close(f1,ref_pipe);
			exit(1);
		}
	}
	l=strlen(g_ref_seq);
	if(l<=0)
	{
		fprintf(stderr,"[error] reference '%s' has an empty sequence\n",refname.c_str());
		ka_close(f1,ref_pipe);
		exit(1);
	}
	fprintf(stderr,"[ref] %s: %d bp\n",refname.c_str(),l);
	arg.ref_len = l;
	for(i=0;i<=l;i++) 
	{
		if(i>=K)
		{
			kmer_ref[h].push_back(i-K);
			h-=last_bit*t[g_ref_seq[i-K]];
		}
		h*=4;
		h+=t[g_ref_seq[i]];	
	}
	return l;
}
void char_to_num(char *query,int l)
{
	char t[200];
	t['A'] = 0,t['C'] = 1,t['G'] = 2,t['T'] = 3,t['N'] = 3;
	for(int i=0;i<l;i++)
	{
		query[i] = t[query[i]];
	}
}
char* num_to_char(char *query,int l)
{
	char t[5];
	t[0] = 'A',t[1] = 'C',t[2] = 'G',t[3] = 'T';
	for(int i=0;i<l;i++)
	{
		query[i] = t[query[i]];
	}
	return query;
}
void reverse_comp(char *query,int l)
{
	char query2[300];
	char oper = 252;
	for(int i=0;i<l;i++)
	{
		query2[l-i-1] = (~query[i])^oper;
	}
	for(int i=0;i<l;i++)
	{
		query[i] = query2[i];
	}
}


/* One SAM record as text.  samout.cpp buckets these by position and emits
   them in order, so nothing here knows or cares whether the final file is a
   BAM or a SAM. */
/* An unaligned read carries make_cigar()'s {-1,'*'} sentinel, and its POS /
   PNEXT / TLEN are whatever the previous alignment attempt left behind. */
bool read_aligned(const read &r)
{
	if(r.map_f==-1){return false;}
	if(r.cigar_string.empty()){return false;}
	if(r.cigar_string[0].len < 0){return false;}
	return true;
}
string sam_record(read &r)
{
	char buf[512];
	string out;
	snprintf(buf,sizeof(buf),"%s",r.query_name+1);
	out = buf;
	if(!r.pat.empty()){out += "|" + r.pat;}
	snprintf(buf,sizeof(buf),"\t%d\t%s\t%d\t%d\t",r.sam_f,r.ref.c_str(),r.pos,r.map_q);
	out += buf;
	/* make_cigar() marks "this read did not align" with a single {-1,'*'}
	   entry.  Printing that verbatim produced "-1*", which is not a legal
	   CIGAR - htslib rejects the whole record.  Emit the SAM spelling. */
	if(r.cigar_string.empty() || r.cigar_string[0].len < 0)
	{
		out += "*";
	}
	else
	{
		for(size_t i=0;i<r.cigar_string.size();i++)
		{
			snprintf(buf,sizeof(buf),"%d%c",r.cigar_string[i].len,r.cigar_string[i].type);
			out += buf;
		}
	}
	snprintf(buf,sizeof(buf),"\t=\t%d\t%d\t",r.p_next,r.tlen);
	out += buf;
	out += num_to_char(r.query,r.query_len);
	char_to_num(r.query,r.query_len);
	out += "\t";
	out += r.query_qual;
	/* empty slots used to be printed anyway, which left a bare \t\t in the
	   record and made the line unparseable */
	for(size_t i=0;i<r.tags.size();i++)
	{
		if(!r.tags[i].empty()){out += "\t" + r.tags[i];}
	}
	if(!r.pat.empty()){out += "\tXV:Z:" + r.pat;}
	return out;
}
void make_sam(read &r1, read &r2,int flag,string gene_name)
{
	r1.sam_f=r2.sam_f=3;
	r1.sam_f+=64;
	r2.sam_f+=128;
	r1.sam_f+=16+(!flag*16);
	r2.sam_f+=16+((flag)*16);
	r1.map_q=60;
	r2.map_q=60;
	r1.pos = r1.ref_st;
	r2.pos = r2.ref_st;
	r1.p_next = r2.pos;
	r2.p_next = r1.pos;
	r1.dist = min(r1.dist,r2.dist);
	r2.dist = min(r1.dist,r2.dist);
	if(r1.pos>=r2.pos)
	{
		r1.tlen = -(abs(r1.pos-r2.pos)+r1.query_len);
		r2.tlen = (abs(r1.pos-r2.pos)+r1.query_len);
	}
	else
	{
		r1.tlen = (abs(r1.pos-r2.pos)+r2.query_len);
		r2.tlen = -(abs(r1.pos-r2.pos)+r2.query_len);
	}
	r1.ref=(gene_name);
	r2.ref=(gene_name);
	r1.del_size = max((r1.ref_ed-r1.ref_st),(r2.ref_ed-r2.ref_st));
	r2.del_size = max((r1.ref_ed-r1.ref_st),(r2.ref_ed-r2.ref_st));
	r1.tags[2]="XD:i:"+to_string(r1.dist);
	r2.tags[2]="XD:i:"+to_string(r2.dist);
	r1.tags[3]="XL:i:"+to_string(r1.del_size);
	r2.tags[3]="XL:i:"+to_string(r2.del_size);
}
string rgb(string hex)
{
	char t[20];
	strcpy(t,hex.c_str());
	for(int i=1;i<=6;i++)
	{
		if(t[i]<='9')
		{
			t[i]-='0';
		}
		else
		{
			t[i]-='A';
			t[i]+=10;
		}
	}
	int r = t[2]+t[1]*16;
	int g = t[4]+t[3]*16;
	int b = t[6]+t[5]*16;
	string res=to_string(r)+","+to_string(g)+","+to_string(b);
	return res;
}
/* Geometry half of long_inner_check(): 1 when that function would record this
   pair as a read-pair large deletion, with no side effects.  Used by
   --ld_confirm to decide whether the pair has to be re-aligned first. */
int mate_gap_call(const read &r1,const read &r2)
{
	int l_p2,r_p2;
	if(r1.ref_st<r2.ref_st){l_p2=r1.ref_ed;r_p2=r2.ref_st;}
	else{l_p2=r2.ref_ed;r_p2=r1.ref_st;}
	if(!(r_p2-l_p2>arg.min_gap)){return 0;}
	if(r_p2 < arg.targetst || arg.targeted < l_p2){return 0;}
	int cutp = arg.targetst+arg.cle;
	if(arg.cut_between && !(l_p2 < cutp && cutp < r_p2)){return 0;}
	int ddist=10000;
	ddist = min(ddist,abs(r1.ref_st-cutp));
	ddist = min(ddist,abs(r1.ref_ed-cutp));
	ddist = min(ddist,abs(r2.ref_st-cutp));
	ddist = min(ddist,abs(r2.ref_ed-cutp));
	return ddist<=arg.dist_cut;
}
int long_inner_check(read &r1,read &r2,vector<read> *split,result &res)
{
	int inner_d,l_p,r_p,l_p2,r_p2;
	if(r1.ref_st<r2.ref_st)
	{
		inner_d=r2.ref_ed-r1.ref_st;
		l_p = r1.ref_st;
		r_p = r2.ref_ed;
		l_p2 = r1.ref_ed;
		r_p2 = r2.ref_st;
	}
	else
	{
		r_p = r1.ref_ed;
		l_p = r2.ref_st;
		r_p2 = r1.ref_st;
		l_p2 = r2.ref_ed;
		inner_d=r1.ref_ed-r2.ref_st;
		int tt;
	}
	/* Read-pair category (2) of the paper: "all read pairs are mapped without
	   splitting and the cleavage site passes between reads".

	   arg.min_gap (400 by default) is NOT the paper's 100 bp large-deletion
	   definition.  The only thing a non-split pair tells us is the distance
	   between the two mates, and with a ~300 bp library that distance is a few
	   hundred bp on a perfectly normal WT fragment.  Measured on the stress
	   simulation, dropping the floor to 100 quadruples the false large deletion
	   rate (0.14 % -> 0.61 %) and costs 10 pp of sensitivity, because ordinary
	   WT pairs start entering this branch; on PRR5 (a line with NO real large
	   deletions) it takes the large-deletion rate from 8.8 % to 40.1 %.  Keep it
	   at roughly (max insert size) + (smallest deletion worth calling). */
	if(r_p2-l_p2>arg.min_gap)
	{
		if(r_p2 < arg.targetst || arg.targeted < l_p2){
			r1.overlap = 0;
			r2.overlap = 0;
			r1.map_f=r2.map_f=-1;return -1;
		}
		/* The test above only asks that the gap OVERLAP the protospacer, which
		   is weaker than what the paper says.  With a large cleavage_offset the
		   cut can sit outside a gap that still overlaps the target, and those
		   pairs were being called large deletions.  On PRR5 (cle=16) requiring
		   the cut to fall strictly between the mates removes 54 % of the
		   mut_only calls and costs nothing on the simulation. */
		if(arg.cut_between)
		{
			int cutp = arg.targetst+arg.cle;
			if(!(l_p2 < cutp && cutp < r_p2))
			{
				r1.overlap = 0;
				r2.overlap = 0;
				r1.map_f=r2.map_f=-1;return -1;
			}
		}
		int ddist=10000;
		ddist = min(ddist,abs(r1.ref_st-arg.targetst-arg.cle));
//		ddist = min(ddist,abs(r1.ref_st-arg.targeted));
		ddist = min(ddist,abs(r1.ref_ed-arg.targetst-arg.cle));
//		ddist = min(ddist,abs(r1.ref_ed-arg.targeted));
		ddist = min(ddist,abs(r2.ref_st-arg.targetst-arg.cle));
//		ddist = min(ddist,abs(r2.ref_st-arg.targeted));
		ddist = min(ddist,abs(r2.ref_ed-arg.targetst-arg.cle));
//		ddist = min(ddist,abs(r2.ref_ed-arg.targeted));
		r1.map_f=2;r2.map_f=2;
		r1.ref_st=r2.ref_st=l_p;r1.ref_ed=r2.ref_ed=r_p;
		r1.del_size=r2.del_size=r_p-l_p;
		r1.dist=r2.dist=ddist;
		if(ddist>arg.dist_cut){r1.map_f=r2.map_f=-1;return -1;}
		else{
			r1.tags[2]="XD:i:"+to_string(r1.dist);
			r2.tags[2]="XD:i:"+to_string(r2.dist);
			r1.tags[3]="XL:i:"+to_string(r1.del_size);
			r2.tags[3]="XL:i:"+to_string(r2.del_size);
			r1.pat=r2.pat="LDEL";
			split[0].push_back(r1);
			split[1].push_back(r2);
			res.long_d++;
			return inner_d;
		}
	}
	res.normal++;
	/* Single denominator.
	   The old test used the pair's OUTER coordinates (l_p..r_p), so two mates
	   sitting either side of the cut passed it even though neither of them
	   sequenced the cut - the variant hides in the unsequenced inner gap.
	   Those pairs then aligned perfectly and were counted as WT.  A pair is
	   only informative when one of the mates actually covers the cut window,
	   or ends in a clip right at it; everything else is counted separately as
	   `not covered' and stays out of every rate.
	   Large deletions are unaffected: they leave through the branch above. */
	{
	int cut = arg.targetst+arg.cle;
	int cov1 = (r1.ref_st<cut-13 && cut+13<r1.ref_ed) || r1.sclip;
	int cov2 = (r2.ref_st<cut-13 && cut+13<r2.ref_ed) || r2.sclip;
	if(!(cov1||cov2))
	{
		if(l_p<cut-13 && cut+13<r_p)
		{
			res.tg_nocov++;
			r1.pat=r2.pat="NOCOV";
		}
		else
		{
			r1.pat=r2.pat="OUT";
		}
	}
	else
	{
				res.tg_align_0++;
		if(f_out!=NULL)
		{
			if((r1.ins|r1.del|r1.miss|r1.sclip)|(r2.ins|r2.del|r2.miss|r2.sclip))
			{
				if(r1.ins|r2.ins)
				{
					res.tg_ins++;
					r1.pat=r2.pat="INS";
					clip_rescue_observe(r1,r2);
					samout_write(f_out2,r1);
					samout_write(f_out2,r2);
				}
				else if(r1.del|r2.del)
				{
					res.tg_del++;
					r1.pat=r2.pat="DEL";
					clip_rescue_observe(r1,r2);
					samout_write(f_out2,r1);
					samout_write(f_out2,r2);
				}
				/* Checked before the substitution case: a terminal clip at the
				   cut site is stronger evidence than the mismatch flag, and a
				   clipped read often carries both. */
				else if(r1.sclip|r2.sclip)
				{
					res.tg_clip++;
					r1.pat=r2.pat="CLIP";
					/* held back until the sample is complete - see clip rescue */
					if(f_out3!=NULL && !clip_rescue_defer(r1,r2))
					{
						samout_write(f_out3,r1);
						samout_write(f_out3,r2);
					}
				}
				else
				{
					res.tg_miss++;
					r1.pat=r2.pat="SUB";
					samout_write(f_out2,r1);
					samout_write(f_out2,r2);
				}
			}
			else
			{
				r1.pat=r2.pat="WT";
				samout_write(f_out,r1);
				samout_write(f_out,r2);
			}
		}
		else
		{
			r1.pat=r2.pat="WT";
		}
	}
	}
	if(!(r_p<arg.targetst || l_p > arg.targeted)){
				res.tg_align_1++;
			if((r1.ins|r1.del|r1.miss)|(r2.ins|r2.del|r2.miss))
			{
				if(r1.ins|r2.ins)
				{
					res.tg_ins1++;
				}
				else if(r1.del|r2.del)
				{
					res.tg_del1++;
				}
				else
				{
					res.tg_miss1++;
				}
			}

	}
	if((l_p<arg.targetst && arg.targeted < r_p)){
		res.tg_align_2++;
			if((r1.ins|r1.del|r1.miss)|(r2.ins|r2.del|r2.miss))
			{
				if(r1.ins|r2.ins)
				{
					res.tg_ins2++;
				}
				else if(r1.del|r2.del)
				{
					res.tg_del2++;
				}
				else
				{
					res.tg_miss2++;
				}
			}
		}
	/* FIX: the original fell off the end of this int-returning function.
	   The caller (align_pair) branches on the result:
	       res == -1  -> discard the pair
	       res >= 400 -> force a full re-alignment
	   This path is the "normal pair, no large deletion" case, so it must
	   trigger neither.  0 is the only value satisfying both. */
	return 0;
}

/* ============================================================ clip rescue
   A small deletion close to the end of a read leaves only a few bases on one
   side of it.  side_align() finds those bases, but prune_short_ends() drops
   any end block shorter than MIN_match (13 bp) - deliberately, because a
   2-11 bp coincidental match a few bp downstream is how PRR5 grew phantom
   deletions out of a terminal insertion.  The price is that a real 1-10 bp
   deletion with a 5-12 bp anchor became a soft clip and was reported as a
   "terminal insertion": about 10 % of all small deletions on the noPCR
   simulation (HPRT1 / TRAC / BCL11A).

   Whether a short anchor is believable cannot be decided from the read alone,
   but it can be decided from the sample.  Every small indel that IS seen with
   full (>= 13 bp) anchors on both sides is collected per sample; a terminal
   clip pair is held back until the whole sample has been read, and a clipped
   mate is then re-assigned to a known allele when

     - the allele was seen in at least CLIP_RESCUE_MIN_SUPPORT read pairs with
       full anchors (any allele seen at all counts for the ambiguity test),
     - the read, laid onto that allele, matches it with no more mismatches than
       its aligned part already had against the reference (the clipped bases
       are explained exactly),
     - at least CLIP_RESCUE_MIN_ANCHOR read bases sit on the short side of the
       junction, and no known allele of the other kind (insertion vs deletion)
       explains the read as well - such a read is genuinely ambiguous.

   A phantom deletion is never supported by full-anchor reads (those show the
   insertion), so it cannot be rescued into; a real deletion always is, because
   most of its reads are not at a read edge.  Deletions of large_del_len and up
   are not collected (they go through the split-read path).  The rescued read
   gets a proper CIGAR (e.g. 7S143M -> 7M10D143M) and an XR:Z: tag naming the
   allele; --clip_rescue 0 turns all of this off. */
#define CLIP_RESCUE_MIN_SUPPORT 3
#define CLIP_RESCUE_MIN_ANCHOR  5
#define CLIP_RESCUE_MIN_CLIP    5

struct ClipAllele { char kind; int p; int d; string s; };   /* D: p,d   I: p,s */
struct ClipPending { read s1,s2,f1,f2; int has_full; };
static unordered_map<string,int>        g_cr_support;
static unordered_map<string,ClipAllele> g_cr_allele;
static vector<ClipPending>              g_cr_pending;
static int g_cr_waiting_full = 0;

static inline char cr_base(int c){ return (c>=0&&c<=3) ? "ACGT"[c] : 'N'; }

static string cr_norm_del(int p,int d,ClipAllele &a)
{
	/* left-normalise so equivalent placements in a repeat share one key */
	while(p>0 && g_ref_seq[p-1]==g_ref_seq[p+d-1]){p--;}
	a.kind='D';a.p=p;a.d=d;a.s="";
	return "D\t"+to_string(p)+"\t"+to_string(d);
}
static string cr_norm_ins(int p,string s,ClipAllele &a)
{
	while(p>0 && !s.empty() && g_ref_seq[p-1]==s.back()){ s = string(1,g_ref_seq[p-1]) + s.substr(0,s.size()-1); p--; }
	a.kind='I';a.p=p;a.d=0;a.s=s;
	return "I\t"+to_string(p)+"\t"+s;
}

/* small indels at the cut (+/-13) in one read, as normalised keys */
static void cr_events(const read &r, vector<pair<string,ClipAllele>> &ev)
{
	if(r.cigar_string.empty() || r.cigar_string[0].len<0){return;}
	int cut = arg.targetst+arg.cle;
	int p = r.ref_st-1, q = 0;
	for(size_t i=0;i<r.cigar_string.size();i++)
	{
		const cigar &c = r.cigar_string[i];
		if(c.type=='M'){p+=c.len;q+=c.len;}
		else if(c.type=='S'){q+=c.len;}
		else if(c.type=='D')
		{
			if(c.len<arg.large_del_len && !(p+c.len < cut-13 || cut+13 < p))
			{
				ClipAllele a; string k=cr_norm_del(p,c.len,a); ev.push_back(make_pair(k,a));
			}
			p+=c.len;
		}
		else if(c.type=='I')
		{
			if(!(p < cut-13 || cut+13 < p) && q+c.len<=r.query_len)
			{
				string s; for(int j=0;j<c.len;j++){s+=cr_base(r.query[q+j]);}
				ClipAllele a; string k=cr_norm_ins(p,s,a); ev.push_back(make_pair(k,a));
			}
			q+=c.len;
		}
	}
}

/* called for every pair written as a small INS / DEL */
void clip_rescue_observe(const read &r1,const read &r2)
{
	if(!arg.clip_rescue){return;}
	vector<pair<string,ClipAllele>> ev;
	cr_events(r1,ev);
	cr_events(r2,ev);
	sort(ev.begin(),ev.end(),[](const pair<string,ClipAllele>&x,const pair<string,ClipAllele>&y){return x.first<y.first;});
	for(size_t i=0;i<ev.size();i++)
	{
		if(i>0 && ev[i].first==ev[i-1].first){continue;}   /* once per pair */
		g_cr_support[ev[i].first]++;
		g_cr_allele[ev[i].first]=ev[i].second;
	}
}

/* called instead of writing a terminal-clip pair; returns 0 when rescue is off
   (the caller then writes the pair itself) */
int clip_rescue_defer(const read &r1,const read &r2)
{
	if(!arg.clip_rescue){return 0;}
	ClipPending p; p.s1=r1; p.s2=r2; p.has_full=0;
	g_cr_pending.push_back(p);
	g_cr_waiting_full=1;
	return 1;
}
int clip_rescue_waiting_full(){ return g_cr_waiting_full; }
void clip_rescue_attach_full(const read *r1,const read *r2)
{
	if(!g_cr_waiting_full || g_cr_pending.empty()){g_cr_waiting_full=0;return;}
	if(r1&&r2){ g_cr_pending.back().f1=*r1; g_cr_pending.back().f2=*r2; g_cr_pending.back().has_full=1; }
	g_cr_waiting_full=0;
}

static inline char cr_allele_base(const ClipAllele &a,int x)
{
	if(x<a.p){return g_ref_seq[x];}
	if(a.kind=='D'){return g_ref_seq[x+a.d];}
	int n=(int)a.s.size();
	if(x<a.p+n){return a.s[x-a.p];}
	return g_ref_seq[x-n];
}

/* Try to explain one clipped mate by a known allele.  On success the read is
   rewritten in place and 'D' or 'I' returned; otherwise 0. */
static char cr_try(read &r,string &tag)
{
	if(!r.sclip || r.cigar_string.size()!=2){return 0;}
	const cigar &c0=r.cigar_string[0], &c1=r.cigar_string[1];
	int lead;
	if(c0.type=='S'&&c1.type=='M'){lead=1;}
	else if(c0.type=='M'&&c1.type=='S'){lead=0;}
	else{return 0;}
	int L=r.query_len, clip=lead?c0.len:c1.len, m=lead?c1.len:c0.len;
	if(clip<CLIP_RESCUE_MIN_CLIP || clip+m!=L){return 0;}
	int x0 = r.ref_st-1;                         /* ref of the first M base */
	if(x0<0 || x0+m>arg.ref_len){return 0;}
	int mmR=0;
	for(int i=0;i<m;i++){ if(g_ref_seq[x0+i]!=cr_base(r.query[(lead?clip:0)+i])){mmR++;} }

	int best=-1,best_sup=0,seenD=0,seenI=0; int best_as=0;
	string best_key;
	for(auto &kv : g_cr_support)
	{
		int sup=kv.second;
		/* every observed allele takes part in the ambiguity test; only one
		   with CLIP_RESCUE_MIN_SUPPORT pairs behind it can be chosen */
		const ClipAllele &a=g_cr_allele[kv.first];
		int n=(a.kind=='I')?(int)a.s.size():0;
		int alen=(a.kind=='D')?arg.ref_len-a.d:arg.ref_len+n;
		int as;
		if(lead)
		{
			int xe=x0+m-1;                           /* anchor: last read base */
			if(a.kind=='D'){ if(xe<a.p+a.d){continue;} as=(xe-a.d)-(L-1); }
			else           { if(xe<a.p){continue;}     as=(xe+n)-(L-1); }
		}
		else
		{
			if(x0>=a.p){continue;}                   /* anchor: first read base */
			as=x0;
		}
		if(as<0 || as+L>alen){continue;}
		if(!(as < a.p+n && as+L > a.p)){continue;}   /* read must reach the event */
		int mmA=0;
		for(int i=0;i<L && mmA<=mmR;i++){ if(cr_allele_base(a,as+i)!=cr_base(r.query[i])){mmA++;} }
		if(mmA>mmR){continue;}
		/* this allele explains the read; it may only be CHOSEN when the read
		   holds the whole event with enough bases on both sides, but any
		   explaining allele counts for the ambiguity test below */
		if(a.kind=='D'){seenD=1;}else{seenI=1;}
		int left=a.p-as;                             /* read bases before the junction */
		int right=L-left-n;                          /* read bases after the event */
		if(sup<CLIP_RESCUE_MIN_SUPPORT){continue;}
		if(left<CLIP_RESCUE_MIN_ANCHOR || right<CLIP_RESCUE_MIN_ANCHOR){continue;}
		/* ties (two placements of the same insertion with equal support) go
		   to the smaller key: g_cr_support is an unordered_map, whose
		   iteration order differs between libstdc++ (native) and libc++
		   (WebAssembly), so "first one seen" made the two builds disagree */
		if(sup>best_sup || (sup==best_sup && best>=0 && kv.first<best_key))
		{ best_sup=sup; best=1; best_key=kv.first; best_as=as; }
	}
	/* A read that a known insertion explains just as well as a known deletion
	   is exactly the PRR5 situation (inserted bases that repeat downstream of
	   the cut): leave it a clip rather than pick a side. */
	if(best<0 || (seenD && seenI)){return 0;}
	const ClipAllele &a=g_cr_allele[best_key];
	int n=(a.kind=='I')?(int)a.s.size():0;
	int left=a.p-best_as;
	vector<cigar> nc;
	cigar t;
	t.type='M';t.len=left;nc.push_back(t);
	if(a.kind=='D'){t.type='D';t.len=a.d;}else{t.type='I';t.len=n;}
	nc.push_back(t);
	t.type='M';t.len=L-left-n;nc.push_back(t);
	r.cigar_string=nc;
	r.ref_st=best_as+1;
	r.ref_ed=r.ref_st+(L-n)+(a.kind=='D'?a.d:0);
	r.pos=r.ref_st;
	r.sclip=0;r.clip_len=0;r.clip_side=0;
	if(a.kind=='D'){r.del=1;}else{r.ins=1;}
	tag = (a.kind=='D') ? ("D"+to_string(a.d)+"@"+to_string(a.p+1)) : ("I"+to_string(n)+"@"+to_string(a.p+1));
	tag += "/n"+to_string(best_sup);
	r.tags[0]="XR:Z:clip_rescue:"+tag;
	return a.kind;
}

static void cr_fix_pair(read &r1,read &r2)
{
	r1.p_next=r2.pos; r2.p_next=r1.pos;
	int lo=min(r1.pos,r2.pos), hi=max(r1.ref_ed,r2.ref_ed);
	int span=hi-lo;
	if(r1.pos<=r2.pos){r1.tlen=span;r2.tlen=-span;}else{r1.tlen=-span;r2.tlen=span;}
}

/* end of a sample: decide every held-back pair, write it, fix the counters */
void clip_rescue_finish(result &res,SamOut *full)
{
	if(!arg.clip_rescue){return;}
	int n_del=0,n_ins=0,n_tot=(int)g_cr_pending.size();
	for(size_t i=0;i<g_cr_pending.size();i++)
	{
		ClipPending &p=g_cr_pending[i];
		vector<cigar> o1=p.s1.cigar_string,o2=p.s2.cigar_string;
		int op1=p.s1.pos,op2=p.s2.pos;
		string t1,t2;
		char k1=cr_try(p.s1,t1), k2=cr_try(p.s2,t2);
		string pat="CLIP";
		if(k1=='I'||k2=='I'){pat="INS";}
		else if(k1=='D'||k2=='D'){pat="DEL";}
		if(k1||k2){cr_fix_pair(p.s1,p.s2);}
		p.s1.pat=p.s2.pat=pat;
		if(pat=="CLIP")
		{
			if(f_out3!=NULL){samout_write(f_out3,p.s1);samout_write(f_out3,p.s2);}
		}
		else
		{
			res.tg_clip--;
			if(pat=="INS"){res.tg_ins++;n_ins++;}else{res.tg_del++;n_del++;}
			if(f_out2!=NULL){samout_write(f_out2,p.s1);samout_write(f_out2,p.s2);}
		}
		if(p.has_full && full!=NULL)
		{
			/* the whole-alignment copy may come from a later re-alignment; it
			   only takes over the rewrite when it is the same alignment */
			read *f[2]={&p.f1,&p.f2}; read *s[2]={&p.s1,&p.s2};
			vector<cigar> *o[2]={&o1,&o2}; int op[2]={op1,op2}; char k[2]={k1,k2};
			for(int j=0;j<2;j++)
			{
				f[j]->pat=pat;
				bool same = k[j] && f[j]->pos==op[j] && f[j]->cigar_string.size()==o[j]->size();
				for(size_t q=0; same && q<o[j]->size(); q++)
				{
					if(f[j]->cigar_string[q].len!=(*o[j])[q].len || f[j]->cigar_string[q].type!=(*o[j])[q].type){same=false;}
				}
				if(same)
				{
					f[j]->cigar_string=s[j]->cigar_string;f[j]->pos=s[j]->pos;
					f[j]->ref_st=s[j]->ref_st;f[j]->ref_ed=s[j]->ref_ed;f[j]->tags[0]=s[j]->tags[0];
				}
			}
			if(k1||k2){cr_fix_pair(p.f1,p.f2);}
			samout_write(full,p.f1);samout_write(full,p.f2);
		}
	}
	if(n_tot)
	{
		printf(" clip rescue: %d of %d terminal-clip pairs matched a small indel seen with full anchors (%d deletion, %d insertion)\n",
		       n_del+n_ins,n_tot,n_del,n_ins);
	}
	g_cr_pending.clear();g_cr_support.clear();g_cr_allele.clear();g_cr_waiting_full=0;
}
