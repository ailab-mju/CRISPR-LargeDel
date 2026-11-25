#include<stack>
#include<stdio.h>
#include<limits>
#include<list>
#include<unordered_map>
#include<algorithm>
#include<string>
#include<string.h>
#include<vector>
#include<stdlib.h>
#include<functional>
#include "arg_parse.h"
#include "align.h"
#include "align_wrap.h"
#define MIN_match 13
/* Large-deletion length threshold for a split read: a D of at least this many
   bp makes the read a large-deletion candidate (map_f=2) instead of a small
   indel.  Was a compile-time #define (50); now set at run time from argument
   #22 / --large_del_len (default 50 as in the original; 100 = the paper's
   definition, the default during 2026-10).  MIN_skip stays
   the internal name. */
#define MIN_skip (arg.large_del_len)
/* A terminal insertion cannot become an 'I' in the CIGAR: make_cigar() only
   emits 'I' when two aligned blocks overlap in reference coordinates, which
   needs aligned sequence on BOTH sides of the inserted bases.  At the edge of
   a read the inserted bases end up in a soft clip instead and the pair used to
   be reported as WT.  A clip of at least MIN_CLIP bases whose boundary sits
   within CLIP_WIN of the cut site is now flagged as insertion evidence. */
#define MIN_CLIP 5
#define CLIP_WIN 13
using namespace std;
extern char g_ref_seq[100000],view_q[1000];
extern Args arg;
extern int skip_read;
void side_align(vector<data> &res,read &r)
{
	if(res.size()==0){return;}
	data tmp;
	char t[5];
	t[0] = 'A',t[1] = 'C',t[2] = 'G',t[3] = 'T';
	int i,l=res.size(),last=-1,p,miss=0,min,mini,j,sw=0;
	if(arg.targetst+arg.cle-20 <= res[0].st && res[0].st<=arg.targetst+arg.cle+20)
	{
		for(i=res[0].st2;i<=res[0].ed2;i++)
		{
			if(i-res[0].st2>20){break;}
			p=res[0].st+(i-res[0].st2);
			if(g_ref_seq[p]!=t[r.query[i]])
			{
				miss++;
				last=i;
			}
		}
		min=miss;
		miss=0;
		/* FIX: 'last' stays -1 when the first scan found no mismatch.
		   The original only tested this AFTER the refinement loops, so
		   r.query[-1] / g_ref_seq[...-1] were read out of bounds.
		   No mismatch near the cut site => nothing to re-align. */
		if(last==-1){return;}
		for(i=-1;i>=-30;i--)
		{
			miss=0;
			for(j=0;j<=last;j++)
			{
				if(t[r.query[j]]!=g_ref_seq[res[0].st+i+j])
				{
					miss++;
				}
			}
			if(miss<min)
			{
				sw=1;
				min = miss;
				mini = i;
				tmp.st2 = 0;
				tmp.ed2 = last;
				tmp.miss=min;
				tmp.st = res[0].st+mini;
				tmp.ed = res[0].st+mini+last;
			}
		}
		if(sw){
		res[0].st = res[0].st+(tmp.ed2+1-res[0].st2);
		res[0].st2 = tmp.ed2+1;
		res.push_back(tmp);
		for(i=res.size()-1;i>=1;i--)
		{
			res[i] = res[i-1];
			}
		res[0] = tmp;
		}
//		printf("%d %d\n",min,mini);
//		printf("%d %d %d %d\n",tmp.st,tmp.ed,tmp.st2,tmp.ed2);

	}
	else if(arg.targetst+arg.cle-20 <= res[l-1].ed && res[l-1].ed<=arg.targetst+arg.cle+20)
	{
		for(i=res[l-1].ed2;i>=res[l-1].st2;i--)
		{
			if(res[l-1].ed2-i>20){break;}
			p=(res[l-1].ed-(res[l-1].ed2-i));
			if(g_ref_seq[p]!=t[r.query[i]])
			{
				miss++;
				last=i;
			}
		}
		min=miss;
		miss=0;
		/* FIX: 'last' stays -1 when the first scan found no mismatch.
		   The original only tested this AFTER the refinement loops, so
		   r.query[-1] / g_ref_seq[...-1] were read out of bounds.
		   No mismatch near the cut site => nothing to re-align. */
		if(last==-1){return;}
		for(i=1;i<=30;i++)
		{
			miss=0;
			for(j=last;j<=res[l-1].ed2;j++)
			{
				if(t[r.query[j]]!=g_ref_seq[res[l-1].st+i+j])
				{
					miss++;
				}
			}
			if(miss<min)
			{
				sw=1;
				min = miss;
				mini = i;
				tmp.st2 = last;
				tmp.ed2 = res[l-1].ed2;
				tmp.miss=min;
				tmp.st = res[l-1].ed+mini-(res[l-1].ed2-last);
				tmp.ed = res[l-1].ed+mini;
			}
		}
		if(sw){
		res[l-1].ed = res[l-1].ed-(res[l-1].ed2-last)-1;
		res[l-1].ed2 = last-1;
		res.push_back(tmp);}

//		printf("%d %d\n",min,mini);
//		printf("%d %d %d %d\n",tmp.st,tmp.ed,tmp.st2,tmp.ed2);

	}
	else
	{
		return;
	}
	if(last==-1){return;}

//	printf("%d\n",last);

}
/* A reference gap is only believable when both flanking blocks are real
   alignments.  side_align() re-places the first/last block up to 30 bp away to
   shave off mismatches and pushes the result without any length check, so a
   2-3 bp coincidental anchor downstream of the cut is enough to fabricate a
   large 'D'.  That is how a terminal insertion whose novel sequence happens to
   match a little way downstream gets reported as a deletion.  Blocks shorter
   than MIN_match at either end are dropped here; make_cigar() then turns the
   leftover into a soft clip, which classify_clip() picks up as an insertion.
   Seeded chunks always pass chunk_align()'s own MIN_match test, so nothing but
   side_align()'s output is affected. */
static void prune_short_ends(vector<data> &res)
{
	while(res.size()>1 && res.back().ed2-res.back().st2+1 < MIN_match){res.pop_back();}
	while(res.size()>1 && res[0].ed2-res[0].st2+1 < MIN_match){res.erase(res.begin());}
}
/* Terminal soft clip next to the cut site.  A soft clip consumes no reference,
   so the boundary is ref_st for a leading clip and ref_ed for a trailing one. */
static void classify_clip(read &r)
{
	int cut = arg.targetst+arg.cle;
	r.sclip=0;r.clip_len=0;r.clip_side=0;
	if(r.cigar_string.empty()){return;}
	cigar first = r.cigar_string.front();
	cigar last  = r.cigar_string.back();
	if(first.type=='S' && first.len>=MIN_CLIP && abs(r.ref_st-cut)<=CLIP_WIN)
	{
		r.sclip=1;r.clip_len=first.len;r.clip_side='L';
	}
	if(last.type=='S' && r.cigar_string.size()>1 && last.len>=MIN_CLIP && abs(r.ref_ed-cut)<=CLIP_WIN)
	{
		if(last.len>r.clip_len){r.clip_len=last.len;r.clip_side='R';}
		r.sclip=1;
	}
}
void make_cigar(vector<data> &res, int K, vector<cigar> &cigar_string,int query_len);
/* mismatches of a read against the reference where its CIGAR puts it */
int placed_mm(const read &r)
{
	static const char tb[4]={'A','C','G','T'};
	int p=r.ref_st-1,q=0,mm=0;
	for(size_t i=0;i<r.cigar_string.size();i++)
	{
		const cigar &c=r.cigar_string[i];
		if(c.len<0){return 1<<20;}
		if(c.type=='M')
		{
			for(int k=0;k<c.len;k++,p++,q++)
			{
				if(p<0||p>=arg.ref_len||q>=r.query_len){mm++;continue;}
				if(g_ref_seq[p]!=tb[(int)r.query[q]&3]){mm++;}
			}
		}
		else if(c.type=='D'){p+=c.len;}
		else if(c.type=='I'||c.type=='S'){q+=c.len;}
	}
	return mm;
}
/* --split_wt_margin.  A split read's evidence for a large deletion is its short
   side, often only 13-20 bp read across the gap.  In a repeat family the same bases sit in
   other copies, and one sequencing error that breaks the k-mers on the straight-on
   diagonal lets the chain jump to a copy kb away.  So the short side is also read
   straight on from the block next to it (+-3 bp, for a small indel at the join):
   if the far placement saves no more than split_wt_margin mismatches, the read is
   aligned straight on instead and classified like any other read.  Returns 1 when
   it did that (r is then fully classified). */
static int split_reads_on(read &r,vector<data> &res)
{
	int g=-1;
	for(size_t i=1;i<res.size();i++){ if(res[i].st-res[i-1].ed-1>=MIN_skip){g=(int)i;break;} }
	if(g<0){return 0;}
	int left=res[g-1].ed2+1, right=r.query_len-res[g].st2;
	int tail_right = right<=left;
	int b0 = tail_right?g:0, b1 = tail_right?(int)res.size():g;   /* short-side blocks [b0,b1) */
	int diag = tail_right ? res[g-1].st-res[g-1].st2 : res[g].st-res[g].st2;
	static const char tb[4]={'A','C','G','T'};
	/* far=1: mismatches where the chain put the short side; far=0: read on, shifted */
	auto mm_at=[&](int shift,int far)->int{
		int mm=0;
		for(int b=b0;b<b1;b++)
			for(int q=res[b].st2;q<=res[b].ed2;q++)
			{
				int p = far ? res[b].st+(q-res[b].st2) : q+diag+shift;
				if(p<0||p>=arg.ref_len){return 1<<20;}
				int c=(int)r.query[q];
				if(c<0||c>3){continue;}
				if(g_ref_seq[p]!=tb[c]){mm++;}
			}
		return mm;
	};
	int mm_far = mm_at(0,1);
	int best=1<<20,best_s=0;
	static const int order[7]={0,-1,1,-2,2,-3,3};
	for(int k=0;k<7;k++){ int m=mm_at(order[k],0); if(m<best){best=m;best_s=order[k];} }
	if(best-mm_far > arg.split_wt_margin){return 0;}
	/* rebuild the blocks with the short side on that diagonal */
	vector<data> nr;
	data t; t.miss=0;
	if(tail_right)
	{
		for(int b=0;b<g;b++){nr.push_back(res[b]);}
		int q0=res[g-1].ed2+1, q1=res.back().ed2;
		if(best_s==0){nr.back().ed2=q1;nr.back().ed=q1+diag;}
		else{t.st2=q0;t.ed2=q1;t.st=q0+diag+best_s;t.ed=q1+diag+best_s;nr.push_back(t);}
	}
	else
	{
		int q0=res[0].st2, q1=res[g].st2-1;
		if(best_s==0){t=res[g];t.st2=q0;t.st=q0+diag;nr.push_back(t);}
		else{t.st2=q0;t.ed2=q1;t.st=q0+diag+best_s;t.ed=q1+diag+best_s;nr.push_back(t);nr.push_back(res[g]);}
		for(size_t b=g+1;b<res.size();b++){nr.push_back(res[b]);}
	}
	res=nr;
	r.cigar_string.clear();
	make_cigar(res,0,r.cigar_string,r.query_len);
	/* the mismatch flags align_query() set were for the old blocks */
	int cut=arg.targetst+arg.cle,miss=0,miss_tg=0,match_tg=0;
	for(size_t b=0;b<res.size();b++)
		for(int q=res[b].st2;q<=res[b].ed2;q++)
		{
			int p=res[b].st+(q-res[b].st2);
			if(p<0||p>=arg.ref_len){continue;}
			if(tb[(int)r.query[q]&3]!=g_ref_seq[p]){miss++;if(cut-5<=p&&p<=cut+5){miss_tg++;}}
			else if(cut-15<=p&&p<=cut+15){match_tg++;}
		}
	r.tg_miss=miss_tg;r.miss=(miss_tg>0);r.tg_match=match_tg;
	r.tags[1]="XM:i:"+to_string(miss)+"\tXT:i:"+to_string(miss_tg)+"\tXW:Z:read_on:"+to_string(mm_far)+"/"+to_string(best);
	read_classify(r,res);
	return 1;
}
void read_classify(read &r,vector<data> &res)
{
	int i;
	r.ins=r.del=0;
	r.sclip=0;r.clip_len=0;r.clip_side=0;
	r.tags[0]="";
	int m=0,pos=0,dist=50000,skip=-1,n_del=0;
	if(r.cigar_string.empty() || res.empty()|| r.cigar_string[0].len==-1){r.map_f=-1;r.ref_st = 0;r.ref_ed = 0; return;}
	vector<int> match;
	cigar now;
	pos = res[0].st+1;
	if(pos<40){r.map_f=-1;r.ref_st=0;r.ref_ed=0;return;}
	r.ref_st = pos;
	for(i=0;i<r.cigar_string.size();i++)
	{
		dist = min(dist, abs(pos-arg.targetst-arg.cle));
		now = r.cigar_string[i];
		switch(now.type){
			case 'M':
				m+=now.len;
				pos+=now.len;
				match.push_back(now.len);
				break;
			case 'D':
				if(now.len>=MIN_skip){n_del++;}
				if(!(pos+now.len<arg.targetst+arg.cle-13|| arg.targetst+arg.cle+13<pos)){r.del=1;}
				skip = max(skip,now.len);
				pos+=now.len;
				break;
			case 'I':
				/* An insertion consumes QUERY, not reference.  The pre-2026
				   code moved the reference cursor backwards here, which was
				   harmless for a 1-3 bp I and corrupts ref_ed badly once the
				   I gets long.  ins_fix=0 restores the old behaviour. */
				if(arg.ins_fix)
				{
					if(!(pos<arg.targetst+arg.cle-13|| arg.targetst+arg.cle+13<pos)){r.ins=1;}
				}
				else
				{
					if(!(pos<arg.targetst+arg.cle-13|| arg.targetst+arg.cle+13<pos-now.len)){r.ins=1;}
					pos-=now.len;
				}
				break;
			default:
				break;
			}

	}
	r.ref_ed = pos;
	classify_clip(r);
	if(r.sclip)
	{
		r.tags[0]="XC:i:"+to_string(r.clip_len)+"\tXS:A:"+string(1,r.clip_side);
	}
	dist = min(dist, abs(pos-arg.targetst));
	dist = min(dist, abs(pos-arg.targeted));
	r.dist = 50000;
	sort(match.begin(),match.end(), greater<int>() );
	if(n_del>=2){r.map_f=-1;return;}
	if(r.ref_ed < arg.targetst || arg.targeted < r.ref_st){
		r.overlap = 0;
	}
	else{
		if(dist<arg.dist_cut){
			r.overlap = 1;
		}
		else
		{
			r.overlap = 0;
		}
	}
	/* Anchor quality.
	   side_align() re-places the first/last aligned block anywhere within 30 bp
	   that lowers its mismatch count, and prune_short_ends() only requires the
	   survivor to be MIN_match (13 bp) long.  A short block full of mismatches
	   therefore fabricates a split-read deletion out of nothing - that is how
	   PRR5 grew a phantom 26D out of 11 bp of coincidental homology.

	   Score every aligned block as (length - 4 x mismatches) and keep the worst
	   one.  The long body block scores in the hundreds, so in practice this is a
	   test on the anchor: at the default 13 a 13 bp anchor must be perfect and a
	   17 bp anchor may carry one mismatch.  Cost on the stress simulation is
	   0.3 pp of small-deletion sensitivity (large deletions actually gain 1 pp);
	   on PRR5 it removes 40 % of the false large deletions and lifts the median
	   surviving anchor from 28 bp to 44 bp.  The score is always reported as
	   XA:i: so it can be filtered on afterwards. */
	int anchor_score = 1000;
	{
		char tb[4]={'A','C','G','T'};
		for(size_t bi=0;bi<res.size();bi++)
		{
			int L=res[bi].ed2-res[bi].st2+1;
			if(L<=0){continue;}
			int mm=0;
			for(int qi=res[bi].st2;qi<=res[bi].ed2;qi++)
			{
				int p=res[bi].st+(qi-res[bi].st2);
				if(p<0||p>=arg.ref_len){continue;}
				int c=(int)r.query[qi];
				if(c<0||c>3){continue;}
				if(g_ref_seq[p]!=tb[c]){mm++;}
			}
			int sc=L-4*mm;
			if(sc<anchor_score){anchor_score=sc;}
		}
	}
	if(r.tags.size()>4){r.tags[4]="XA:i:"+to_string(anchor_score);}
	if(match.size()==1){r.map_f = 1;return;}
	if(arg.anchor_min>0 && skip>=MIN_skip && anchor_score<arg.anchor_min)
	{
		r.map_f=-1;r.ref_st=0;r.ref_ed=0;return;
	}
	if(skip>=MIN_skip && match[1] >= MIN_match && r.overlap ==1)
	{
		r.dist = dist;
		if(abs(arg.ref_len-r.ref_ed)<40)
		{
			r.map_f=-1;r.ref_st=0;r.ref_ed=0;return;
		}
		if(arg.split_wt_margin>=0 && split_reads_on(r,res)){return;}
		r.map_f=2;
		return;
	}
	else if(skip>=MIN_skip && match[1] >=MIN_match && r.overlap==0)
	{
		r.map_f=-1;
		return;
	}
	r.map_f=1;
}
int pre_chunk(char *query, unordered_map<long long int,vector<int>> &kmer_ref2,int l)
{
	int flag=0,i;
	int K=17;
	long long int h=0,last_bit=(1LL<<(2*(K-1))),h1=0,hh1=0,hh2=0;
	for(i=0;i<K;i++)
	{
		h*=4;
		h+=query[i];
	}
	for(i=l-K;i<l;i++)
	{
		h1*=4;
		h1+=query[i];
	}
	hh1 = h;
	hh2 = h1;
	if(kmer_ref2[h].size()!=0&&kmer_ref2[h][0]>arg.targeted+50){flag=kmer_ref2[h][0]; return flag;}
	if(kmer_ref2[h1].size()!=0&&kmer_ref2[h1][0]+K<arg.targetst-50){flag=kmer_ref2[h1][0]-l+K; return flag;}
	if(kmer_ref2[h].size()==0||kmer_ref2[h1].size()==0){return -1;}
	if(abs(l-(kmer_ref2[h1][0]-kmer_ref2[h][0]+K)) <= 3)
	{
        flag=-1;
        //flag=kmer_ref2[h][0];
		return flag;
	}
	else
	{
		flag=-1;
		return flag;
	}
}
void LIS(vector<data> &kmer_list,char *query, int K, unordered_map<long long int,vector<int>> &kmer_ref,int l,vector<int> &d, vector<int> &via)
{
	d.clear();
	via.clear();
	int i;
	long long int h=0,last_bit=(1LL<<(2*(K-1)));
	vector<data> res;
	data tmp;
	for(i=0;i<=l;i++)
	{
		if(i>=K)
		{
			for(int j=0;j<kmer_ref[h].size();j++)
			{
				tmp.st = kmer_ref[h][j],tmp.ed = kmer_ref[h][j]+K-1; //ref 위치
				tmp.st2 = i-K,tmp.ed2 = i-1; //query 위치
				kmer_list.push_back(tmp);
			}
			h-=last_bit*query[i-K];
		}
		h*=4;
		h+=query[i];
	}
	sort(kmer_list.begin(),kmer_list.end(),dataComparator());
	d.push_back(-1);
	for(i=0;i<kmer_list.size();i++)
	{
		if(d.back()<kmer_list[i].st2)
		{
			d.push_back(kmer_list[i].st2);
			via.push_back(d.size()-2);
		}
		else{
			auto it = lower_bound(d.begin(),d.end(),kmer_list[i].st2);
			*it = kmer_list[i].st2;
			via.push_back((int)(it-d.begin()));
		}
	}
}

void make_cigar(vector<data> &res, int K, vector<cigar> &cigar_string,int query_len)
{
	if(res.empty()){return;}
	cigar tmp;
	int match=0,r_match=0,skip=0;
	if(res.size()==0){tmp.len = -1; tmp.type='*'; cigar_string.push_back(tmp);return;}
	if(res[0].st2>0){ tmp.len = res[0].st2; tmp.type='S'; match = tmp.len;cigar_string.push_back(tmp);}
	tmp.len = res[0].ed2-res[0].st2+1;
	if(tmp.len>0){
	match+=tmp.len;
	tmp.type= 'M';
	r_match+=tmp.len;
	cigar_string.push_back(tmp);
	}
	int ins_total=0;
	for(int i=1;i<res.size();i++)
	{
		tmp.len=0;
		tmp.type='Z';

		/* A gap in QUERY coordinates is inserted sequence; an overlap in
		   REFERENCE coordinates is the pre-2026 code's own way of spotting an
		   insertion.  Both can be true of one junction, so they are summed
		   into a single I - two adjacent I ops are not a legal CIGAR and
		   htslib rejects the record. */
		int qgap = 0;
		if(arg.ins_fix)
		{
			qgap = res[i].st2-res[i-1].ed2-1;
			if(qgap<0){qgap=0;}
		}
		int ovl = 0;
		if(res[i].st<=res[i-1].ed)
		{
			ovl = res[i-1].ed-res[i].st+1;
			res[i].st2+=ovl;
			res[i].st+=ovl;
		}
		if(qgap+ovl>0)
		{
			tmp.len = qgap+ovl;
			tmp.type='I';
			match += tmp.len;
			ins_total += qgap;
		}
		else if(res[i].st>res[i-1].ed+1)
		{
			tmp.len = res[i].st-(res[i-1].ed+1);
			tmp.type='D';
		}
		if(tmp.len!=0)
		{
			cigar_string.push_back(tmp);
		}
		tmp.len = res[i].ed2-res[i].st2+1;
		tmp.type='M';
		r_match+=tmp.len;
		match+=tmp.len;
		cigar_string.push_back(tmp);
	}
	if(query_len > match)
	{
		tmp.len = query_len-match;
		tmp.type = 'S';
		cigar_string.push_back(tmp);
	}
	/* inserted bases are query the reference cannot match, so they must not
	   count against the coverage floor */
	if(r_match < 0.8*(query_len-ins_total))
	{
		cigar_string.clear();
		tmp.len = -1;
		tmp.type='*';
		cigar_string.push_back(tmp);
		res.clear();
	}

	return;
}
void chunk_align(vector<data> &lis_res,vector<data> &lis_res2)
{
	data tmp;
	lis_res2.clear();
	tmp.st=-1;tmp.ed=-1;
	tmp.miss=0;
	for(int i=0;i<lis_res.size();i++)
	{
		if(lis_res[i].st2 >= lis_res[i].ed2){continue;}
		if(tmp.st==-1)
		{
			tmp = lis_res[i];
		}
		else if(tmp.ed-lis_res[i].ed == tmp.ed2 - lis_res[i].ed2)// && tmp.miss<4)
		{
//			tmp.miss++;
			tmp.ed = lis_res[i].ed;
			tmp.ed2 = lis_res[i].ed2;
		}
		else if(tmp.st <= lis_res[i].st && lis_res[i].ed <= tmp.ed)
		{
			continue;
		}
		else
		{
			if(tmp.ed-tmp.st+1>=MIN_match&& tmp.st2<tmp.ed2){lis_res2.push_back(tmp);}
			tmp = lis_res[i];
		}
	}
	if(tmp.ed-tmp.st+1>=MIN_match&& tmp.st2<tmp.ed2){
	lis_res2.push_back(tmp);
	}
}

void align_query(int K, unordered_map<long long int,vector<int>> &kmer_ref, read &r,unordered_map<long long int,vector<int>> &kmer_ref2, int flag2)
{
	int i,pre_flag=-1;
	list<data> chunk;
	vector<data> res,kmer_list;
	vector<int> d,via;
	r.cigar_string.clear();
	r.miss=r.ins=r.del=0;
	r.sclip=0;r.clip_len=0;r.clip_side=0;
	r.tags[0]="";
	pre_flag=pre_chunk(r.query,kmer_ref2,r.query_len);
    //pre_flag=-1;
    if(pre_flag!=-1 && flag2==0)
	{
		data tt;
		cigar tt2;
		tt.st = pre_flag;
		tt.st2=0;
		tt.ed2=r.query_len-1;
		tt.ed = pre_flag+r.query_len-1;
		res.push_back(tt);
		tt2.type='M';
		tt2.len=r.query_len;
		r.cigar_string.push_back(tt2);
		skip_read++;
		read_classify(r,res);
	}
	else if(K>0)// &&strcmp(r.query_name+1,"MN00350:100:000H2TV2H:1:13101:8940:17307")==0)
	{
/*		if(strcmp(r.query_name+1,"MN00254:225:000H2JMTW:1:22104:21400:10375")==0){
			printf("asdfasdfasdf");
		}*/
		LIS(kmer_list,r.query,K,kmer_ref,r.query_len,d,via);
		stack<int> st;
		int l=d.size()-2;
		for(i=via.size()-1;i>=0;i--)
		{
			if(l==0)
			{
				st.push(i);
				break;
			}
			if(l==via[i])
			{
				st.push(i);
				l--;
			}
		}
		data tmp;
		tmp.st = -1,tmp.ed=-1;
		vector<data> lis_res,lis_res2;
		while(!st.empty())
		{
			if(tmp.st==-1)
			{
				tmp.st = kmer_list[st.top()].st;
				tmp.ed = kmer_list[st.top()].ed;
				tmp.st2 = kmer_list[st.top()].st2;
				tmp.ed2 = kmer_list[st.top()].ed2;
				tmp.miss=0;
			}
			else
			{
				if(tmp.ed+1==kmer_list[st.top()].ed && tmp.ed2+1 == kmer_list[st.top()].ed2)
				{
					tmp.ed++;
					tmp.ed2++;
				}
				else
				{
					lis_res.push_back(tmp);
					tmp.st=tmp.ed=tmp.st2=tmp.ed2=-1;
				}
			}
			st.pop();
		}
		if(tmp.st!=-1){
		lis_res.push_back(tmp);
		}
		chunk_align(lis_res,lis_res2);
		lis_res = lis_res2;
		for(i=0;i<lis_res2.size();i++){
			if(lis_res2[i].st2<K-3)
			{
				lis_res2[i].st-=lis_res2[i].st2;
				lis_res2[i].st2=0;
			}
			if(r.query_len-lis_res2[i].ed2<K-3)
			{
				lis_res2[i].ed+=r.query_len-1-lis_res2[i].ed2;
				lis_res2[i].ed2=r.query_len-1;
			}
//			printf("%d %d %d %d|| ",lis_res2[i].st,lis_res2[i].ed,lis_res2[i].st2,lis_res2[i].ed2);
		}

		tmp.st2=-1;tmp.ed2=-1;
		for(i=0;i<lis_res2.size();i++)
		{
			if(tmp.ed2>=lis_res2[i].st2)
			{
				int gap = tmp.ed2-lis_res2[i].st2+1;
				lis_res2[i].st+=gap,lis_res2[i].st2+=gap;
			}
			tmp=lis_res2[i];
		}
		chunk_align(lis_res2,lis_res);
		char t[5];
	        t[0] = 'A',t[1] = 'C',t[2] = 'G',t[3] = 'T';
		for(i=0;i<lis_res.size();i++)
		{
			int miss=0;
			/* with ins_fix the blocks may legitimately be apart in query
			   space; never let this walk-back chew into the inserted bases */
			int qfloor = (arg.ins_fix && i>0) ? lis_res[i-1].ed2+1 : 0;
			/* FIX: the ref cursor must stay in range too; st==0 read g_ref_seq[-1]. */
			while(lis_res[i].st2>qfloor && lis_res[i].st>0)
			{
				if(t[r.query[lis_res[i].st2-1]] == g_ref_seq[lis_res[i].st-1])
				{
					lis_res[i].st2--;
					lis_res[i].st--;
				}
				else
				{
					if(miss<=0){
					lis_res[i].st2--;
					lis_res[i].st--;
					lis_res[i].miss++;
					miss++;
					}
					else
					{
						break;
					}
				}
			}
			miss=0;
			/* FIX: symmetric guard on the right; ed+1 could run past the reference. */
			int qceil = (arg.ins_fix && i+1<(int)lis_res.size()) ? lis_res[i+1].st2-1 : r.query_len-1;
			while(lis_res[i].ed2<qceil && lis_res[i].ed+1<arg.ref_len)
			{
				if(t[r.query[lis_res[i].ed2+1]] == g_ref_seq[lis_res[i].ed+1])
				{
					lis_res[i].ed2++;
					lis_res[i].ed++;
				}
				else
				{
					if(miss<=0){
					lis_res[i].ed2++;
					lis_res[i].ed++;
					miss++;

					lis_res[i].miss++;
					}
					else
					{
						break;
					}

				}
			}
		}
		tmp.ed2=-1;
		for(i=0;i<lis_res.size();i++)
		{
			if(tmp.ed2>=lis_res[i].st2)
			{
				int gap = tmp.ed2-lis_res[i].st2+1;
				lis_res[i].st+=gap,lis_res[i].st2+=gap;
			}
			tmp=lis_res[i];
		}
		chunk_align(lis_res,lis_res2);
		lis_res = lis_res2;
		for(i=1;i<lis_res.size();i++){
			 int gap = lis_res[i].st2-lis_res[i-1].ed2-1;
			/* A gap in QUERY coordinates means the read carries bases the
			   reference does not: an insertion.  The pre-2026 code slides the
			   next block back over them, which turns every inserted base into
			   a mismatch (and then the miss>=8 guard throws the alignment
			   away), or clears the chain outright once the gap passes K+K/2.
			   A deletion has no query gap at all, which is why deletions never
			   hit either guard.  With ins_fix the gap is left in place and
			   make_cigar() writes it as an I. */
			if(arg.ins_fix && gap>0)
			{
				int rgap = lis_res[i].st-lis_res[i-1].ed-1;
				if(rgap<gap)
				{
					if(arg.ins_max>0 && gap-(rgap>0?rgap:0) > arg.ins_max)
					{
						lis_res.clear();
						return;
					}
					continue;           /* leave both cursors alone */
				}
			}
                        if(gap>K+K/2)
                        {
                                lis_res.clear();
                                return;
                        }
                        lis_res[i].st2-=gap;
                        lis_res[i].st-=gap;
		}
		chunk_align(lis_res,lis_res2);
		lis_res = lis_res2;
		int pos=0,miss=0,miss_tg=0,match_tg=0;
		for(i=0;i<lis_res.size();i++)
		{
			pos=0;
			for(int j=lis_res[i].st;j<=lis_res[i].ed;j++)
			{
				if(t[r.query[lis_res[i].st2+pos]]!=g_ref_seq[j])
				{
					miss++;
					if(arg.targetst+arg.cle-5<=j && j<=arg.targetst+arg.cle+5)
					{
						miss_tg++;
					}
				}
				else
				{
					if(arg.targetst+arg.cle-15<=j && j<=arg.targetst+arg.cle+15)
					{
						match_tg++;
					}
				}

				pos++;
			}
		}
		r.tg_miss= miss_tg;
		if(miss_tg>0){r.miss=1;}
		r.tg_match= match_tg;
		r.tags[1]= "XM:i:"+to_string(miss)+"\tXT:i:"+to_string(miss_tg);
		//if(miss-miss_tg>=8){lis_res.clear();}//last_modify
		if(miss>=8){lis_res.clear();}//last_modify
		for(i=1;i<lis_res.size();i++)
		{
				while(1)
				{
					/* FIX: this left-shift had no lower bound.  In a homopolymer
					   or short tandem repeat the condition can keep holding and
					   all four cursors walk past 0 into out-of-bounds reads.
					   Stop at the array boundary instead. */
					if(lis_res[i-1].ed2<=0 || lis_res[i-1].ed<=0 || lis_res[i].st<=1 || lis_res[i].st2<=0){break;}
					if((t[r.query[lis_res[i-1].ed2]] != g_ref_seq[lis_res[i-1].ed]) && (t[r.query[lis_res[i-1].ed2]] == g_ref_seq[lis_res[i].st-1]))
					{
						lis_res[i-1].ed2--;
						lis_res[i-1].ed--;
						lis_res[i].st--;
						lis_res[i].st2--;
					}
					else
					{
						break;
					}
				}
		}
		side_align(lis_res,r);
		prune_short_ends(lis_res);
		make_cigar(lis_res,K,r.cigar_string,r.query_len);
		read_classify(r,lis_res);

	}
	return;
}


