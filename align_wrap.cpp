#include<stdio.h>
#include<thread>
#include<iostream>
#include<list>
#include<stdlib.h>
#include<unordered_map>
#include<algorithm>
#include<string>
#include<string.h>
#include<vector>
#include<time.h>
#include<stdlib.h>
#include "align.h"
#include "func.h"
#include "report.h"
#include "post.h"
#include "align_wrap.h"
#include "samout.h"
int skip_read,is_skip,skip;//mm=-1
string reff,dir;
char view_q[1000];
/* FIX: post.cpp declared these as `extern int`.  The mismatched type was
   benign on little-endian targets but still UB; both sides are int now. */
int mut_only=0,ctr_only=0,mut_common=0,ctr_common=0;
Args arg;
SamOut *f_out,*f_out2,*f_out3;
char g_ref_seq[100000],bf[1000];
/* FIX: the original ignored each fgets() result and unconditionally wrote
   buf[strlen(buf)-1].  Two ways that goes wrong:
     - at EOF fgets() leaves the buffer untouched, and because char_to_num()
       encodes 'A' as the byte 0x00, strlen() of a previously parsed read can
       be 0 -> buf[-1] is written, one byte before the array.
     - a truncated/blank final line gives the same underflow.
   Now every read is checked and the trailing newline is stripped only when
   it is actually there.  Well-formed FASTQ parses identically. */
inline int read_fastq(FILE *in,read &rr)
{
	if(fgets(rr.query_name,250,in)==NULL){return 0;}
	int l = strlen(rr.query_name);
	if(l>0 && rr.query_name[l-1]=='\n'){rr.query_name[--l]='\0';}
	for(int i=0;i<l;i++)
	{
		if(rr.query_name[i]==' ')
		{
			rr.query_name[i]='\0';
			break;
		}
	}
	if(fgets(rr.query,250,in)==NULL){return 0;}
	l = strlen(rr.query);
	if(l>0 && rr.query[l-1]=='\n'){rr.query[--l]='\0';}
	if(l==0){return 0;}
	if(fgets(bf,250,in)==NULL){return 0;}
	rr.strand=bf[0];
	if(fgets(rr.query_qual,250,in)==NULL){return 0;}
	int lq = strlen(rr.query_qual);
	if(lq>0 && rr.query_qual[lq-1]=='\n'){rr.query_qual[lq-1]='\0';}
	return 1;
}
/* per-sample tallies of --ld_confirm, printed at the end of process() */
long long g_confirm_tested=0,g_confirm_gap=0,g_confirm_split=0,g_confirm_vetoed=0,g_confirm_fallback=0,g_confirm_dropped=0;
int align_pair(int K, unordered_map<long long int,vector<int>> &kmer_ref, unordered_map<long long int,vector<int>> &kmer_ref2, read &r1, read &r2,vector<read> *split,int rev,result &res_cont)
{
	r1.map_f=-1,r2.map_f=-1;
	r1.pat=r2.pat="UNCALLED";   /* never let a tag leak from the previous pair; */
	                            /* overwritten as soon as the pair is classified */
	align_query(K,kmer_ref,r1,kmer_ref2,0);
	if(r1.map_f==-1){return 0;}
	align_query(K,kmer_ref,r2,kmer_ref2,0);
	if(r2.map_f==-1){return 0;}
	int confirmed=0; /* --ld_confirm re-aligned this pair already */
	if(r1.map_f!=-1 && r2.map_f!=-1)
	{
		if(r1.map_f==1 && r2.map_f==1)
		{
			make_sam(r1,r2,rev,arg.genename);
			/* --ld_confirm: a pair that would be called here is re-aligned
			   without pre_chunk first.  cli/ runs that same re-alignment only
			   after the call is recorded, so a mate pre_chunk put in the wrong
			   repeat copy still made the call.  A mate the re-alignment loses
			   keeps its first-pass placement only if that passes the aligner's
			   own limit (< 8 mismatches, align_query); otherwise the pair has
			   nothing to call from. */
			if(arg.ld_confirm && mate_gap_call(r1,r2))
			{
				read o1=r1,o2=r2;
				read *rr[2]={&r1,&r2}; read *oo[2]={&o1,&o2};
				g_confirm_tested++;
				align_query(K,kmer_ref,r1,kmer_ref2,1);
				align_query(K,kmer_ref,r2,kmer_ref2,1);
				int lost=0,kept_first=0;
				for(int j=0;j<2;j++)
				{
					if(rr[j]->map_f!=-1){continue;}
					int mm=placed_mm(*oo[j]);
					if(mm<8){*rr[j]=*oo[j];kept_first++;}else{lost=1;}
				}
				if(lost){g_confirm_dropped++;return 0;}
				confirmed=1;
				if(kept_first){g_confirm_fallback++;}
				if(r1.map_f==2||r2.map_f==2){g_confirm_split++;}
				else
				{
					make_sam(r1,r2,rev,arg.genename);
					if(mate_gap_call(r1,r2)){g_confirm_gap++;}else{g_confirm_vetoed++;}
				}
			}
		}
		if(r1.map_f==1 && r2.map_f==1)
		{
			int st1=r1.ref_st,ed1=r1.ref_ed,st2=r2.ref_st,ed2=r2.ref_ed;
			int res = long_inner_check(r1,r2,split,res_cont);
			if(res==-1){return 0;}
			if(res>=arg.min_gap)
			{
				if(confirmed)
				{
					/* already the final alignments: only undo the pair-level
					   coordinates long_inner_check() wrote into them */
					r1.map_f=r2.map_f=1;
					r1.ref_st=st1;r1.ref_ed=ed1;r2.ref_st=st2;r2.ref_ed=ed2;
				}
				else
				{
					align_query(K,kmer_ref,r1,kmer_ref2,1);
					align_query(K,kmer_ref,r2,kmer_ref2,1);
				}
			}
		}
		if(r1.map_f==2 || r2.map_f==2){
			if(!confirmed)
			{
				align_query(K,kmer_ref,r1,kmer_ref2,1);
				align_query(K,kmer_ref,r2,kmer_ref2,1);
			}
			if(r1.map_f!=-1&&r2.map_f!=-1){
			make_sam(r1,r2,rev,arg.genename);
			r1.pat=r2.pat="LDEL";
			split[0].push_back(r1);
			split[1].push_back(r2);
			}
		}
	}
	return 1;
}
void process(char *l_fastq, char *r_fastq, vector<read> *split, unordered_map<long long int,vector<int>> &kmer_ref,int K,unordered_map<long long int,vector<int>> &kmer_ref2,string output,result &res_cont)
{
	FILE *f2,*f3;
	SamOut *f4;
	read r1,r2;
	double inp_time=0,align_t=0,t1=0,t2=0;
	int total_align=0,total_read=0;
	long long unaligned_skipped=0;
	for(int i=0;i<5;i++){r1.tags.push_back("");r2.tags.push_back("");}
	int p2=0,p3=0;
	f2 = ka_open((dir+string(l_fastq)).c_str(),&p2);
	f3 = ka_open((dir+string(r_fastq)).c_str(),&p3);
	if(f2==NULL||f3==NULL)
	{
		fprintf(stderr,"[error] cannot open '%s' or '%s'\n",
			(dir+string(l_fastq)).c_str(),(dir+string(r_fastq)).c_str());
		exit(1);
	}
	f4 = samout_open(output);
//	mm_tot=mm_cnt=0;
	t2=clock();

	while(read_fastq(f2,r1)&&read_fastq(f3,r2))
	{
		inp_time+=(clock()-t2)/CLOCKS_PER_SEC;
		int pair_res=0;
		total_read++;
		t1=clock();
		parse_read(r1);
		parse_read(r2);
		reverse_comp(r1.query,r1.query_len);
		is_skip=0;
		pair_res=align_pair(K,kmer_ref,kmer_ref2,r1,r2,split,1,res_cont);    /// r1: rev   r2: normal
		if(pair_res==0){
			reverse_comp(r2.query,r2.query_len);
			reverse_comp(r1.query,r1.query_len);
			pair_res=align_pair(K,kmer_ref,kmer_ref2,r1,r2,split,0,res_cont); // r1 : normal r2: rev
		}
		/* align_pair() can return 1 after a re-alignment has invalidated one of
		   the mates; those records used to be written out with the "-1*"
		   sentinel as their CIGAR and stale POS/TLEN from the previous
		   attempt.  Keep them out of the file rather than emitting a record no
		   parser accepts. */
		if(pair_res==1 && !(read_aligned(r1) && read_aligned(r2)))
		{
			unaligned_skipped++;
			if(clip_rescue_waiting_full()){clip_rescue_attach_full(NULL,NULL);}
		}
		else if(clip_rescue_waiting_full())
		{
			/* a held-back terminal-clip pair: its whole-alignment records are
			   written by clip_rescue_finish(), with the final verdict */
			if(f4!=NULL&&pair_res==1){clip_rescue_attach_full(&r1,&r2);}
			else{clip_rescue_attach_full(NULL,NULL);}
		}
		else if(f4!=NULL&&pair_res==1)
		{
			samout_write(f4,r1);
			samout_write(f4,r2);
		}
		total_align+=pair_res;
		skip+=is_skip;
		t2=clock();
	}
	printf("total_align: %d\n total_read %d\n alignment_rate %.2lf%%\n",total_align,total_read,(double)((double)total_align*100/(double)total_read));
	if(unaligned_skipped)
	{
		printf(" incomplete pairs not written to the whole-alignment output: %lld\n",unaligned_skipped);
	}
	res_cont.t_align = total_align;
	if(arg.ld_confirm)
	{
		printf(" [detect] ld_confirm: %lld read-pair calls re-aligned first: %lld still a mate gap, %lld became a split read, %lld no longer a call, %lld dropped (a mate has no placement with < 8 mismatches); %lld used a first-pass mate\n",
			g_confirm_tested,g_confirm_gap,g_confirm_split,g_confirm_vetoed,g_confirm_dropped,g_confirm_fallback);
	}
	g_confirm_tested=g_confirm_gap=g_confirm_split=g_confirm_vetoed=g_confirm_fallback=g_confirm_dropped=0;
	clip_rescue_finish(res_cont,f4);
	//printf("align_t : %lf\ninp_time : %lf\n",align_t/CLOCKS_PER_SEC,inp_time);
	ka_close(f2,p2);
	ka_close(f3,p3);
	samout_close(f4);   /* NULL when the whole-alignment output is disabled */
}
/* ------------------------------------------------------------------ summary
   Per pattern read pair counts for one sample.

   tg_align_0 counts the pairs that span the cut site +/-13 bp, and
   tg_ins/tg_del/tg_miss/tg_clip are the mutated subset of exactly those
   pairs, so WT is simply the remainder.  tg_clip is kept out of the legacy
   `small_indel rate' line.  tg_align_0 now requires one of the mates to
   actually cover the cut window, so every row shares one denominator; pairs
   that merely straddle the cut are reported separately as `not covered'.  Large deletions are counted after the
   control/treated cluster comparison, i.e. only the sample specific ones.
   The denominator is the same one the reported rates already use
   (tg_align_0 + sample specific large deletions), so the five rows add up
   to 100%.  Clusters shared with the other sample are reported separately
   because they are excluded from the rate as PCR artefacts.             */
struct pattern_counts
{
	int wt,ins,del,sub,clip,ldel,lcommon,nocov,total,aligned;
};
static pattern_counts count_patterns(const result &r,int only_cnt,int common_cnt)
{
	pattern_counts p;
	p.ins = r.tg_ins;
	p.del = r.tg_del;
	p.sub = r.tg_miss;
	p.clip = r.tg_clip;
	p.nocov = r.tg_nocov;
	p.ldel = only_cnt;
	p.lcommon = common_cnt;
	p.wt = r.tg_align_0-(r.tg_ins+r.tg_del+r.tg_miss+r.tg_clip);
	if(p.wt<0){p.wt=0;}
	p.total = r.tg_align_0+only_cnt;
	p.aligned = r.t_align;
	return p;
}
static double pat_pct(int n,int d)
{
	if(d<=0){return 0.0;}
	return (double)n*100.0/(double)d;
}
static void print_pattern_row(const char *name,int n,int total)
{
	printf("  %-30s %12d %9.2lf%%\n",name,n,pat_pct(n,total));
}
void print_pattern_table(const char *label,const result &r,int only_cnt,int common_cnt)
{
	pattern_counts p = count_patterns(r,only_cnt,common_cnt);
	printf("\n###### %s sample mutation pattern summary\n",label);
	printf("  %-30s %12s %10s\n","Pattern","Read pairs","Percent");
	printf("  ---------------------------------------------------------\n");
	print_pattern_row("WT (no mutation)",p.wt,p.total);
	print_pattern_row("Small insertion",p.ins,p.total);
	print_pattern_row("Small deletion",p.del,p.total);
	print_pattern_row("Terminal insertion (clipped)",p.clip,p.total);
	print_pattern_row("Substitution",p.sub,p.total);
	print_pattern_row("Large deletion",p.ldel,p.total);
	printf("  ---------------------------------------------------------\n");
	print_pattern_row("Total analyzed",p.total,p.total);
	printf("  %-30s %12d %10s\n","Not covered (excluded)",p.nocov,"-");
	printf("  %-30s %12d %10s\n","Large del (common, excluded)",p.lcommon,"-");
	printf("  %-30s %12d %10s\n","Aligned pairs (all positions)",p.aligned,"-");
}
static void write_pattern_rows(FILE *f,const char *label,const result &r,int only_cnt,int common_cnt)
{
	pattern_counts p = count_patterns(r,only_cnt,common_cnt);
	fprintf(f,"%s\twt\t%d\t%.2lf\n",label,p.wt,pat_pct(p.wt,p.total));
	fprintf(f,"%s\tsmall_insertion\t%d\t%.2lf\n",label,p.ins,pat_pct(p.ins,p.total));
	fprintf(f,"%s\tsmall_deletion\t%d\t%.2lf\n",label,p.del,pat_pct(p.del,p.total));
	fprintf(f,"%s\tterminal_insertion_clipped\t%d\t%.2lf\n",label,p.clip,pat_pct(p.clip,p.total));
	fprintf(f,"%s\tsubstitution\t%d\t%.2lf\n",label,p.sub,pat_pct(p.sub,p.total));
	fprintf(f,"%s\tlarge_deletion\t%d\t%.2lf\n",label,p.ldel,pat_pct(p.ldel,p.total));
	fprintf(f,"%s\ttotal_analyzed\t%d\t%.2lf\n",label,p.total,pat_pct(p.total,p.total));
	fprintf(f,"%s\tnot_covered_excluded\t%d\tNA\n",label,p.nocov);
	fprintf(f,"%s\tlarge_deletion_common_excluded\t%d\tNA\n",label,p.lcommon);
	fprintf(f,"%s\taligned_pairs_all_positions\t%d\tNA\n",label,p.aligned);
}
void write_pattern_tsv(const string &path,const result &ctr,int ctr_o,int ctr_c,const result &mut,int mut_o,int mut_c)
{
	FILE *f = fopen(path.c_str(),"w");
	if(f==NULL)
	{
		fprintf(stderr,"[warning] could not write '%s'\n",path.c_str());
		return;
	}
	fprintf(f,"sample\tpattern\tread_pairs\tpercent\n");
	write_pattern_rows(f,"control",ctr,ctr_o,ctr_c);
	write_pattern_rows(f,"treated",mut,mut_o,mut_c);
	fclose(f);
	printf("\nmutation pattern summary written to %s\n",path.c_str());
}
void print_result(string type,result r)
{
	//printf("%s large del size : %d\n",type.c_str(),r.long_d);
	//printf("%s normal : %d %lf\n",type.c_str(),r.normal,(double)((double)r.normal/(double)r.t_align));
	//if(type=="mut"){
	//printf("%s sample Large Deletion Rate (removed false-postive) %.2lf\n",type.c_str(),(double)((double)r.tg_align_0/(double)r.t_align));
	//}
//	printf("%s tg_align_1 : %d %lf\n",type.c_str(),r.tg_align_1,(double)((double)r.tg_align_1/(double)r.t_align));

	if(type=="ctr"){

	/* FIX: this passed type.c_str() as an extra vararg for a format with a
	   single %.2lf.  On x86-64 SysV the double still landed in xmm0 so it
	   happened to print correctly; wasm passes varargs in order, so printf
	   read the char* as the double and printed 0.00. */
	printf("Control sample Large Deletion Rate (wihtout filtering) : %.2lf%%\n",(double)((double)r.tg_align_0*100/(double)r.t_align));
	printf("Control sample small_indel rate: %.2lf%%\n",(double)((double)(r.tg_ins+r.tg_del+r.tg_miss))*100/((double)(r.tg_align_0+r.only)));}
	//printf("%s skip %d\n",type.c_str(),skip);
}
align_wrap::align_wrap(char *argv[])
{
	vector<read> split_ctr[2],split_mut[2];
	unordered_map<long long int,vector<int>> kmer_ref,kmer_ref2;
	result ctr_r={0,},mut_r={0,};
	double t1 = clock(),t2=clock();
	status = 0;
	int arg_rc = Arg_Parse(&arg,argv);
	if(arg_rc!=ARG_OK)
	{
		status = (arg_rc==ARG_USAGE)?0:1;
		return;
	}
	for(int i=1;argv[i]!=NULL;i++)
	{
		printf("---%s\n",argv[i]);
	}
	/// make ref
	make_ref(arg.refname,arg.K,kmer_ref);  //k-mer reference
	make_ref(arg.refname,17,kmer_ref2); //17-mer reference unique assume
	//printf("ref time %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	printf("#######k-mer hashtable finished\n");fflush(stdout);
	/// ref end
	t1=clock();
	skip=0;

	//ctr align
	if(arg.ctrsite!=""){f_out=samout_open(arg.ctrsite);f_out2=samout_open(arg.ctrsite_sm);f_out3=samout_open(arg.ctrsite_clip);}
	printf("#######Control sample alignment start\n");
	process(arg.ctrfastq_l,arg.ctrfastq_r,split_ctr,kmer_ref,arg.K,kmer_ref2,arg.resctr,ctr_r);
	printf("######Control sample alignment end\n");fflush(stdout);
	printf("Control sample alignment time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	print_result("ctr",ctr_r);
	t1 = clock();
	skip=0;

	//mut align
		/* the control's three files have to be finished before the same globals
	   are pointed at the treated sample's */
	samout_close(f_out);samout_close(f_out2);samout_close(f_out3);
	f_out=f_out2=f_out3=NULL;
	if(arg.mutsite!=""){f_out=samout_open(arg.mutsite);f_out2=samout_open(arg.mutsite_sm);f_out3=samout_open(arg.mutsite_clip);}
	printf("######Cas9 Mutation sample alignment start\n");
	process(arg.mutfastq_l,arg.mutfastq_r,split_mut,kmer_ref,arg.K,kmer_ref2,arg.resmut,mut_r);
	printf("######Cas9 Mutation sample alignment end\n");fflush(stdout);
	printf("Cas9 Mutation sample alignment time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	print_result("mut",mut_r);
	t1 = clock();
	samout_close(f_out);samout_close(f_out2);samout_close(f_out3);
	f_out=f_out2=f_out3=NULL;
	printf("######post process start\n");
	post_process(split_ctr,split_mut,arg,ctr_r.t_align,mut_r.t_align);
	printf("######post process end\n");fflush(stdout);


	printf("post time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	printf("whole time : %lf\n",(clock()-t2)/CLOCKS_PER_SEC);

	//printf("mut tg_align_2 : %d %lf--%.3lf\n",mut_r.tg_align_2,(double)((double)mut_r.tg_align_2/(double)mut_r.t_align),(double)((double)post_process)/(double)(mut_r.tg_align_2+mut_only));
	//printf("mut tg_align_0 : %d %lf--%.3lf\n",mut_r.tg_align_0,(double)((double)mut_r.tg_align_0/(double)mut_r.t_align),(double)((double)mut_only)/(double)(mut_r.tg_align_0+mut_only));
	//printf("mut tg_align_1 : %d %lf--%.3lf\n",mut_r.tg_align_1,(double)((double)mut_r.tg_align_1/(double)mut_r.t_align),(double)((double)mut_only)/(double)(mut_r.tg_align_1+mut_only));
	//printf("Cas9 tg_align_0 : %d %lf--%.3lf\n",mut_r.tg_align_0,(double)((double)mut_r.tg_align_0/(double)mut_r.t_align),(double)((double)mut_only)/(double)(mut_r.tg_align_0+mut_only));
	printf("Cas9 Mutation sample Large Deletion Rate (removed false-postive) : %.2lf%%\n",(double)((double)mut_only)*100/(double)(mut_r.tg_align_0+mut_only));
	printf("Cas9 Mutation sample small_indel rate: %.2lf%%\n",((double)(mut_r.tg_ins+mut_r.tg_del+mut_r.tg_miss))*100/((double)mut_r.tg_align_0+mut_only));

	/* Per pattern counts.  Printed after the existing rate lines so that any
	   downstream parser keyed on those lines keeps working unchanged. */
	print_pattern_table("Control",ctr_r,ctr_only,ctr_common);
	print_pattern_table("Cas9 Mutation",mut_r,mut_only,mut_common);
	write_pattern_tsv(arg.mutsummary,ctr_r,ctr_only,ctr_common,mut_r,mut_only,mut_common);
	{
		/* HTML report - same numbers as the table above, plus the size
		   distributions tallied as the SAM records were written */
		pattern_counts a = count_patterns(ctr_r,ctr_only,ctr_common);
		pattern_counts b = count_patterns(mut_r,mut_only,mut_common);
		ka_pat pc = {a.wt,a.ins,a.del,a.sub,a.clip,a.ldel,a.lcommon,a.nocov,a.total,a.aligned};
		ka_pat pm = {b.wt,b.ins,b.del,b.sub,b.clip,b.ldel,b.lcommon,b.nocov,b.total,b.aligned};
		report_finish(pc,pm);
	}
	fflush(stdout);

//	printf("mut small_indel : %d %d %.4lf %d %.4lf %d %.4lf %d %.4lf\n",mut_r.tg_align_0+mut_only,mut_r.tg_ins,(double)((double)mut_r.tg_ins)/(double)(mut_r.tg_align_0+mut_only), mut_r.tg_del, (double)((double)mut_r.tg_del)/(double)(mut_r.tg_align_0+mut_only), mut_r.tg_ins+mut_r.tg_del,(double)((double)mut_r.tg_ins+mut_r.tg_del)/(double)(mut_r.tg_align_0+mut_only),mut_r.tg_ins+mut_r.tg_del+mut_r.tg_miss,((double)(mut_r.tg_ins+mut_r.tg_del+mut_r.tg_miss))/((double)mut_r.tg_align_0+mut_only));
}
align_wrap::~align_wrap(){
}
