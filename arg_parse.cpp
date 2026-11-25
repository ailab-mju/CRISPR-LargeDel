#include<stdlib.h>
#include<stdio.h>
#include<string>
#include<string.h>
#include<errno.h>
#include<limits.h>
#include<dirent.h>
#include<vector>
#include "arg_parse.h"
#include "samout.h"
using namespace std;

string g_cmdline;

/* argv is NULL terminated both for main() and for the wasm entry point, so the
   count can be recovered without changing every caller's signature. */
static int argv_count(char *argv[])
{
	int n=0;
	if(argv==NULL){return 0;}
	while(argv[n]!=NULL){n++;}
	return n;
}

void print_usage(const char *prog)
{
	if(prog==NULL||prog[0]=='\0'){prog="kmer_align";}
	printf(
"CRISPR-LargeDel : large deletion / small indel caller for Cas9 amplicon data\n"
"\n"
"USAGE\n"
"  %s <ref.fa> <ctr_R1.fastq> <ctr_R2.fastq> <mut_R1.fastq> <mut_R2.fastq>\n"
"               <target_start> <target_end> <cleavage_offset>\n"
"               <gene_name> <sample_name> <K> <out_dir> [--option value ...]\n"
"\n"
"  Easier: use the wrapper script, which fills in paths and sorts the output.\n"
"    bash run_kmer_align.sh <gene> <sample> <target_start> <target_end>\n"
"                           <cleavage_offset> <ctr_prefix> <mut_prefix> <K>\n"
"\n"
"REQUIRED ARGUMENTS (positional, in this order)\n"
"   1  ref.fa            reference amplicon FASTA, one sequence, < 100000 bp\n"
"   2  ctr_R1.fastq      control sample, read 1\n"
"   3  ctr_R2.fastq      control sample, read 2\n"
"   4  mut_R1.fastq      Cas9 treated sample, read 1\n"
"   5  mut_R2.fastq      Cas9 treated sample, read 2\n"
"   6  target_start      protospacer start, 0-based offset into ref.fa\n"
"   7  target_end        protospacer end, must be > target_start\n"
"   8  cleavage_offset   cut site measured from target_start (e.g. 3 for\n"
"                        SpCas9 with a 20 bp protospacer, i.e. 17 bp in)\n"
"   9  gene_name         used as the SAM reference name and file prefix\n"
"  10  sample_name       second half of the output file prefix\n"
"  11  K                 k-mer length for seeding, 5..31 (11 works well for\n"
"                        150 bp reads; a second 17-mer index is always built)\n"
"  12  out_dir           existing, writable directory for the results\n"
"\n"
"OPTIONAL ARGUMENTS\n"
"  Give any of these by name, in any order and anywhere after the program\n"
"  name, as --name value or --name=value; unnamed ones keep their default.\n"
"  The older positional form still works: the same options in the order\n"
"  listed are arguments 13-26 (shown as [#13] ... [#26]). If an option is\n"
"  given both ways the two values must agree.\n"
"\n"
"  --dist_cut N          [#13]\n"
"                        max distance (bp) between the cut site and the\n"
"                        nearest read end of a large-deletion call (split\n"
"                        read or read pair); the bigger the deletion, the\n"
"                        farther its reads sit from the cut. Default 8000:\n"
"                        on the simulations 4000 and up keeps every deletion\n"
"                        up to 7 kb, 1000 lost most of those >= 3 kb. Small\n"
"                        indels do not depend on it. 1000 = the default\n"
"                        before 2026-10-23.\n"
"  --full_sam 0|1        [#14]\n"
"                        1 = also write the two whole-alignment SAM files\n"
"                        (default), 0 = skip them. They are ~110 MB per\n"
"                        150k read pairs and nothing downstream reads them.\n"
"  --out_fmt FMT         [#15]\n"
"                        bam (default) = coordinate-sorted BAM + .bai only,\n"
"                        sam = plain SAM, both = write each output twice.\n"
"                        Sorting and indexing happen inside kmer_align, so no\n"
"                        samtools call is needed.\n"
"  --min_gap N           [#16]\n"
"                        minimum inner gap (bp) between the two mates before\n"
"                        the pair is called a large deletion on geometry alone\n"
"                        (the paper's read-pair category 2). Default 400.\n"
"                        NOTE this is an INSERT SIZE floor, not the paper's\n"
"                        100 bp large-deletion definition: with a ~300 bp\n"
"                        library a gap below ~400 bp is an ordinary fragment.\n"
"                        Lowering it to 100 makes normal WT pairs look deleted\n"
"                        (measured: 4x more false calls, 10 pp less sensitivity).\n"
"  --anchor_min N        [#17]\n"
"                        minimum anchor quality for a SPLIT read large deletion,\n"
"                        scored as (aligned block length) - 4 x (mismatches) on\n"
"                        the worst block. Default 13 (= MIN_match, i.e. a 13 bp\n"
"                        anchor must be perfect, a 17 bp anchor may carry one\n"
"                        mismatch). 0 disables the test. Raise to 16 for\n"
"                        artefact heavy amplicons.\n"
"  --cut_between 0|1     [#18]\n"
"                        1 (default) = the cleavage site must fall strictly\n"
"                        between the two mates for a geometry call, as the\n"
"                        paper specifies; 0 = only require the gap to overlap\n"
"                        the protospacer (the pre-2026 behaviour).\n"
"  --ins_fix 0|1         [#19]\n"
"                        1 (default) = a gap in query coordinates between two\n"
"                        seed blocks is recorded as an I in the CIGAR. 0 =\n"
"                        the pre-2026 behaviour, which slid the next block\n"
"                        back over the inserted bases, turned every one of\n"
"                        them into a mismatch and then discarded the read.\n"
"                        With 0, insertions of 9 bp and up are largely lost\n"
"                        (measured: 46%% of 9-20 bp insertions seen, against\n"
"                        97%% with 1) and no insertion longer than 16 bp gets\n"
"                        a size. Deletions, wild type and the wrong-call rate\n"
"                        are identical either way. Use 0 only to reproduce\n"
"                        pre-2026 results.\n"
"  --ins_max N           [#20]\n"
"                        largest insertion (bp) to accept; a junction needing\n"
"                        a longer I is discarded instead of aligned. Default\n"
"                        0 = no cap. Only meaningful with --ins_fix 1.\n"
"  --report 0|1          [#21]\n"
"                        1 (default) = also write <gene>_<sample>_report.html,\n"
"                        a self-contained HTML report (indel size spectrum,\n"
"                        composition, large deletion sizes, and insertions of\n"
"                        8 bp and up checked against the control's raw reads),\n"
"                        plus _report.tsv for combining runs later with\n"
"                          %s report -o all.html [--compare DIR...] -- DIR...\n"
"                        0 = skip both (and the extra pass over the FASTQs).\n"
"  --large_del_len N     [#22]\n"
"                        minimum deletion length (bp) inside ONE read for that\n"
"                        read to be called a split-read large deletion; shorter\n"
"                        deletions stay small indels. Default 50, the value\n"
"                        the original tool used (a compile-time constant,\n"
"                        stored internally as MIN_skip): on the stress\n"
"                        simulation 50-99 bp deletions were detected 69%% at\n"
"                        50 and 40%% at 100. 100 = the paper's large-deletion\n"
"                        definition (the default during 2026-10). Independent of\n"
"                        --min_gap, which governs read-pair calls.\n"
"  --clip_rescue 0|1     [#23]\n"
"                        1 (default) = a read whose small deletion sits within\n"
"                        13 bp of its end (so it only shows a terminal soft\n"
"                        clip) is re-assigned to that deletion when the same\n"
"                        allele is seen with full anchors in >= 3 other pairs\n"
"                        of the sample and explains the clipped bases exactly;\n"
"                        insertions the same way. Without it ~10%% of 1-10 bp\n"
"                        deletions were reported as terminal insertions.\n"
"                        0 = the behaviour before 2026-10-23.\n"
"  --common_mode MODE    [#24]\n"
"                        how a large deletion is judged to be shared with the\n"
"                        control (an artefact) rather than CRISPR derived.\n"
"                        local (default) = per read pair: the control has > 3\n"
"                        pairs with both breakpoints within max(20 bp, 2%% of\n"
"                        the span) AND at least half the treated sample's rate\n"
"                        there. Deterministic, independent of read order.\n"
"                        kmeans = the paper's cluster-centroid comparison; its\n"
"                        only/common split for 2 kb+ deletions moves by up to\n"
"                        30 pp when 1%% of the reads change.\n"
"  --ld_confirm 0|1      [#25]\n"
"                        1 (default) = a read-pair (mate gap) call is recorded only\n"
"                        after the re-alignment without pre_chunk, which follows\n"
"                        such a pair anyway, still shows the gap; a pair that\n"
"                        no longer does is counted like any other. 0 = off.\n"
"  --split_wt_margin N   [#26]\n"
"                        a split read is a large deletion only when the bases on\n"
"                        its short side fit the far block with more than N fewer\n"
"                        mismatches than reading straight on (+-3 bp); otherwise\n"
"                        it is aligned straight on and counted like any other\n"
"                        read. Default 2; -1 = off.\n"
"\n"
"EXAMPLE (the simulated TRAC amplicon in test_data/; out_dir must exist)\n"
"  mkdir -p result/TRAC_SIM_11\n"
"  %s test_data/TRAC.fa \\\n"
"      test_data/TRAC_SIM_ctr_1.fastq.gz test_data/TRAC_SIM_ctr_2.fastq.gz \\\n"
"      test_data/TRAC_SIM_mut_1.fastq.gz test_data/TRAC_SIM_mut_2.fastq.gz \\\n"
"      2952 2971 17 TRAC SIM 11 result/TRAC_SIM_11\n"
"\n"
"  with SAM output and the paper's 100 bp split-read cutoff instead of the\n"
"  default 50 (into an existing folder out/):\n"
"  %s ref.fa c_R1.fastq.gz c_R2.fastq.gz m_R1.fastq.gz m_R2.fastq.gz \\\n"
"      2952 2971 17 TRAC S1 11 out --out_fmt sam --large_del_len 100\n"
"\n"
"OUTPUT (written into out_dir, prefixed <gene_name>_<sample_name>_, as\n"
" coordinate-sorted .bam + .bam.bai unless out_fmt says otherwise)\n"
"  _mut_only      large deletions seen only in the treated sample\n"
"                     -> these are the CRISPR derived calls\n"
"  _mut_common    large deletions also present in the control\n"
"                     -> treated as PCR artefacts and excluded\n"
"  _ctr_only      large deletions seen only in the control\n"
"  _ctr_common    control side of the common large deletions\n"
"  _site_mut      reads spanning the cut site with no mutation (WT)\n"
"  _site_mut_sm   reads spanning the cut site with a small indel or\n"
"                     substitution within +/-13 bp\n"
"  _site_mut_clip reads whose terminal soft clip (>=5 bp) starts within\n"
"                     +/-13 bp of the cut site -> insertion at the read edge\n"
"  _site_ctr      same two files for the control sample\n"
"  _site_ctr_sm\n"
"  _ctr / _mut        every aligned pair (omitted when full_sam is 0)\n"
"  _mutation_summary.tsv  per pattern read pair counts and percentages,\n"
"                         the same table that is printed at the end of the run\n"
"  _report.html   self-contained HTML report (unless report is 0)\n"
"  _report.tsv    the same numbers, for `%s report` to combine runs\n"
,prog,prog,prog,prog,prog);
}

static int parse_int(const char *s,const char *name,int pos,int *out)
{
	char *end=NULL;
	long v;
	if(s==NULL||s[0]=='\0')
	{
		fprintf(stderr,"[argument error] #%d %s is empty\n",pos,name);
		return 0;
	}
	errno=0;
	v = strtol(s,&end,10);
	if(end==s||*end!='\0')
	{
		fprintf(stderr,"[argument error] #%d %s: '%s' is not an integer\n",pos,name,s);
		return 0;
	}
	if(errno==ERANGE||v<INT_MIN||v>INT_MAX)
	{
		fprintf(stderr,"[argument error] #%d %s: '%s' is out of range\n",pos,name,s);
		return 0;
	}
	*out=(int)v;
	return 1;
}

static int check_readable(const char *path,const char *name,int pos)
{
	FILE *f;
	if(path==NULL||path[0]=='\0')
	{
		fprintf(stderr,"[argument error] #%d %s is empty\n",pos,name);
		return 0;
	}
	f=fopen(path,"r");
	if(f==NULL)
	{
		fprintf(stderr,"[argument error] #%d %s: cannot open '%s' for reading (%s)\n",pos,name,path,strerror(errno));
		return 0;
	}
	fclose(f);
	return 1;
}

/* The directory must already exist - the wrapper script mkdir -p's it.  Probing
   with the real summary file avoids leaving a stray temp file behind. */
static int check_outdir(const string &d,const string &probe)
{
	DIR *dp;
	FILE *f;
	if(d.empty())
	{
		fprintf(stderr,"[argument error] #12 out_dir is empty\n");
		return 0;
	}
	dp=opendir(d.c_str());
	if(dp==NULL)
	{
		fprintf(stderr,"[argument error] #12 out_dir: cannot open directory '%s' (%s)\n",d.c_str(),strerror(errno));
		fprintf(stderr,"                 create it first, e.g.  mkdir -p %s\n",d.c_str());
		return 0;
	}
	closedir(dp);
	f=fopen(probe.c_str(),"a");
	if(f==NULL)
	{
		fprintf(stderr,"[argument error] #12 out_dir: '%s' is not writable (%s)\n",d.c_str(),strerror(errno));
		return 0;
	}
	fclose(f);
	return 1;
}

static int check_name(const string &s,const char *name,int pos)
{
	if(s.empty())
	{
		fprintf(stderr,"[argument error] #%d %s is empty\n",pos,name);
		return 0;
	}
	if(s.find('/')!=string::npos)
	{
		fprintf(stderr,"[argument error] #%d %s: '%s' must not contain '/'\n",pos,name,s.c_str());
		return 0;
	}
	return 1;
}

int Arg_Parse(Args *arg,char *argv_in[])
{
	/* Named options (--name value / --name=value) for every OPTIONAL argument.
	   They are taken out of argv first; if the 12 required arguments are all
	   there, the positional list is then extended up to the highest named
	   position, filling any gap with that position's default.  The checks
	   further down therefore see exactly what the old positional form would
	   have given them, and both forms produce identical runs. */
	struct NamedOpt { const char *name; int pos; const char *def; };
	static const NamedOpt k_named[] = {
		{"dist_cut",13,"8000"},{"full_sam",14,"1"},{"out_fmt",15,"bam"},
		{"min_gap",16,"400"},{"anchor_min",17,"13"},{"cut_between",18,"1"},
		{"ins_fix",19,"1"},{"ins_max",20,"0"},{"report",21,"1"},
		{"large_del_len",22,"50"},{"clip_rescue",23,"1"},{"common_mode",24,"local"},
		{"ld_confirm",25,"1"},{"split_wt_margin",26,"2"}};
	const int n_named = (int)(sizeof(k_named)/sizeof(k_named[0]));
	static vector<char*> pos_argv;
	static vector<string> named_val;        /* index = position, "" = not given */
	int argc_in = argv_count(argv_in);
	int named_ok = 1, named_max = 0;
	pos_argv.clear(); named_val.assign(27,string());
	for(int i=0;i<argc_in;i++)
	{
		const char *a = argv_in[i];
		if(i==0 || strncmp(a,"--",2)!=0 || strcmp(a,"--help")==0){ pos_argv.push_back(argv_in[i]); continue; }
		const char *eq = strchr(a,'=');
		string key = eq ? string(a+2,eq-(a+2)) : string(a+2);
		int k=-1;
		for(int j=0;j<n_named;j++){ if(key==k_named[j].name){k=j;break;} }
		if(k<0)
		{
			string known;
			for(int j=0;j<n_named;j++){ known += (j?", --":"--"); known += k_named[j].name; }
			fprintf(stderr,"[argument error] unknown option '%s'\n                 known: %s\n",a,known.c_str());
			named_ok=0; continue;
		}
		const char *v = eq ? eq+1 : ((i+1<argc_in) ? argv_in[i+1] : NULL);
		if(!eq && v!=NULL){ i++; }
		if(v==NULL || v[0]=='\0' || (strncmp(v,"--",2)==0))
		{
			fprintf(stderr,"[argument error] --%s needs a value, e.g. --%s %s\n",k_named[k].name,k_named[k].name,k_named[k].def);
			named_ok=0; if(v!=NULL && !eq && strncmp(v,"--",2)==0){ i--; } continue;
		}
		int p = k_named[k].pos;
		if(!named_val[p].empty() && named_val[p]!=v)
		{
			fprintf(stderr,"[argument error] --%s given twice (%s, %s)\n",k_named[k].name,named_val[p].c_str(),v);
			named_ok=0; continue;
		}
		named_val[p]=v; if(p>named_max){named_max=p;}
	}
	if(named_ok && named_max>0 && (int)pos_argv.size()>=13)
	{
		for(int p=13;p<=26;p++)
		{
			if(named_val[p].empty()){ continue; }
			if(p<(int)pos_argv.size() && named_val[p]!=pos_argv[p])
			{
				const char *nm=""; for(int j=0;j<n_named;j++){ if(k_named[j].pos==p){nm=k_named[j].name;} }
				fprintf(stderr,"[argument error] %s given twice with different values (#%d = %s, --%s = %s)\n",
				        nm,p,pos_argv[p],nm,named_val[p].c_str());
				named_ok=0;
			}
		}
		for(int p=(int)pos_argv.size();p<=named_max;p++)
		{
			const char *def=""; for(int j=0;j<n_named;j++){ if(k_named[j].pos==p){def=k_named[j].def;} }
			if(named_val[p].empty()){ named_val[p]=def; }
			pos_argv.push_back((char*)named_val[p].c_str());
		}
	}
	pos_argv.push_back(NULL);
	char **argv = pos_argv.data();
	int argc = argv_count(argv);
	int ok=1;
	const char *prog = (argc>0)?argv[0]:"kmer_align";

	if(argc>1&&(strcmp(argv[1],"-h")==0||strcmp(argv[1],"--help")==0||strcmp(argv[1],"help")==0))
	{
		print_usage(prog);
		return ARG_USAGE;
	}
	if(argc<13)
	{
		if(argc>1)
		{
			fprintf(stderr,"[argument error] %d argument(s) given, 12 required.\n\n",argc-1);
		}
		print_usage(prog);
		return (argc<=1)?ARG_USAGE:ARG_ERROR;
	}

	g_cmdline.clear();
	for(int i=0;i<argc_in;i++){ if(i){g_cmdline+=" ";} g_cmdline+=argv_in[i]; }

	arg->refname = argv[1];
	arg->ctrfastq_l = argv[2];
	arg->ctrfastq_r = argv[3];
	arg->mutfastq_l = argv[4];
	arg->mutfastq_r = argv[5];
	arg->genename = string(argv[9]);
	arg->samplename = string(argv[10]);
	arg->outdir = string(argv[12]);

	ok &= check_readable(arg->refname,"ref.fa",1);
	ok &= check_readable(arg->ctrfastq_l,"ctr_R1.fastq",2);
	ok &= check_readable(arg->ctrfastq_r,"ctr_R2.fastq",3);
	ok &= check_readable(arg->mutfastq_l,"mut_R1.fastq",4);
	ok &= check_readable(arg->mutfastq_r,"mut_R2.fastq",5);
	ok &= parse_int(argv[6],"target_start",6,&arg->targetst);
	ok &= parse_int(argv[7],"target_end",7,&arg->targeted);
	ok &= parse_int(argv[8],"cleavage_offset",8,&arg->cle);
	ok &= check_name(arg->genename,"gene_name",9);
	ok &= check_name(arg->samplename,"sample_name",10);
	ok &= parse_int(argv[11],"K",11,&arg->K);

	if(ok)
	{
		if(arg->targetst<0)
		{
			fprintf(stderr,"[argument error] #6 target_start: must be >= 0, got %d\n",arg->targetst);ok=0;
		}
		if(arg->targeted<=arg->targetst)
		{
			fprintf(stderr,"[argument error] #7 target_end: must be > target_start (%d), got %d\n",arg->targetst,arg->targeted);ok=0;
		}
		if(arg->cle<0||arg->cle>arg->targeted-arg->targetst)
		{
			fprintf(stderr,"[argument error] #8 cleavage_offset: must be within 0..%d (target length), got %d\n",arg->targeted-arg->targetst,arg->cle);ok=0;
		}
		if(arg->K<5||arg->K>31)
		{
			fprintf(stderr,"[argument error] #11 K: must be 5..31, got %d\n",arg->K);ok=0;
		}
	}

	arg->resctronly = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_ctr_only";
	arg->resctrcom = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_ctr_common";
	arg->resmutonly = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_mut_only";
	arg->resmutcom = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_mut_common";
	arg->resctr = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_ctr";
	arg->resmut = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_mut";
	arg->ctrsite = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_ctr";
	arg->mutsite = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_mut";
	arg->ctrsite_sm = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_ctr_sm";
	arg->mutsite_sm = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_mut_sm";
	arg->ctrsite_clip = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_ctr_clip";
	arg->mutsite_clip = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_mut_clip";
	arg->mutsummary = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_mutation_summary.tsv";

	if(ok)
	{
		ok &= check_outdir(arg->outdir,arg->mutsummary);
	}

	if(argc<14)
	{
		arg->dist_cut = 8000;
	}
	else
	{
		ok &= parse_int(argv[13],"dist_cut",13,&arg->dist_cut);
		if(ok&&arg->dist_cut<0)
		{
			fprintf(stderr,"[argument error] #13 dist_cut: must be >= 0, got %d\n",arg->dist_cut);ok=0;
		}
	}
	/* Optional argv[14]: write the full per-read alignment SAMs (default 1).
	   The browser build passes 0 - those two files are by far the largest
	   output (~110 MB per 150k pairs) and nothing downstream reads them. */
	if(argc>=15)
	{
		int full_sam=1;
		ok &= parse_int(argv[14],"full_sam",14,&full_sam);
		if(ok&&full_sam!=0&&full_sam!=1)
		{
			fprintf(stderr,"[argument error] #14 full_sam: must be 0 or 1, got %d\n",full_sam);ok=0;
		}
		if(ok&&full_sam==0)
		{
			arg->resctr = "";
			arg->resmut = "";
		}
	}

	/* Defaults for the large-deletion knobs; overridden by argv[16..18]. */
	arg->min_gap     = 400;
	arg->anchor_min  = 13;
	arg->cut_between = 1;
	/* Insertion handling; overridden by argv[19..20]. */
	arg->ins_fix     = 1;
	arg->ins_max     = 0;
	arg->report      = 1;
	/* Split-read large-deletion length; overridden by argv[22]. */
	arg->large_del_len = 50;
	arg->clip_rescue   = 1;
	arg->common_mode   = 0;
	arg->ld_confirm = 1;
	arg->split_wt_margin = 2;

	if(argc>=16)
	{
		const char *f = argv[15];
		if(strcmp(f,"bam")==0){samout_set_format(OUT_BAM);}
		else if(strcmp(f,"sam")==0){samout_set_format(OUT_SAM);}
		else if(strcmp(f,"both")==0){samout_set_format(OUT_BOTH);}
		else
		{
			fprintf(stderr,"[argument error] #15 out_fmt: expected bam|sam|both, got '%s'\n",f);
			ok=0;
		}
	}
	if(argc>=17)
	{
		ok &= parse_int(argv[16],"min_gap",16,&arg->min_gap);
		if(ok&&arg->min_gap<1)
		{
			fprintf(stderr,"[argument error] #16 min_gap: must be >= 1, got %d\n",arg->min_gap);ok=0;
		}
	}
	if(argc>=18)
	{
		ok &= parse_int(argv[17],"anchor_min",17,&arg->anchor_min);
		if(ok&&arg->anchor_min<0)
		{
			fprintf(stderr,"[argument error] #17 anchor_min: must be >= 0, got %d\n",arg->anchor_min);ok=0;
		}
	}
	if(argc>=19)
	{
		ok &= parse_int(argv[18],"cut_between",18,&arg->cut_between);
		if(ok&&arg->cut_between!=0&&arg->cut_between!=1)
		{
			fprintf(stderr,"[argument error] #18 cut_between: must be 0 or 1, got %d\n",arg->cut_between);ok=0;
		}
	}
	if(argc>=20)
	{
		ok &= parse_int(argv[19],"ins_fix",19,&arg->ins_fix);
		if(ok&&arg->ins_fix!=0&&arg->ins_fix!=1)
		{
			fprintf(stderr,"[argument error] #19 ins_fix: must be 0 or 1, got %d\n",arg->ins_fix);ok=0;
		}
	}
	if(argc>=21)
	{
		ok &= parse_int(argv[20],"ins_max",20,&arg->ins_max);
		if(ok&&arg->ins_max<0)
		{
			fprintf(stderr,"[argument error] #20 ins_max: must be >= 0, got %d\n",arg->ins_max);ok=0;
		}
	}
	if(argc>=22)
	{
		ok &= parse_int(argv[21],"report",21,&arg->report);
		if(ok&&arg->report!=0&&arg->report!=1)
		{
			fprintf(stderr,"[argument error] #21 report: must be 0 or 1, got %d\n",arg->report);ok=0;
		}
	}
	if(argc>=23)
	{
		ok &= parse_int(argv[22],"large_del_len",22,&arg->large_del_len);
		if(ok&&arg->large_del_len<1)
		{
			fprintf(stderr,"[argument error] #22 large_del_len: must be >= 1, got %d\n",arg->large_del_len);ok=0;
		}
	}
	if(argc>=24)
	{
		ok &= parse_int(argv[23],"clip_rescue",23,&arg->clip_rescue);
		if(ok&&arg->clip_rescue!=0&&arg->clip_rescue!=1)
		{
			fprintf(stderr,"[argument error] #23 clip_rescue: must be 0 or 1, got %d\n",arg->clip_rescue);ok=0;
		}
	}
	if(argc>=25)
	{
		if(strcmp(argv[24],"local")==0){arg->common_mode=0;}
		else if(strcmp(argv[24],"kmeans")==0){arg->common_mode=1;}
		else{fprintf(stderr,"[argument error] #24 common_mode: expected local|kmeans, got '%s'\n",argv[24]);ok=0;}
	}
	if(argc>=26)
	{
		ok &= parse_int(argv[25],"ld_confirm",25,&arg->ld_confirm);
		if(ok&&arg->ld_confirm!=0&&arg->ld_confirm!=1)
		{
			fprintf(stderr,"[argument error] #25 ld_confirm: must be 0 or 1, got %d\n",arg->ld_confirm);ok=0;
		}
	}
	if(argc>=27)
	{
		ok &= parse_int(argv[26],"split_wt_margin",26,&arg->split_wt_margin);
		if(ok&&arg->split_wt_margin<-1)
		{
			fprintf(stderr,"[argument error] #26 split_wt_margin: must be >= -1, got %d\n",arg->split_wt_margin);ok=0;
		}
	}
	ok &= named_ok;

	if(ok && !samout_have_bam() && samout_format()!=OUT_SAM)
	{
		fprintf(stderr,"[note] this build has no htslib - writing SAM instead of BAM\n");
	}

	if(!ok)
	{
		fprintf(stderr,"\n");
		print_usage(prog);
		return ARG_ERROR;
	}
	return ARG_OK;
}
