#include<string>
#ifndef AG
#define AG

/* Arg_Parse() return codes. */
#define ARG_OK     0   /* everything parsed and validated -> run           */
#define ARG_USAGE  1   /* usage was requested / too few args -> exit 0     */
#define ARG_ERROR  2   /* an argument is invalid -> message printed, exit 1*/

struct Args
{
	char *refname,*ctrfastq_l,*ctrfastq_r,*mutfastq_l,*mutfastq_r,*tmpfile,*repeatfile,*view;
	int targetst,targeted,K,dist_cut,viewer,cle,ref_len;
	/* large deletion calling knobs - see print_usage() for what they mean */
	int min_gap;      /* pair inner-gap floor for a geometry call   (400) */
	int anchor_min;   /* split-read anchor score floor, 0 = off      (13) */
	int cut_between;  /* require the cut to fall between the mates    (1) */
	int ins_fix;      /* record a query-side gap as an I, not overwrite (1) */
	int ins_max;      /* cap on that I in bp, 0 = no cap               (0) */
	int report;       /* write <gene>_<sample>_report.html + .tsv      (1) */
	int large_del_len;/* split-read D length that counts as large (50) -> MIN_skip */
	int common_mode;  /* 0 = local breakpoint test (default), 1 = k-means (paper / pre 2026-10-23) */
	int clip_rescue;  /* re-assign terminal clips to small indels seen in full reads (1) */
	/* detection build: switches from comparison_results/cores/cli_detect, on by default */
	int ld_confirm;   /* record a read-pair call only after the re-alignment
	                     that already follows it agrees                    (1) */
	int split_wt_margin; /* split read: the short side must fit the deletion
	                     by more than N mismatches over reading on; -1 = off (2) */
	std::string dir,resctrcom,resmutcom,resmutonly,resctronly,resctr,resmut,ctrsite,mutsite,ctrsite_sm,mutsite_sm,ctrsite_clip,mutsite_clip,outdir,genename,samplename,mutsummary;
};
#endif
extern std::string g_cmdline;   /* for the @PG header line */
void print_usage(const char *prog);
int Arg_Parse(Args *arg,char *argv[]);
