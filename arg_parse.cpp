#include<stdlib.h>
#include<string>
#include<string.h>
#include "arg_parse.h"
using namespace std;
void Arg_Parse(Args *arg,char *argv[])
{
	arg->refname = argv[1];
	arg->ctrfastq_l = argv[2];
	arg->ctrfastq_r = argv[3];
	arg->mutfastq_l = argv[4];
	arg->mutfastq_r = argv[5];
	arg->targetst = atoi(argv[6]);
	arg->targeted = atoi(argv[7]);
	arg->cle = atoi(argv[8]);
//	cle=arg->cle;
//	target_st=arg->targetst;
//	target_ed=arg->targeted;
	arg->genename = string(argv[9]);
	arg->samplename = string(argv[10]);
	arg->K= atoi(argv[11]);
	arg->outdir=string(argv[12]);
//	reff = string(arg->genename);
	arg->resctronly = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_ctr_only.sam";
	arg->resctrcom = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_ctr_common.sam";
	arg->resmutonly = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_mut_only.sam";
	arg->resmutcom = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_mut_common.sam";
	arg->resctr = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_ctr.sam";
	arg->resmut = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_mut.sam";
	arg->ctrsite = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_ctr.sam";
	arg->mutsite = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_mut.sam";
	arg->ctrsite_sm = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_ctr_sm.sam";
	arg->mutsite_sm = 	arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_site_mut_sm.sam";
	if(argv[13]==NULL)
	{
		arg->dist_cut = 1000;
	}
	else
	{
		arg->dist_cut=atoi(argv[13]);
	}
}
