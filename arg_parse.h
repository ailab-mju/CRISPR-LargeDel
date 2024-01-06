#include<string>
#ifndef AG
#define AG
struct Args
{
	char *refname,*ctrfastq_l,*ctrfastq_r,*mutfastq_l,*mutfastq_r,*tmpfile,*repeatfile,*view;
	int targetst,targeted,K,dist_cut,viewer,cle,ref_len;
	std::string dir,resctrcom,resmutcom,resmutonly,resctronly,resctr,resmut,ctrsite,mutsite,ctrsite_sm,mutsite_sm,outdir,genename,samplename;
};
#endif
void Arg_Parse(Args *arg,char *argv[]);
