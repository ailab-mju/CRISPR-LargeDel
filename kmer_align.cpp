#include<stdio.h>
#include<string.h>
#include "align_wrap.h"
#include "report.h"


int main(int argc, char *argv[])
{
	/* `kmer_align report ...` combines the _report.tsv of earlier runs */
	if(argc>1 && strcmp(argv[1],"report")==0){ return report_main(argc-1,argv+1); }
	/* Arg_Parse() prints the usage block itself when the arguments are
	   missing or malformed, so nothing is echoed before it has had a look. */
	align_wrap t = align_wrap(argv);
	return t.status;
}
