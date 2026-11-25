#include<stdio.h>
#include "align_wrap.h"


int main(int argc, char *argv[])
{
	for(int i=1;i<argc;i++)
	{
		printf("---%s\n",argv[i]);
	}
	align_wrap t = align_wrap(argv);
}
