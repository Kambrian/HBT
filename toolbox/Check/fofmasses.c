#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc, char** argv)
{
	CATALOGUE Cat;
		
	HBTInt Nsnap=0;
	HBTInt i;
	FILE *fp;
	char buf[1024];

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	
	sprintf(buf,"%s/anal/fofmass.%d", SUBCAT_DIR, (int)Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<Cat.Ngroups;i++)
	fprintf(fp,HBTIFMT"\n",Cat.Len[i]);
	fclose(fp);
	
	return 0;
}
