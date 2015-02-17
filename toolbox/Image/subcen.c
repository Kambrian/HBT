#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

//~ #define SUBFIND_DIR "/home/kambrain/data/8213/subcatS"

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

int main()
{
	//~ CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	int Nsnap=59,i,j,subid,pid;
	char buf[1024];
	FILE *fpsub;

	logfile=stdout;

	sprintf(buf,"%s/anal/subcen_%03d",SUBCAT_DIR,Nsnap);
	if((fpsub=fopen(buf,"w"))==NULL)
	{
		printf("error: file open failed for %s!\n",buf);
		exit(1);
	}
	
	//~ load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	//~ fresh_ID2Index(&Cat,-1);
	fresh_ID2Index(&SubCat,-2);
	
	
for(i=0;i<SubCat.GrpOffset_Sub[100];i++)
{
		pid=SubCat.PSubArr[i][0];
		fprintf(fpsub,"%g\t%g\t%g\n",Pdat.Pos[pid][0],Pdat.Pos[pid][1],Pdat.Pos[pid][2]);
}
	fclose(fpsub);
return 0;
}
