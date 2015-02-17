#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

//~ #define SUBFIND_DIR "/home/kambrain/data/6702DM/subcatS"

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

#define NSUB 2

int main()
{
	//~ CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	HBTInt Nsnap,i,j,subid,pid;
	HBTInt sublist[NSUB]={7,88};//,259,444,1278,4752,4832,5310};
	//~ int sublist[NSUB]={12,476,683,2409,2501,2823,2824}; //correspondence:  2409:1372; 2501:1436; 2823,2824: 1648;
	//~ int sublist[NSUB]={ 2571,2576,3265,3268};
	//~ int sublist[NSUB]={6311,6317};
	//~ int sublist[NSUB]={43863};
	//~ int sublist[NSUB]={61120,70041};
	char buf[1024];
	FILE *fpsub;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/subparticles",SUBCAT_DIR);
	mkdir(outputdir,0755);	
	
	logfile=stdout;


	Nsnap=MaxSnap-1;
	//~ load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	//~ fresh_ID2Index(&Cat,-1);
	fresh_ID2Index(&SubCat,-2);
	free_PIDHash();
	
for(i=0;i<NSUB;i++)
{
	subid=sublist[i];
	sprintf(buf,"%s/subparticle_%03d_%d.bin",outputdir,Nsnap,subid);
	if((fpsub=fopen(buf,"w"))==NULL)
	{
		printf("error: file open failed for %s!\n",buf);
		exit(1);
	}
	fwrite(&SubCat.SubLen[subid],sizeof(HBTInt),1,fpsub);
	for(j=0;j<SubCat.SubLen[subid];j++)
	{
		pid=SubCat.PSubArr[subid][j];
		fwrite(Pdat.Pos[pid],sizeof(HBTReal),3,fpsub);
	}
	fwrite(&SubCat.SubLen[subid],sizeof(HBTInt),1,fpsub);
	fclose(fpsub);
}
return 0;
}
