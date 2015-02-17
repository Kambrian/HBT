//extract fof data for raw massfunction estimates
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc,char **argv)
{
	CATALOGUE Cat;
SUBCATALOGUE SubCat;
	FILE *fpg,*fps;
	char buf[1024];
	int Nsnap=80;
	int i,j,pid;
	double partmass;
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/massfun",SUBCAT_DIR);	
	
	logfile=stdout;
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise set Nsnap inside\n",argv[0]);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	sprintf(buf,"%s/grptab_%03d",outputdir,Nsnap);
	if(!(fpg=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
	sprintf(buf,"%s/sublen_%03d",outputdir,Nsnap);
	if(!(fps=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	//~ load_particle_data(Nsnap,snapdir);
	//~ fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCatA,-2);
	//~ partmass=header.mass[1];
	//~ printf("dm %g,baryon %g, fraction %g,Nsplitter %d, Nquasi %d\n",partmass,header.mass[0],header.mass[0]/header.mass[1],SubCatA.Nsplitter,SubCatA.NQuasi);
	//~ partmass=0.000058;
	for(i=0;i<Cat.Ngroups;i++)
	{
		fwrite(Cat.Len+i,sizeof(int),1,fpg);
		fwrite(SubCat.GrpOffset_Sub+i,sizeof(int),1,fpg);
		fwrite(SubCat.GrpLen_Sub+i,sizeof(int),1,fpg);
	}
	for(i=0;i<SubCat.Nsubs;i++)
	{
		fwrite(SubCat.SubLen+i,sizeof(int),1,fps);
	}
	fclose(fpg);
	fclose(fps);
	
	return 0;
}
