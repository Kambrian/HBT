//this file produces raw data output
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
	FILE *fp;
	char buf[1024];
	int Nsnap=70;
	int i,j,pid;
	double partmass;
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/massfun",SUBCAT_DIR);	
	
	logfile=stdout;
	if(argc!=2)
	{
	printf("usage:%s [Snap]\n",argv[0]);
	exit(1);
	}
	Nsnap=atoi(argv[1]);

	sprintf(buf,"%s/sub_mass_%03d",outputdir,Nsnap);
	if(!(fp=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	load_particle_data(Nsnap,snapdir);
	fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCat,-2);
	partmass=header.mass[1];
	printf("dm %g,baryon %g, fraction %g,Nsplitter %d, Nquasi %d\n",
			partmass,header.mass[0],header.mass[0]/header.mass[1],SubCat.Nsplitter,SubCat.NQuasi);
	for(i=0;i<SubCat.Nsubs;i++)
	{
		fprintf(fp,"%g\t",SubCat.SubLen[i]*partmass);
		fprintf(fp,"%g\t%g\t%g\n",SubCat.Property[i].CoM[0],SubCat.Property[i].CoM[1],SubCat.Property[i].CoM[2]);
	}
	fclose(fp);
			
	sprintf(buf,"%s/group_offset_%03d",outputdir,Nsnap);
	if(!(fp=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
		for(i=0;i<SubCat.Ngroups;i++)
	{
		fprintf(fp,"%d\t%d\t%f\n",SubCat.GrpOffset_Sub[i],SubCat.GrpLen_Sub[i],Cat.Len[i]*partmass);
	}
	fclose(fp);
	
	return 0;
}
