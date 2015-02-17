#include <stdio.h>
#include <stdlib.h>
#include "subfindread.c"

int main()
{
	char inputdir[512]="/home/kambrain/data/6702/subcatS";
	char outputdir[1024]="/home/kambrain/data/6702/subcat/anal";
	int Nsnap=99,i;
	FILE *fp;
	char buf[1024];
	
	sprintf(buf,"%s/SubFind_mass_%03d",outputdir,Nsnap);
	if(!(fp=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
	
	subread(Nsnap,inputdir);
	for(i=0;i<Nsubhalos;i++)
	{
		fprintf(fp,"%d,%d,%g,%g,%g\n",SubLen[i],SubParentHalo[i],SubPos[i*3],SubPos[i*3+1],SubPos[i*3+2]);
	}
	fclose(fp);
	
	sprintf(buf,"%s/SubFindgroup_offset_%03d",outputdir,Nsnap);
	if(!(fp=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
		for(i=0;i<Ngroups;i++)
	{
		fprintf(fp,"%d\t%d\t%g\t%g\t%g\n",NsubPerHalo[i],FirstSubOfHalo[i],Halo_M_Crit200[i],Halo_R_Crit200[i],Halo_R_Mean200[i]);
	}
	fclose(fp);
	return 0;
	
	
}
