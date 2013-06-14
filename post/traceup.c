#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
//~ #include "subfindread.c"

int main(int nargc,char **argv)
{
//~ char subcatSdir[512]="/SANdisk1/wenting/mergetree/SIM6702/subcatS";
CATALOGUE Cat;
SUBCATALOGUE SubCat;
int i,Nsnap,pid,pid0,subid,fofid,Nsplitter,Npro;
int *sp2pro;
float dr,rvir0;
FILE *fp;
char buf[1024];
	
	logfile=stdout;
	if(nargc!=3)
	{printf("usage: traceup <Nsnap> <subid>\n");
	return 1;}
	Nsnap=atoi(argv[1]);
	subid=atoi(argv[2]);
	sprintf(buf, "%s/anal/trace/traceup/trace_%03d_%d.dat",SUBCAT_DIR,Nsnap,subid);
	if(!(fp= fopen(buf, "w")))
	{
		printf("can't open file `%s'\n", buf);
		exit(1);
	}
	
	Nsplitter=0;
	do{
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	if(subid>=SubCat.Nsubs)
	{
	load_sp2pro(Nsnap+1,&Npro,&Nsplitter,&sp2pro,SUBCAT_DIR);
	if(Npro!=SubCat.Nsubs)
	{
		printf("error loading sp2pro, Npro mismatch %d=%d\n",Npro,SubCat.Nsubs);
		exit(1);
	}
	subid=sp2pro[subid];
	free_sp2pro(sp2pro,Npro,Nsplitter);
	}
	fofid=SubCat.HaloChains[subid].HostID;
	if(SubCat.SubLen[subid])
	{
	rvir0=comoving_virial_radius(Cat.Len[fofid]);
	dr=distance(SubCat.Property[subid].CoM,SubCat.Property[SubCat.GrpOffset_Sub[fofid]].CoM);
	fprintf(fp,"%d:G(%d,%d),S(%d,%d),R%d,D(%f,%f)\n",Nsnap,fofid,Cat.Len[fofid],subid,SubCat.SubLen[subid],SubCat.SubRank[subid],dr,dr/rvir0);
	}
	else
	fprintf(fp,"%d:G(%d,%d),S(%d,%d),R%d,D(%f,%f)\n",Nsnap,fofid,Cat.Len[fofid],subid,SubCat.SubLen[subid],SubCat.SubRank[subid],0.,0.);
	fflush(fp);
	subid=SubCat.HaloChains[subid].ProSubID;
	//~ Nsplitter=SubCat.Nsplitter;//this is the number of splitters from last snap
	free_catalogue(&Cat);
	for(i=0;i<SubCat.Nsubs;i++)
	if(SubCat.SubLen[i])
		free(SubCat.PSubArr[i]);
	free_sub_catalogue(&SubCat);
	Nsnap--;
	}while(subid>=0&&Nsnap>=0);
	fclose(fp);
	return 0;
}	
