//to search for the halo at a given location
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
	SUBCATALOGUE SubCat;
	SRCCATALOGUE SrcCat;
	int Nsubs,*pro2dest,Npro,Nsplitter,*sp2pro;
	float x[3]={149564.703, 173700.344, 148429.359},r=250.0; //CoM and rvir for halo 6702S51G78
	float *com,d;
	int Nsnap=0,i;
	int grpid,subid,pid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	for(i=0;i<Cat.Ngroups;i++)
	{
		subid=SubCat.GrpOffset_Sub[i];
		com=SubCat.Property[subid].CoM;
		d=distance(x,com);
		if(d<r)
		printf("ID: %d, Len: %d, Nsubs: %d, d/r: %f\n",i,Cat.Len[i],SubCat.GrpLen_Sub[i],d/r);
	}
	
	free_catalogue(&Cat);
	erase_sub_catalogue(&SubCat);
	return 0;
}
