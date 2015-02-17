#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

SUBCATALOGUE SubCatTmp;

int main()
{
//~ char outputdir[512]="/SANdisk5/kambrain/8213/fofpot";
char outputdir[512]="/home/kambrain/data/8213/fofpot";
	
int Nsnap;

CATALOGUE Cat;

int *fofPID,fofLen,recordsize,haloid,i;
float *pot,*fofpot, scale;
//~ double *pot,*fofpot, scale;
//~ int id_minpot; 
//~ float rvir;

char buf[1024]; FILE *fp; 
	
	logfile=stdout;	
	Nsnap=0;
	sprintf(buf,"%s/potfof.b20.%d.%04d",outputdir,RUN_NUM,snaplist[Nsnap]);
	if(!(fp=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}	
 	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	pot=mymalloc(sizeof(*pot)*Cat.Nids);
	scale=1./(1.+header.ztp);
	#pragma omp parallel for  private(haloid,fofPID,fofpot,fofLen,i)  schedule(dynamic) num_threads(6) 
	for(haloid=0;haloid<Cat.Ngroups;haloid++)
	{
		fofPID=Cat.PIDorIndex+Cat.Offset[haloid];
		fofpot=pot+Cat.Offset[haloid];
		fofLen=Cat.Len[haloid];
		tree_tree_allocate(TREE_ALLOC_FACTOR*fofLen,fofLen);
		maketree(fofLen,fofPID);
			for(i=0;i<fofLen;i++)
			{
				fofpot[i]=tree_treeevaluate_potential(fofPID[i],fofPID);
			    fofpot[i]=(PMass*fofpot[i]+PMass/SofteningHalo)*G/scale;	
			 }
		//~ id_minpot=fofPID[Dmin_of_vec(fofpot,fofLen)];
		 //~ printf("MinPot:%g,%g,%g\n",Pdat.Pos[id_minpot][0],Pdat.Pos[id_minpot][1],Pdat.Pos[id_minpot][2]);
		 //~ printf("fofcen:%g,%g,%g\n",Cat.HaloCen[0][haloid],Cat.HaloCen[1][haloid],Cat.HaloCen[2][haloid]);
		 //~ rvir=pow(G*fofLen*PMass/100/header.Hz/header.Hz,1.0/3)/scale;
		 //~ printf("rvir:%g,foflen:%d\n",rvir,fofLen);
	}
		
	recordsize=sizeof(float)*Cat.Nids;
	fwrite(&recordsize,sizeof(int),1,fp);
	fwrite(pot,sizeof(float),Cat.Nids,fp);
	fwrite(&recordsize,sizeof(int),1,fp);
	fclose(fp);
	return 0;
}
