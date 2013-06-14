#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"


#define NBIN 120
#define Factor_relax 3


int main()
{
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	
//~ CATALOGUE Cat;
SUBCATALOGUE SubCat;
int Nsnap=99,i,j,k,pid,subid,hostsubid;

int Ngas,Ndm,Nother;
float *cen,rvir,rcen,r_this[NBIN],r_other[NBIN],r_back[NBIN],redges[NBIN+1];
int n_this[NBIN],n_other[NBIN],n_back[NBIN];
//~ int *PartHost;


char buf[1024];
FILE *fpsubprof;
logfile=stdout;
sprintf(outputdir,"%s/anal/tidal_tree",SUBCAT_DIR);	

	//~ load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	load_particle_data(Nsnap,snapdir);
	//~ fresh_ID2Index(&Cat,-1); 	
	fresh_ID2Index(&SubCat,-2);
	Ngas=header.npart[0];	Ndm=header.npart[1];	Nother=header.npart[2]+header.npart[3]+header.npart[4]+header.npart[5];
	//~ PartHost=mymalloc(sizeof(int)*Ndm);
	//~ for(i=0;i<Ndm;i++)
		//~ PartHost[i]=-1;
//~ for(subid=0;subid<SubCat.Nsubs;subid++)
//~ {
	//~ for(i=0;i<SubCat.SubLen[subid];i++)
		//~ PartHost[SubCat.PSubArr[subid][i]]=subid;
//~ }
	sprintf(buf,"%s/subprof_dm_%03d.%01dvir",outputdir,Nsnap,Factor_relax);
	myfopen(fpsubprof,buf,"w");


for(subid=0;subid<SubCat.Nsubs;subid++)
{
	for(i=0;i<NBIN;i++)	
	{
		n_this[i]=0;
		//~ n_other[i]=0;
		//~ n_back[i]=0;
		//~ r_this[i]=0;
		//~ r_other[i]=0;
		//~ r_back[i]=0;
	}
	if(SubCat.SubLen[subid])
	{
	cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
	rvir=pow(G*SubCat.SubLen[subid]*header.mass[1]/100/header.Hz/header.Hz,1.0/3)/header.time;
	rvir*=Factor_relax;
	for(i=0;i<=NBIN;i++)
	redges[i]=rvir*i/NBIN;	
	
	tree_tree_allocate(TREE_ALLOC_FACTOR*SubCat.SubLen[subid],SubCat.SubLen[subid]);
	maketree(SubCat.SubLen[subid],SubCat.PSubArr[subid]);
	tree_count_bin(cen,redges,NBIN,n_this,SubCat.PSubArr[subid]);
	tree_tree_free();
	}
	
	for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%d\t",n_this[i]);//number of particles that belong to this sub
	fprintf(fpsubprof,"\n");

}
		
fclose(fpsubprof);
return 0;
}
