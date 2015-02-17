//to get the min-pot and max-density center of halos
//output the id of minpot and maxdens particles
// ToDo: sph_dens() just runs too slow.... shit, got to improve this!!.......
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

CATALOGUE Cat;
int *IDrhomax,*IDpotmin;
void get_CenID(int grpid);

int main(int argc,char **argv)
{
	int Nsnap,i,j,pid,grpid;
	FILE *fp;
	char buf[1024];

	char outputdir[1024];
	
	sprintf(outputdir,"%s/profile/",SUBCAT_DIR);
	mkdir(outputdir,0755);
	logfile=stdout;
	
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
		sprintf(buf,"%s/profile/haloCenID_%03d",SUBCAT_DIR,Nsnap);
		myfopen(fp,buf,"w");
		load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
		load_particle_data(Nsnap,SNAPSHOT_DIR);
		fresh_ID2Index(&Cat,-1); 
		IDrhomax=mymalloc(sizeof(int)*Cat.Ngroups);
		IDpotmin=mymalloc(sizeof(int)*Cat.Ngroups);
		#ifdef HALO_PARA
		#pragma omp parallel for 
		#endif
		for(grpid=0;grpid<Cat.Ngroups;grpid++)
			get_CenID(grpid);

		fwrite(&Cat.Ngroups,sizeof(int),1,fp);
		fwrite(IDrhomax,sizeof(int),Cat.Ngroups,fp);
		fwrite(IDpotmin,sizeof(int),Cat.Ngroups,fp);
		fwrite(&Cat.Ngroups,sizeof(int),1,fp);
		fclose(fp);
		myfree(IDpotmin);
		myfree(IDrhomax);
		free_catalogue(&Cat);
	}
	
return 0;
}

void get_CenID(int grpid)
{
	int i,*PIndex,irhomax,ipotmin;
	double *rho,*pot,potmin,rhomax;
	float hguess;
	
	PIndex=Cat.PIDorIndex+Cat.Offset[grpid];
	fprintf(logfile,"calculating densities ...\n");fflush(logfile);
	tree_tree_allocate(TREE_ALLOC_FACTOR*Cat.Len[grpid],Cat.Len[grpid]);
	maketree(Cat.Len[grpid],PIndex,Pdat.Pos);
	hguess=guess_ngb_range_halo(1);
	potmin=0.;ipotmin=0;rhomax=0.;irhomax=0;
	rho=mymalloc(sizeof(double)*Cat.Len[grpid]);
	pot=mymalloc(sizeof(double)*Cat.Len[grpid]);
	#pragma omp parallel for private(i) firstprivate(hguess) schedule(dynamic)
	for(i=0;i<Cat.Len[grpid];i++)
	{
		//~ clock_t T[3];
		//~ T[0]=clock();
		rho[i]=sph_density(Pdat.Pos[PIndex[i]],&hguess,PIndex,Pdat.Pos);
		//~ T[1]=clock();
		pot[i]=tree_treeevaluate_potential(Pdat.Pos[PIndex[i]],PIndex,Pdat.Pos);
		//~ T[2]=clock();
		//~ printf("%d:%f,%f\n",i,(double)(T[1]-T[0]),(double)(T[2]-T[1]));
		//~ printf("%ld\n", CLOCKS_PER_SEC);
	}
	for(i=0;i<Cat.Len[grpid];i++)
	{
		if(potmin>pot[i])
		{
		ipotmin=i;
		potmin=pot[i];
		}
		if(rhomax<rho[i])
		{
		irhomax=i;
		rhomax=rho[i];
		}
	}
	tree_tree_free();
	#ifndef PID_ORDERED
	IDrhomax[grpid]=Pdat.PID[PIndex[irhomax]];
	IDpotmin[grpid]=Pdat.PID[PIndex[ipotmin]];
	#else
	IDrhomax[grpid]=PIndex[irhomax]+1;
	IDpotmin[grpid]=PIndex[ipotmin]+1;
	#endif
}

void load_haloCenID(int *IDrhomax, int *IDpotmin, int Nsnap,int Ngroups)
{  //note these loaded are IDs, need to fresh to index before use.
	int N;
	char buf[1024];
	FILE *fp;
	
	if(Ngroups)
	if(!IDrhomax||!IDpotmin)
	{
		printf("error: allocate IDrhomax and IDpotmin first!\n");
		exit(1);
	}
	sprintf(buf,"%s/profile/haloCenID_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	fread(&N,sizeof(int),1,fp);
	if(N!=Ngroups)
	{
		printf("error loading %s\n 1: %d!=%d\n", buf, N, Ngroups);
		exit(1);
	}
	fread(IDrhomax,sizeof(int),Ngroups,fp);
	fread(IDpotmin,sizeof(int),Ngroups,fp);
	fread(&N,sizeof(int),1,fp);
	if(N!=Ngroups)
	{
		printf("error loading %s\n 2: %d!=%d\n", buf, N, Ngroups);
		exit(2);
	}
	fclose(fp);
}
