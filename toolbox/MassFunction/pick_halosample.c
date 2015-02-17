//to pick up a sample of halos 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>


#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "history_vars.h"
#include "history_proto.h"

#define NSAMPLE 30
#define NMASS 3   //bin number for Mhost



EVOLUTIONCAT EvoCat;
HALOSIZE *halosize;
float PartMass;
void load_evocat_raw(EVOLUTIONCAT *EvoCat)
{
	load_history_pre(EvoCat,SUBCAT_DIR);
	EvoCat->PartMass=header.mass[1];;
}
void select_groups(int *grpid,int Ngroups)
{
	int i,j,k,l;
	int Mvir[NMASS],mvir[NMASS];
	Mvir[0]=3e5/PartMass;
	Mvir[1]=1.08e4/PartMass;
	Mvir[2]=1.05e3/PartMass;
	//~ Mvir[3]=1.02e2/PartMass;
	mvir[0]=0.2e5/PartMass;
	mvir[1]=0.92e4/PartMass;
	mvir[2]=0.95e3/PartMass;
	//~ mvir[3]=0.98e2/PartMass;
	//~ int Mvir[3],mvir[3];
	//~ Mvir[0]=3e4/PartMass;
	//~ Mvir[1]=1.2e3/PartMass;
	//~ Mvir[2]=1.1e2/PartMass;
	//~ mvir[0]=0.3e4/PartMass;
	//~ mvir[1]=0.8e3/PartMass;
	//~ mvir[2]=0.9e2/PartMass;
	j=0;k=0;l=0;
	for(i=0;i<Ngroups;i++)
	{
		if(halosize[i].flag_badvir[0]||
		   halosize[i].flag_fakehalo||
		   halosize[i].mass>halosize[i].Mvir[0]*1.5)
		continue;
		if(halosize[i].Mvir[0]<Mvir[j]&&halosize[i].Mvir[0]>mvir[j])
		{
			grpid[k]=i;
			printf("%.1e,",halosize[i].Mvir[0]*PartMass);fflush(stdout);
			k++;
			l++;
			if(l==NSAMPLE/NMASS)
			{
				printf("\n");
				l=0;
				j++;
				if(j==NMASS) break;
			}
		}
	}
	if(j<NMASS)
	{
		printf("halos not enough\n");
		exit(1);
	}
}

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
char buf[1024];FILE *fp;
int grpid[NSAMPLE],mainid[NSAMPLE];
int Ngroups,Nsnap,Nsubs;
int * Sub2Hist,i;

logfile=stdout;

Nsnap=MaxSnap-1;
Ngroups=read_Ngroups(Nsnap);
halosize=mymalloc(sizeof(HALOSIZE)*Ngroups);
load_halo_size(halosize,Ngroups,Nsnap);
load_particle_header(Nsnap,SNAPSHOT_DIR);
PartMass=header.mass[1];
select_groups(grpid,Ngroups);
for(i=0;i<NSAMPLE;i++)
	mainid[i]=read_mainsubid(MaxSnap-1,grpid[i]);
printf("loading histories...\n");fflush(stdout);	
load_evocat_raw(&EvoCat);
printf("loading sub2hist...\n");fflush(stdout);
load_sub2hist(MaxSnap-1,&Sub2Hist,&Nsubs,SUBCAT_DIR);
printf("writing...\n");fflush(stdout);
int Hid[NSAMPLE];
for(i=0;i<NSAMPLE;i++)
Hid[i]=Sub2Hist[mainid[i]];
sprintf(buf,"%s/anal/massfun/samplehaloid.txt",SUBCAT_DIR);
myfopen(fp,buf,"w");
int HostID;
for(Nsnap=MaxSnap-1;Nsnap>=0;Nsnap--)
{
	for(i=0;i<NSAMPLE;i++)
	{
		HostID=GetHostID(GetMember(&EvoCat,Hid[i],Nsnap));
		fprintf(fp,"%d,",HostID);
	}
	fprintf(fp,"\n");
}
fclose(fp);
return 0;
}

