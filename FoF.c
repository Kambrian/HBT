#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define GRPLENMIN NBOUNDMIN


extern void build_group_catalogue(CATALOGUE *Cat, HBTReal linklength);

int main(int argc,char **argv)
{
CATALOGUE Cat;
HBTInt Nsnap;
HBTReal b,r;

mkdir(GRPCAT_DIR,0755);
logfile=stdout;//redirect BT routines' log info to standard output

Nsnap=atoi(argv[1]);
//b=atof(argv[2]);
b=0.2; //linking length

load_particle_data_bypart(Nsnap,SNAPSHOT_DIR,FLAG_LOAD_POS);//only load position
#ifdef PERIODIC_BDR
r=b*BOXSIZE/cbrt(NP_DM);
#elif NP_DM==NP_GAS||NP_GAS==0
r=b*pow((MP_DM+MP_GAS)/(3*HUBBLE0*HUBBLE0/8./M_PI/G*header.Omega0),1.0/3.0);
#else
#error	 "cannot automatically determine linking length"
#endif
build_group_catalogue(&Cat,r);
free_particle_data();

fprintf(logfile,"saving groups ...\n");fflush(stdout);
load_particle_data_bypart(Nsnap,SNAPSHOT_DIR,FLAG_LOAD_ID);//only load id
save_group_catalogue_HBT(Nsnap,&Cat,GRPCAT_DIR);
free_particle_data();

return 0;
}

static int comp_int(const void *a, const void *b)//in descending order
{
  if(*((HBTInt *) a) > *((HBTInt *)b))
    return -1;

  if(*((HBTInt *) a) < *((HBTInt *)b))
    return +1;

  return 0;
}
static HBTInt *GrpLen;
static int comp_grplen(const void *a, const void *b)//in descending order
{
  if(GrpLen[((struct ParticleGroup *) a)->GrpID] > GrpLen[((struct ParticleGroup *) b)->GrpID])
    return -1;

  if(GrpLen[((struct ParticleGroup *) a)->GrpID] < GrpLen[((struct ParticleGroup *) b)->GrpID])
    return +1;  
//equal grplen, sort in ascending order of grpid
  if(((struct ParticleGroup *) a)->GrpID > ((struct ParticleGroup*) b)->GrpID)
    return 1;
  
  if(((struct ParticleGroup *) a)->GrpID < ((struct ParticleGroup*) b)->GrpID)
    return -1;

  return 0;
}

void build_group_catalogue(CATALOGUE *Cat, HBTReal linklength)
{
HBTInt *PIndex, i;	
struct GroupData GrpData;

PIndex=mymalloc(sizeof(HBTInt)*NP_DM);
for(i=0;i<NP_DM;i++) PIndex[i]=i;
GrpData.Np=NP_DM;
treesearch_linkgrp(linklength,PIndex,&GrpData);
myfree(PIndex);

GrpLen=GrpData.GrpLen;
qsort(GrpData.GrpTags,GrpData.Np,sizeof(struct ParticleGroup),comp_grplen);
qsort(GrpLen,GrpData.Ngrp,sizeof(HBTInt), comp_int);
fprintf(logfile," MaxLen: "HBTIFMT", "HBTIFMT", "HBTIFMT" ... \n", GrpData.GrpLen[0], GrpData.GrpLen[1], GrpData.GrpLen[2]);

for(i=0;i<GrpData.Ngrp;i++)
	if(GrpLen[i]<GRPLENMIN) break;
Cat->Ngroups=i;
fprintf(logfile,HBTIFMT" good groups\n", Cat->Ngroups);fflush(logfile);
Cat->Len=realloc(GrpLen,sizeof(HBTInt)*Cat->Ngroups);
Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
HBTInt Offset=0;
for(i=0;i<Cat->Ngroups;i++)
{
	Cat->Offset[i]=Offset;
	Offset+=GrpLen[i];
}
Cat->Nids=Offset;
Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
for(i=0;i<Cat->Nids;i++)
	Cat->PIDorIndex[i]=GrpData.GrpTags[i].PIndex;

myfree(GrpData.GrpTags);	
}
