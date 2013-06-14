#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

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

  return 0;
}
int main(int argc,char **argv)
{
CATALOGUE Cat;
struct GroupData GrpData;

HBTInt Nsnap=0;
HBTInt i;
HBTInt GrpTags;
HBTReal b,r;

logfile=stdout;//redirect BT routines' log info to standard output

Nsnap=atoi(argv[1]);
b=atof(argv[2]);

load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
load_particle_data(Nsnap,SNAPSHOT_DIR);
fill_PIDHash();
fresh_ID2Index(&Cat,-1);
free_PIDHash();

r=b*pow(MP_DM/(3*HUBBLE0*HUBBLE0/8./M_PI/G*header.Omega0),1.0/3.0);

tree_tree_allocate(TREE_ALLOC_FACTOR*Cat.Len[0],Cat.Len[0]);
maketree(Cat.Len[0],Cat.PIDorIndex,Pdat.Pos);
GrpData.Np=Cat.Len[0];
GrpData.GrpTags=mymalloc(sizeof(struct ParticleGroup)*Cat.Len[0]);
for(i=0;i<Cat.Len[0];i++) GrpData.GrpTags[i].PIndex=Cat.PIDorIndex[i];
treesearch_linkgrp(r,Pdat.Pos,&GrpData);
tree_tree_free();
GrpLen=GrpData.GrpLen;
qsort(GrpData.GrpTags,GrpData.Np,sizeof(struct ParticleGroup),comp_grplen);
qsort(GrpLen,GrpData.Ngrp,sizeof(HBTInt), comp_int);
printf("%d Groups, %d MaxLen, %d, %d \n", GrpData.Ngrp, GrpData.GrpLen[0],GrpData.GrpLen[1],GrpData.GrpLen[2]);
return 0;
}
