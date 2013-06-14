/* To find destination for those subs who just die
 * 
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define CHECK_PID(pid)	if(pid<=NP_GAS||pid>NP_GAS+NP_DM) {printf("wrong pid:%d, out of range (%d,%d]\n",pid,NP_GAS,NP_GAS+NP_DM);	exit(1);}

int main(int argc,char **argv)
{
SUBCATALOGUE SubCatA,SubCatB;
int *pro2dest,Nsubs,*id2sub_base,*id2sub;


if(argc!=3)
{
printf("usage:%s [Mmin] [Mmax]\n",argv[0]);
exit(1);
}

id2sub_base=mymalloc(sizeof(int)*NP_DM);
id2sub=id2sub_base-1-NP_GAS;//pid range: NP_GAS+1 ~ NP_GAS+NP_DM

for(SnapshotNum=0;SnapshotNum<MaxSnap-1;SnapshotNum++)
{
	load_pro2dest(SnapshotNum,&pro2dest,&Nsubs,SUBCAT_DIR);
	load_sub_catalogue(SnapshotNum,&SubCatA,SUBCAT_DIR);
	load_sub_catalogue(SnapshotNum,&SubCatB,SUBCAT_DIR);
	for(i=0;i<NP_DM;i++)
	{
		id2sub_base[i]=-1;
	}
	for(subid=0;subid<SubCatB.Nsubs;subid++)
	{
		for(i=0;i<SubCatB.SubLen[subid];i++)
		{
			pid=SubCatB.PSubArr[subid][i];
			CHECK_PID(pid)
			id2sub[pid]=subid;
		}
	}
	for(subid=0;subid<SubCatA.Nsubs;subid++)
	{
		if(pro2dest[subid]<0)
		{
			for(i=0;i<SubCatA.SubLen[subid];i++)
			{
				pid=SubCatA.PSubArr[subid][i];
				CHECK_PID(pid)
				
			}
		}
		
	}
	
}

return 0;
}

int prepare_ind2halo(CATALOGUE *A)
{int i,haloid,pid;
//	A->ID2Halo=mymalloc(sizeof(int)*NP_DM);

	#pragma omp parallel 
		{
	#pragma omp for
	for(i=0;i<NP_DM;i++)//initialization
	{
		A->ID2Halo[i]=-1;/*return -1 if the PID does not belong to a Halo,
								i.e,we consider the backgroud as a halo with haloid=-1; 
								note that this only make sense when we try to find host for bound structures */
	}
	#pragma omp for private(i,pid,haloid)
	for(haloid=0;haloid<A->Ngroups;haloid++)
	{
		for(i=0;i<A->Len[haloid];i++)
		{
			pid=A->PIDorIndex[A->Offset[haloid]+i];//Pindex ranges [0,NP_DM);
			A->ID2Halo[pid]=haloid;//haloIDs begins from id=0
		}
	}
		}
	return 1; // means id2halo ready;
}
