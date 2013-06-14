#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "globals.c"
#include "mymath.c"
#include "load_group.c"

int main()
{
SUBCATALOGUE SubCatA,SubCatB;
char dir1[1024]="/SANdisk5/kambrain/Sim6702/SubCat3";
char dir2[1024]="/SANdisk5/kambrain/Sim6702/SubCat3/bk";
int i,j,subid;	
	load_sub_catalogue(4,&SubCatA,dir1);
	load_sub_catalogue(4,&SubCatB,dir2);
	
	for(i=0;i<SubCatA.Nsubs;i++)
	{
		if(SubCatA.HostPortion[i]-SubCatB.HostPortion[i])
			printf("0: %d,%f,%f\n",i,SubCatA.HostPortion[i],SubCatB.HostPortion[i]);
	}
	for(i=0;i<SubCatA.Nsubs;i++)
	{
		if(SubCatA.HostWght[i]-SubCatB.HostWght[i])
			printf("1: %d,%f,%f\n",i,SubCatA.HostWght[i],SubCatB.HostWght[i]);
	}
	for(i=0;i<SubCatA.Nsubs;i++)
	{
		for(j=0;j<3;j++)
		{
		if(SubCatA.CoM[i][j]-SubCatB.CoM[i][j])
			printf("(%d: %d),%f,%f\n",i,j,SubCatA.CoM[i][j],SubCatB.CoM[i][j]);
		}
	}	
	printf("%d,%d\n",SubCatA.ProSubID[266],SubCatA.ProSubID[267]);
	printf("%f,%f,%f\n",SubCatA.CoM[0][0],SubCatA.CoM[0][1],SubCatA.CoM[0][2]);
	printf("%f,%f,%f\n",SubCatB.CoM[0][0],SubCatB.CoM[0][1],SubCatB.CoM[0][2]);

	if(SubCatA.Ngroups-SubCatB.Ngroups)
			 printf("ngroups %d,%d\n",SubCatA.SubLen[i],SubCatB.SubLen[i]);
	if(SubCatA.Nsubs-SubCatB.Nsubs)
			 printf("nsubs %d,%d\n",SubCatA.SubLen[i],SubCatB.SubLen[i]);
	if(SubCatA.Nbirth-SubCatB.Nbirth)
			 printf("nbirth %d,%d\n",SubCatA.SubLen[i],SubCatB.SubLen[i]);
	if(SubCatA.Ndeath-SubCatB.Ndeath)
			 printf("ndeath %d,%d\n",SubCatA.SubLen[i],SubCatB.SubLen[i]);
	if(SubCatA.NQuasi-SubCatB.NQuasi)
			 printf("nqua %d,%d\n",SubCatA.SubLen[i],SubCatB.SubLen[i]);
	if(SubCatA.Nsplitter-SubCatB.Nsplitter)
			 printf("nsp %d,%d\n",SubCatA.SubLen[i],SubCatB.SubLen[i]);
	for(subid=0;subid<SubCatA.Nsubs;subid++)
	{
		for(i=0;i<SubCatA.SubLen[subid];i++)
		{
		if(SubCatA.PSubArr[subid][i]-SubCatB.PSubArr[subid][i])
			printf("(%d:%d),%d,%d\n",subid,i,SubCatA.PSubArr[subid][i],SubCatB.PSubArr[subid][i]);
		}
	}
return 0;
}