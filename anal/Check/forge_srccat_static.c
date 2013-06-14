//to forge a srccat for the major merger simulations
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"


int main()
{
	char buf[1024];
	FILE *fp;
	SRCCATALOGUE SrcCat;
	SUBCATALOGUE SubCat;
	HBTInt *PIDs,i,len0,len1;
	float tmp;
	
	logfile=stdout;
	SrcCat.Nsubs=2;
	create_src_cat(&SrcCat);
	load_particle_data(0,SNAPSHOT_DIR);
	SrcCat.PSubArr[0]=mymalloc(sizeof(HBTInt)*NP_DM);
	SrcCat.PSubArr[1]=mymalloc(sizeof(HBTInt)*NP_DM);
	for(i=0,len0=0,len1=0;i<NP_DM;i++)
	{
		if(Pdat.Pos[i][0]<55.)
		{
		SrcCat.PSubArr[0][len0]=Pdat.PID[i];
		len0++;
		}
		else
		{
		SrcCat.PSubArr[1][len1]=Pdat.PID[i];
		len1++;	
		}
	}
	SrcCat.SubLen[0]=len0;
	SrcCat.SubLen[1]=len1;
	SrcCat.SubLen2[0]=0;
	SrcCat.SubLen2[1]=0;
	SrcCat.SubOffset[0]=0;
	SrcCat.SubOffset[1]=len0;
	SrcCat.CoreFrac[0]=CoreFrac0;
	SrcCat.CoreFrac[1]=CoreFrac0;
	SrcCat.Nids=NP_DM;
	SrcCat.NDeathSp=0;
	printf("%d,%d:%d\n",len0,len1,NP_DM);
	
	save_src_catalogue(0,&SrcCat,SUBCAT_DIR);
	
	SubCat.Ngroups=2;
	SubCat.GrpOffset_Sub=mymalloc(sizeof(HBTInt)*2);
	SubCat.GrpLen_Sub=mymalloc(sizeof(HBTInt)*2);
	SubCat.Nsubs=2;	create_sub_cat(&SubCat);
	for(i=0;i<2;i++)
	{
	SubCat.GrpLen_Sub[i]=1;
	SubCat.GrpOffset_Sub[i]=i;
	SubCat.SubLen[i]=0;
	SubCat.SubOffset[i]=0;
	SubCat.SubRank[i]=0;
	SubCat.HaloChains[i].HostID=i;	
	SubCat.HaloChains[i].ProSubID=-1;	
	SubCat.sub_hierarchy[i].nibs=-1;
	SubCat.sub_hierarchy[i].pre=-1;
	SubCat.sub_hierarchy[i].next=-1;
	SubCat.sub_hierarchy[i].sub=-1;
	SubCat.NQuasi=0;
	SubCat.Nsplitter=0;
	}
	save_sub_catalogue(0,&SubCat,SUBCAT_DIR);
	return 0;
}
