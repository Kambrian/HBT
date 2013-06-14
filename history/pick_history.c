/*to pick out the history for a specified subhalo*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

//~ #define GAS_TRACE
#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#ifdef GAS_TRACE
#include "gas_vars.h"
#include "gas_proto.h"
#endif
#include "history_vars.h"
#include "history_proto.h"


struct HistoryShards HistoryRevShard;
EVOLUTIONCAT EvoCat;
HBTInt * (Sub2Hist[MaxSnap]);

HBTInt sum_SubInSub_mass(HBTInt subid,SUBCATALOGUE *SubCat)
{
	HBTInt id,mass;
	mass=SubCat->SubLen[subid];
	id=SubCat->sub_hierarchy[subid].sub;
	while(id>=0) //add-up sub-in-sub mass
	{
	mass+=sum_SubInSub_mass(id,SubCat);
	id=SubCat->sub_hierarchy[id].next;
	}
	return mass;
}
#ifdef GAS_TRACE
HBTReal U_mean_sub(HBTInt subid,GASSUBCAT *GSubCat)
{
	HBTInt i,pid;
	HBTReal Um;
	Um=0;
	if(0==GSubCat->SubLen[subid]) return Um;
	
	#pragma omp parallel for private(i,pid) reduction(+:Um)
	for(i=0;i<GSubCat->SubLen[subid];i++)
	{
		pid=GSubCat->PSubArr[subid][i];
		Um+=Gdat.U[pid];
	}
	Um/=GSubCat->SubLen[subid];
	return Um;
}
#endif
int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
#ifdef GAS_TRACE
//~ GASHALOCAT GCat;
GASSUBCAT GSubCat;
#endif

char buf[1024];FILE *fp;
HBTInt Nsnap,i,SnapInfall,SnapTidal,SnapBirth,SnapEnd,SnapDeath;
HBTInt SnapLoad,subid;

HBTInt NumHist,NumShards,HistID,ShardID,Nsubs;

if(argc!=3)
{
printf("usage:%s [Nsnap] [subid]\n",argv[0]);
exit(1);
}
SnapLoad=atoi(argv[1]);
subid=atoi(argv[2]);
logfile=stdout;

load_evocat_rev(&EvoCat,SUBCAT_DIR);
NumHist=EvoCat.NHist;
load_historyshards(&HistoryRevShard);
NumShards=HistoryRevShard.NumShards;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
HistID=Sub2Hist[SnapLoad][subid];

sprintf(buf,"%s/anal/history_%03d_"HBTIFMT".txt",SUBCAT_DIR,(int)SnapLoad,subid);
myfopen(fp,buf,"w");
SnapBirth=EvoCat.History[HistID].SnapBirth;
SnapDeath=EvoCat.History[HistID].SnapDeath;
fprintf(fp,"Nsnap,sublen,hostlen,subrank,r,r/Rvir\n");
HBTReal partmass;
load_particle_header(0,SNAPSHOT_DIR);
//~ partmass=header.mass[1];
partmass=1;
SubNode *Member;
for(Nsnap=SnapBirth;Nsnap<SnapDeath;Nsnap++)
{
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	Member=GetMember(&EvoCat,HistID,Nsnap);
	HBTReal r,rn;
	if(Member->HostID>=0)
	{
	r=distance(SubCat.Property[Member->SubID].CoM,SubCat.Property[SubCat.GrpOffset_Sub[Member->HostID]].CoM);
	rn=r/comoving_virial_radius(Member->Mhost);
	}
	else
	{
	r=0;
	rn=0;	
	}
	fprintf(fp,""HBTIFMT",%g,%g,"HBTIFMT",%g,%g\n",Nsnap,Member->Mdm*partmass,Member->Mhost*partmass,Member->SubRank,r,rn);
	free_sub_table(&SubCat);
}
fclose(fp);
return 0;
}


