/*selection rules:
 * exclude splitters
 * Mass ratio bin
 * Minimum mass at Dstrp
 * Minimum snapshot_num at Dstrp
 * directly infall to central cluster
 * Shard based history to cut for crossing
 * */
#include "datatypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

//#define GAS_TRACE
#define MAIN_FAMILY

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#ifdef GAS_TRACE
#include "gas_vars.h"
#include "gas_proto.h"
#endif
#include "history_vars.h"
#include "history_proto.h"

#define Mass_Min_at_Accr 300
#define Snap_Min_at_Accr 0   //or 40 to count only z>1
//~ #define Mmin 0.01
//~ #define Mmax 0.1

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
HBTInt is_main_family(HBTInt HistID,HBTInt ShardID,HBTInt level)//return 1 if this shard is connected to the central halo
{
	HBTInt SnapRvir,HostID;
	if(level>100) return -1;  //for safety, deadlock check
	
	if(HistID==Sub2Hist[MaxSnap-1][0])	return 1;
	
	SnapRvir=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	if(SnapRvir<0)	return 0; //no host
	
	HostID=GetMember(&EvoCat,HistID,SnapRvir)->HostID;
	HistID=Sub2Hist[SnapRvir][read_mainsubid(SnapRvir,HostID)];
	if(HistoryRevShard.NBirth[HistID]==0)//host have no host
	{
		if(HistID==Sub2Hist[MaxSnap-1][0])
		return 1;
		else
		return 0;
	}
	for(ShardID=1;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)  //switch to host ShardID
	{
		if(HistoryRevShard.Par[HistID][ShardID].SnapBirth>SnapRvir)
		break;
	}
	ShardID--;
	return is_main_family(HistID,ShardID,level+1);
}

HBTInt is_direct_main_family(HBTInt HistID,HBTInt ShardID)//return 1 if this shard falls directly to the central halo
{
	HBTInt SnapRvir,HostID;

	if(HistID==Sub2Hist[MaxSnap-1][0])	return 0;//do not count the main-branch itself

	SnapRvir=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	if(SnapRvir<0)	return 0; //no host
	
	HostID=GetMember(&EvoCat,HistID,SnapRvir)->HostID;
	HistID=Sub2Hist[SnapRvir][read_mainsubid(SnapRvir,HostID)];
	if(HistID==Sub2Hist[MaxSnap-1][0])
	return 1;
	else
	return 0;
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
HBTInt main(HBTInt argc,char **argv)
{
SUBCATALOGUE SubCat;
#ifdef GAS_TRACE
//~ GASHALOCAT GCat;
GASSUBCAT GSubCat;
#endif
SubNodeExp (*History)[MaxSnap];

char buf[1024];FILE *fp;
HBTInt Nsnap,i,SnapInfall,SnapTidal,SnapBirth,SnapEnd,SnapDeath;
HBTReal Mmin,Mmax,Mratio;//mass ratio at RTidal

HBTInt NumHist,NumShards,HistID,ShardID,subid,hostsubid,Nsubs;
HBTInt *HistList,*ShardList,NList,*NNode;
SubNode *OldNode;

if(argc!=3)
{
printf("usage:%s [Mmin] [Mmax]\n",argv[0]);
exit(1);
}
Mmin=atof(argv[1]);
Mmax=atof(argv[2]);
logfile=stdout;

load_evocat_rev(&EvoCat,SUBCAT_DIR);
NumHist=EvoCat.NHist;
load_historyshards(&HistoryRevShard);
NumShards=HistoryRevShard.NumShards;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
HistList=mymalloc(sizeof(HBTInt)*NumShards);
ShardList=mymalloc(sizeof(HBTInt)*NumShards);
NList=0;
for(HistID=0;HistID<NumHist;HistID++)
{
	if(HistoryRevShard.NBirth[HistID]&&EvoCat.History[HistID].ProHistID<0)//exclude splitters
	{
	//~ ShardID=HistoryRevShard.NBirth[HistID]-1;
	ShardID=0;
	//~ for(ShardID=0;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)
	//~ {
	SnapInfall=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	SnapTidal=HistoryRevShard.Par[HistID][ShardID].SnapTidal;
	if(SnapInfall>=Snap_Min_at_Accr&&SnapInfall<MaxSnap&&SnapTidal>=0) //SnapInfall==-1 if no infall
	{
	Mratio=HistoryRevShard.Par[HistID][ShardID].Mrate[0];
	if(Mratio>Mmin&&Mratio<Mmax&&GetMember(&EvoCat,HistID,SnapTidal)->Mdm>=Mass_Min_at_Accr)
	{
		#ifdef MAIN_FAMILY
		if(is_direct_main_family(HistID,ShardID)==1)
		{
		#endif	
		HistList[NList]=HistID;
		ShardList[NList]=ShardID;
		NList++;
		#ifdef MAIN_FAMILY
		}
		#endif
	}
	}
	//~ }
	}
}
printf("NList="HBTIFMT"\n",NList);
printf("sizeof(SubNodeExp)=%ld\n",sizeof(SubNodeExp));
History=mymalloc(sizeof(SubNodeExp)*MaxSnap*NList);
NNode=calloc(NList,sizeof(HBTInt));
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	printf(""HBTIFMT".",Nsnap);fflush(stdout);
	#ifdef GAS_TRACE
	//~ load_gashalocat(Nsnap,&GCat,GASCAT_DIR);
	load_gassubcat(Nsnap,&GSubCat,GASCAT_DIR);
	load_gas_data(Nsnap,SNAPSHOT_DIR);
	#ifndef GRPINPUT_INDEX
	fresh_gasID2Index(&GSubCat,FRSH_SUBCAT); 	//fofcat of JING's data are originally PIndex rather than PID
	#endif
	#endif
	#pragma omp parallel for private(i,HistID,ShardID,SnapBirth,SnapDeath,SnapEnd,OldNode,hostsubid,subid)
	for(i=0;i<NList;i++)
	{
		HistID=HistList[i];
		ShardID=ShardList[i];
		//~ SnapTidal=HistoryRevShard.Par[HistID][ShardID].SnapTidal;
		SnapBirth=HistoryRevShard.Par[HistID][ShardID].SnapBirth;
		SnapDeath=EvoCat.History[HistID].SnapDeath;
		if(ShardID==HistoryRevShard.NBirth[HistID]-1)//the last piece
		SnapEnd=SnapDeath;
		else
		SnapEnd=HistoryRevShard.Par[HistID][ShardID+1].SnapBirth;
		OldNode=GetMember(&EvoCat,HistID,Nsnap);
		if((Nsnap>=SnapBirth)&&(Nsnap<SnapEnd)&&(OldNode->Mdm)&&(OldNode->HostID>=0))//also exclude those whose hosthalo haven't been bore
		{
			NNode[i]++;
			hostsubid=SubCat.GrpOffset_Sub[OldNode->HostID];
			subid=OldNode->SubID;
			History[i][Nsnap].HostID=OldNode->HostID;
			History[i][Nsnap].subid=subid;
			History[i][Nsnap].Mdm=OldNode->Mdm;
			History[i][Nsnap].Mhost=OldNode->Mhost;
			History[i][Nsnap].Mcen=SubCat.SubLen[hostsubid];
			History[i][Nsnap].Msub=sum_SubInSub_mass(subid,&SubCat);
			#ifdef GAS_TRACE
			History[i][Nsnap].Mgas=GSubCat.SubLen[subid];
			History[i][Nsnap].Mhostgas=GSubCat.SubLen[hostsubid];
			//~ History[i][Nsnap].Mhostgas=GCat.Len[OldNode->HostID];
			History[i][Nsnap].Umsat=U_mean_sub(subid,&GSubCat);
			#endif
			History[i][Nsnap].CoM[0]=SubCat.Property[subid].CoM[0]-SubCat.Property[hostsubid].CoM[0];
			History[i][Nsnap].CoM[1]=SubCat.Property[subid].CoM[1]-SubCat.Property[hostsubid].CoM[1];
			History[i][Nsnap].CoM[2]=SubCat.Property[subid].CoM[2]-SubCat.Property[hostsubid].CoM[2];
			History[i][Nsnap].VCoM[0]=SubCat.Property[subid].VCoM[0]-SubCat.Property[hostsubid].VCoM[0];
			History[i][Nsnap].VCoM[1]=SubCat.Property[subid].VCoM[1]-SubCat.Property[hostsubid].VCoM[1];
			History[i][Nsnap].VCoM[2]=SubCat.Property[subid].VCoM[2]-SubCat.Property[hostsubid].VCoM[2];
			History[i][Nsnap].Chost=OldNode->Chost;
		}
		else
		{
			History[i][Nsnap].Mdm=-1;
		}
	}
	free_sub_table(&SubCat);
	#ifdef GAS_TRACE
	//~ free_gashalocat(&GCat);
	erase_gassubcat(&GSubCat);
	#endif
}
printf("\n");
sprintf(buf,"%s/anal/history_%03d_%03d.dat",SUBCAT_DIR,(HBTInt)(100*Mmin),(HBTInt)(100*Mmax));
myfopen(fp,buf,"w");
fwrite(&EvoCat.PartMass,sizeof(HBTReal),1,fp);
HBTInt flag_gas;
#ifdef GAS_TRACE
flag_gas=1;
#else
flag_gas=0;
#endif
fwrite(&flag_gas,sizeof(HBTInt),1,fp);
fwrite(&NList,sizeof(HBTInt),1,fp);
for(i=0;i<NList;i++)
{
	fwrite(NNode+i,sizeof(HBTInt),1,fp);
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
		if(History[i][Nsnap].Mdm>0)
		{
		fwrite(&Nsnap,sizeof(HBTInt),1,fp);
		fwrite(History[i]+Nsnap,sizeof(SubNodeExp),1,fp);
		}
	}
}
fwrite(&NList,sizeof(HBTInt),1,fp);
fclose(fp);

sprintf(buf,"%s/anal/historypar_%03d_%03d.dat",SUBCAT_DIR,(HBTInt)(100*Mmin),(HBTInt)(100*Mmax));
myfopen(fp,buf,"w");
fwrite(&NList,sizeof(HBTInt),1,fp);
for(i=0;i<NList;i++)
{
	HistID=HistList[i];ShardID=ShardList[i];
	fwrite(HistoryRevShard.Par[HistID]+ShardID,sizeof(struct ShardParam),1,fp);
}
fwrite(&NList,sizeof(HBTInt),1,fp);
fclose(fp);

return 0;
}


