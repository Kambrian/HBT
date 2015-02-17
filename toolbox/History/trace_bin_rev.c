/*selection rules:
 * exclude splitters
 * Mass ratio bin
 * Minimum mass at Dstrp
 * Minimum snapshot_num at Dstrp
 * directly infall to central cluster
 * Shard based history to cut for crossing
 * */

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
#include "gas_vars.h"
#include "proto.h"
#include "gas_proto.h"
#include "history_vars.h"
#include "history_proto.h"

#define Mass_Min_at_Accr 300
#define Snap_Min_at_Accr 0   //or 40 to count only z>1
//~ #define Mmin 0.01
//~ #define Mmax 0.1
struct SubNodeExp
{
	int HostID;
	int subid;
	int Mdm;
	int Mgas;
	int Mhost;//host halo mass
	int Mcen;//host subhalo mass
	int Mhostgas;//host halo gas mass
	int Msub;//subhalo total dm mass (including sub-in-subs);
	float CoM[3];//relative position,comoving
	float VCoM[3];//relative Vel,physical
	float Chost;//host halo concentration
	float Umsat;//sat average thermal energy (Temperature)
};
struct ShardParam
{
	int SnapBirth;
	int SnapTidal;
	int SnapRvir;
	float Mrate[2];
	float Mhost[2];
	float Rhost[2];
	float Kappa[2];
	float j2[2];
	float ErrK[2];
	float Errj2[2];
	float Chost[2];//host concentration
	float Csat[2];//sat concentration
};
struct HistoryShards
{
	int NumHist;
	int NumShards;
	int *NBirth;
	int *HistoryOffset;
	struct ShardParam **Par;
};

struct HistoryShards HistoryRevShard;
SubHist *ClusterHistory;
int * (Sub2Hist[MaxSnap]);

void load_history_rev(int *P2NumHist,SubHist **P2ClusterHistory,char *subcatdir)
{
FILE *fp;
char buf[1024];
int NumHist;
SubHist *ClusterHistory;

  sprintf(buf, "%s/history/ClusterHistoryRev", subcatdir);
  myfopen(fp,buf,"r");
  fread(&NumHist,sizeof(int),1,fp);
  ClusterHistory=mymalloc(sizeof(SubHist)*NumHist);
  fread(ClusterHistory,sizeof(SubHist),NumHist,fp);
  fclose(fp);
  *P2NumHist=NumHist;
  *P2ClusterHistory=ClusterHistory;		
}
void load_historyshards(struct HistoryShards *HistoryShard)
{
int NumHist,i,Offset;
char buf[1024];
FILE *fp;
  	
  sprintf(buf, "%s/history/HistoryShards", SUBCAT_DIR);
  myfopen(fp,buf,"r");
  fread(&NumHist,sizeof(int),1,fp);
  
  HistoryShard->NumHist=NumHist;
  HistoryShard->NBirth=mymalloc(sizeof(int)*NumHist);
  HistoryShard->HistoryOffset=mymalloc(sizeof(int)*NumHist);
  HistoryShard->Par=mymalloc(sizeof(struct ShardParam *)*NumHist);
  
  fread(&HistoryShard->NumShards,sizeof(int),1,fp);
  fread(HistoryShard->NBirth,sizeof(int),NumHist,fp);
  fread(HistoryShard->HistoryOffset,sizeof(int),NumHist,fp);
  for(i=0;i<NumHist;i++)
  {
	HistoryShard->Par[i]=mymalloc(sizeof(struct ShardParam)*HistoryShard->NBirth[i]);  
	fread(HistoryShard->Par[i],sizeof(struct ShardParam),HistoryShard->NBirth[i],fp);
  }
  fclose(fp);	
}			
int sum_SubInSub_mass(int subid,SUBCATALOGUE *SubCat)
{
	int id,mass;
	mass=SubCat->SubLen[subid];
	id=SubCat->sub_hierarchy[subid].sub;
	while(id>=0) //add-up sub-in-sub mass
	{
	mass+=sum_SubInSub_mass(id,SubCat);
	id=SubCat->sub_hierarchy[id].next;
	}
	return mass;
}
int read_mainsubid(int Nsnap,int hostid)
{
FILE *fd;
char buf[1024];
int subid,Ngroups;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(int), 1, fd);
  fseek(fd,sizeof(int)*(2+Ngroups+hostid),SEEK_CUR);
  fread(&subid,sizeof(int),1,fd);
  fclose(fd);
  return subid;
}
int is_main_family(int HistID,int ShardID,int level)//return 1 if this shard is connected to the central halo
{
	int SnapRvir,HostID;
	if(level>100) return -1;
	
	if(HistID==Sub2Hist[MaxSnap-1][0])	return 1;

	SnapRvir=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	if(SnapRvir<0||SnapRvir>=MaxSnap)	return 0; //no host
	
	HostID=ClusterHistory[HistID].Member[SnapRvir].HostID;
	HistID=Sub2Hist[SnapRvir][read_mainsubid(SnapRvir,HostID)];
	if(HistoryRevShard.NBirth[HistID]==0)//host have no host
	{
		if(HistID==Sub2Hist[MaxSnap-1][0])
		return 1;
		else
		return 0;
	}
	for(ShardID=1;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)
	{
		if(HistoryRevShard.Par[HistID][ShardID].SnapBirth>SnapRvir)
		break;
	}
	ShardID--;
	return is_main_family(HistID,ShardID,level+1);
}

int is_direct_main_family(int HistID,int ShardID)//return 1 if this shard falls directly to the central halo
{
	int SnapRvir,HostID;

	if(HistID==Sub2Hist[MaxSnap-1][0])	return 0;//do not count the main-branch itself

	SnapRvir=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	if(SnapRvir<0||SnapRvir>=MaxSnap)	return 0; //no host
	
	HostID=ClusterHistory[HistID].Member[SnapRvir].HostID;
	HistID=Sub2Hist[SnapRvir][read_mainsubid(SnapRvir,HostID)];
	if(HistID==Sub2Hist[MaxSnap-1][0])
	return 1;
	else
	return 0;
}
float U_mean_sub(int subid,GASSUBCAT *GSubCat)
{
	int i,pid;
	float Um;
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

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
//~ GASHALOCAT GCat;
GASSUBCAT GSubCat;
struct SubNodeExp (*History)[MaxSnap];

char buf[1024];FILE *fp;
int Nsnap,i,SnapInfall,SnapTidal,SnapBirth,SnapEnd;
float Mmin,Mmax,Mratio;//mass ratio at RTidal

int NumHist,NumShards,HistID,ShardID,subid,hostsubid,Nsubs;
int *HistList,*ShardList,NList,*NNode;
SubNode *OldNode;

if(argc!=3)
{
printf("usage:%s [Mmin] [Mmax]\n",argv[0]);
exit(1);
}
Mmin=atof(argv[1]);
Mmax=atof(argv[2]);
logfile=stdout;

load_history_rev(&NumHist,&ClusterHistory,SUBCAT_DIR);
load_historyshards(&HistoryRevShard);
NumShards=HistoryRevShard.NumShards;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
HistList=mymalloc(sizeof(int)*NumShards);
ShardList=mymalloc(sizeof(int)*NumShards);
NList=0;
for(HistID=0;HistID<NumHist;HistID++)
{
	if(HistoryRevShard.NBirth[HistID]&&ClusterHistory[HistID].ProHistID<0)//exclude splitters
	{
	//~ ShardID=HistoryRevShard.NBirth[HistID]-1;
	ShardID=0;
	//~ for(ShardID=0;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)
	//~ {
	SnapInfall=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	SnapTidal=HistoryRevShard.Par[HistID][ShardID].SnapTidal;
	if(SnapInfall>=Snap_Min_at_Accr&&SnapInfall<MaxSnap)
	{
	Mratio=HistoryRevShard.Par[HistID][ShardID].Mrate[0];
	if(Mratio>Mmin&&Mratio<Mmax&&ClusterHistory[HistID].Member[SnapTidal].Mdm>=Mass_Min_at_Accr)
	{
		if(is_direct_main_family(HistID,ShardID)==1)
		{
		HistList[NList]=HistID;
		ShardList[NList]=ShardID;
		NList++;
		}
	}
	}
	//~ }
	}
}
printf("NList=%d\n",NList);
printf("sizeof(SubNodeExp)=%ld\n",sizeof(struct SubNodeExp));
History=mymalloc(sizeof(struct SubNodeExp)*MaxSnap*NList);
NNode=calloc(NList,sizeof(int));
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	//~ load_gashalocat(Nsnap,&GCat,GASCAT_DIR);
	load_gassubcat(Nsnap,&GSubCat,GASCAT_DIR);
	load_gas_data(Nsnap,SNAPSHOT_DIR);
	#ifndef GRPINPUT_INDEX
	fresh_gasID2Index(&GSubCat,FRSH_SUBCAT); 	//fofcat of JING's data are originally PIndex rather than PID
	#endif
	for(i=0;i<NList;i++)
	{
		HistID=HistList[i];
		ShardID=ShardList[i];
		//~ SnapTidal=HistoryRevShard.Par[HistID][ShardID].SnapTidal;
		SnapBirth=HistoryRevShard.Par[HistID][ShardID].SnapBirth;
		if(ShardID==HistoryRevShard.NBirth[HistID]-1)//the last piece
		SnapEnd=MaxSnap;
		else
		SnapEnd=HistoryRevShard.Par[HistID][ShardID+1].SnapBirth;
		OldNode=ClusterHistory[HistID].Member+Nsnap;
		if((Nsnap>=SnapBirth)&&(Nsnap<SnapEnd)&&(OldNode->Mdm)&&(OldNode->HostID>=0))//also exclude those whose hosthalo haven't been bore
		{
			NNode[i]++;
			hostsubid=SubCat.GrpOffset_Sub[OldNode->HostID];
			subid=OldNode->SubID;
			History[i][Nsnap].HostID=OldNode->HostID;
			History[i][Nsnap].subid=subid;
			History[i][Nsnap].Mdm=OldNode->Mdm;
			History[i][Nsnap].Mgas=OldNode->Mgas;
			History[i][Nsnap].Mhost=OldNode->Mhost;
			History[i][Nsnap].Mcen=SubCat.SubLen[hostsubid];
			History[i][Nsnap].Mhostgas=GSubCat.SubLen[hostsubid];
			//~ History[i][Nsnap].Mhostgas=GCat.Len[OldNode->HostID];
			History[i][Nsnap].Msub=sum_SubInSub_mass(subid,&SubCat);
			History[i][Nsnap].CoM[0]=SubCat.Property[subid].CoM[0]-SubCat.Property[hostsubid].CoM[0];
			History[i][Nsnap].CoM[1]=SubCat.Property[subid].CoM[1]-SubCat.Property[hostsubid].CoM[1];
			History[i][Nsnap].CoM[2]=SubCat.Property[subid].CoM[2]-SubCat.Property[hostsubid].CoM[2];
			History[i][Nsnap].VCoM[0]=SubCat.Property[subid].VCoM[0]-SubCat.Property[hostsubid].VCoM[0];
			History[i][Nsnap].VCoM[1]=SubCat.Property[subid].VCoM[1]-SubCat.Property[hostsubid].VCoM[1];
			History[i][Nsnap].VCoM[2]=SubCat.Property[subid].VCoM[2]-SubCat.Property[hostsubid].VCoM[2];
			History[i][Nsnap].Chost=OldNode->Chost;
			History[i][Nsnap].Umsat=U_mean_sub(subid,&GSubCat);
		}
		else
		{
			History[i][Nsnap].Mdm=0;
		}
	}
	erase_sub_catalogue(&SubCat);
	//~ free_gashalocat(&GCat);
	erase_gassubcat(&GSubCat);
}
sprintf(buf,"%s/anal/history_%03d_%03d.dat",SUBCAT_DIR,(int)(100*Mmin),(int)(100*Mmax));
myfopen(fp,buf,"w");
fwrite(&NList,sizeof(int),1,fp);
for(i=0;i<NList;i++)
{
	fwrite(NNode+i,sizeof(int),1,fp);
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
		if(History[i][Nsnap].Mdm)
		{
		fwrite(&Nsnap,sizeof(int),1,fp);
		fwrite(History[i]+Nsnap,sizeof(struct SubNodeExp),1,fp);
		}
	}
}
fwrite(&NList,sizeof(int),1,fp);
fclose(fp);

sprintf(buf,"%s/anal/historypar_%03d_%03d.dat",SUBCAT_DIR,(int)(100*Mmin),(int)(100*Mmax));
myfopen(fp,buf,"w");
fwrite(&NList,sizeof(int),1,fp);
for(i=0;i<NList;i++)
{
	HistID=HistList[i];ShardID=ShardList[i];
	fwrite(HistoryRevShard.Par[HistID]+ShardID,sizeof(struct ShardParam),1,fp);
}
fwrite(&NList,sizeof(int),1,fp);
fclose(fp);

return 0;
}


