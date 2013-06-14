/*extending the history to associate host halos all along each history;
 * for switch of host halos during cross, associate hosts according 
 * to distances scaled by host Rvir,choosing as the host to which the relative
 * distance is smaller*/
//to be done: check the case for HostID<0 when extending history between crossing
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
#include "gas_vars.h"
#include "proto.h"
#include "gas_proto.h"
#include "history_vars.h"
#include "history_proto.h"
typedef float TYPE_XYZ[3];
typedef struct
{
	int mass;//fof mass
	int Mvir[3];
	float Rvir[3];//[tophat,c200,b200],comoving
	int flag_badvir[3];
	int flag_fakehalo;//set to 1 when halo is not self-bound,to 0 otherwise.
} HALOSIZE;
struct ShardParam
{
	int SnapBirth;
	int SnapTidal;
	int SnapRvir;
	float Mrate[2];//at Rtidal, Rvir
	float Mhost[2];
	float Rhost[2];
	float Kappa[2];
	float j2[2];
	float ErrK[2];//error estimation due to discrete output
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
void load_halo_size(HALOSIZE *halosize,int Ngroups,int Nsnap)
{
	int i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,13*4L,SEEK_CUR);
		fread(&halosize[i].mass,sizeof(int),1,fp);
		fread(halosize[i].Mvir,sizeof(int),3,fp);
		fread(halosize[i].Rvir,sizeof(float),3,fp);
		fread(halosize[i].flag_badvir,sizeof(int),3,fp);
		fread(&halosize[i].flag_fakehalo,sizeof(int),1,fp);
		if(halosize[i].flag_badvir[0])
		{
			halosize[i].Rvir[0]=comoving_virial_radius(halosize[i].mass);
			halosize[i].Mvir[0]=halosize[i].mass;
		}
	}
	fclose(fp);
}
int load_halo_concentration(float *halocon,int Nsnap)
{
	int grpid,Ngroups,Ngroups2,halostatus;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_param_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	fread(&Ngroups,sizeof(int),1,fp);
	for(grpid=0;grpid<Ngroups;grpid++)
	{
		fread(&halostatus,sizeof(int),1,fp);
		fseek(fp,4*sizeof(float),SEEK_CUR);
		//~ fread(halostatus+grpid,sizeof(int),1,fp);
		fread(halocon+grpid,sizeof(float),1,fp);
		if(halostatus!=0) halocon[grpid]=-1; //set concentration to -1 if fitting not successful
		fseek(fp,sizeof(float)*2,SEEK_CUR);
	}	
	fread(&Ngroups2,sizeof(int),1,fp);
	if(Ngroups2!=Ngroups)
	{
		printf("error:Ngroups=%d,%d do not match when loading \n %s\n" 
			"probably file corruption or different file format\n",
			Ngroups,Ngroups2,buf);
		exit(-1);
	}
	fclose(fp);
	return Ngroups;
}
int read_mainsubid(int Nsnap,int hostid);
int read_Ngroups(int Nsnap);
int get_newinfall(int SnapBegin,int HistID);
int read_subpos(int Nsnap,TYPE_XYZ *(*snappos));
int read_subvel(int Nsnap,TYPE_XYZ *(*snapvel));
void history_copy_host(SubNode HistoryTo[],SubNode HistoryFrom[],int Nsnap);
void extend_history(int HistID,int SnapInfall);
void create_historyshards(struct HistoryShards *HistoryShard,int NumHist);
void split_history(int HistID);
void save_historyshards(struct HistoryShards *HistoryShard);
void get_shard_param(int HistID);

struct HistoryShards HistoryRevShard;
SubHist *ClusterHistory,*ClusterHistoryRev;
HALOSIZE *(halosize[MaxSnap]);
float *(halocon[MaxSnap]);
TYPE_XYZ *(SubCoM[MaxSnap]),*(SubVCoM[MaxSnap]);
int * (Sub2Hist[MaxSnap]);
float scaleF[MaxSnap];
float PartMass;
int main(int argc,char **argv)
{
//~ SUBCATALOGUE SubCat;
//~ CATALOGUE Cat;
//~ GASSUBCAT GSubCat;

//~ char gasdir[1024]=GASCAT_DIR;
char buf[1024];FILE *fp;
int Nsnap,i,SnapInfall,Nsubs,Ngroups;
int NumHist,HostID,HistID,HHistID;

logfile=stdout;

load_history(&NumHist,&ClusterHistory,SUBCAT_DIR);
load_history(&NumHist,&ClusterHistoryRev,SUBCAT_DIR);
create_historyshards(&HistoryRevShard,NumHist);
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
Ngroups=read_Ngroups(Nsnap);
halosize[Nsnap]=mymalloc(sizeof(HALOSIZE)*Ngroups);
load_halo_size(halosize[Nsnap],Ngroups,Nsnap);
halocon[Nsnap]=mymalloc(sizeof(float)*Ngroups);
load_halo_concentration(halocon[Nsnap],Nsnap);
read_subpos(Nsnap,SubCoM+Nsnap);
read_subvel(Nsnap,SubVCoM+Nsnap);
load_particle_header(Nsnap,SNAPSHOT_DIR);
scaleF[Nsnap]=header.time;
}
PartMass=header.mass[1];

for(HistID=0;HistID<NumHist;HistID++)
{
	SnapInfall=ClusterHistory[HistID].SnapEnter+1;
	if(SnapInfall<MaxSnap)
	{
	HostID=ClusterHistory[HistID].Member[SnapInfall].HostID;
	if(HostID<0)
	{
		printf("error: hostid<0,quasi-halo as infall?\n");
		exit(1);
	}
	else
	{
	HHistID=Sub2Hist[SnapInfall][read_mainsubid(SnapInfall,HostID)];
	for(Nsnap=ClusterHistory[HistID].SnapBirth;Nsnap<SnapInfall;Nsnap++)//extend pre-infall part
		history_copy_host(ClusterHistoryRev[HistID].Member,ClusterHistory[HHistID].Member,Nsnap);
	extend_history(HistID,SnapInfall);//cover gaps during branch-crossing
	}
	}
	split_history(HistID);	
	get_shard_param(HistID);
}

  sprintf(buf, "%s/history/ClusterHistoryRev", SUBCAT_DIR);
  myfopen(fp,buf,"w");
  fwrite(&NumHist,sizeof(int),1,fp);
  fwrite(ClusterHistoryRev,sizeof(SubHist),NumHist,fp);
  fclose(fp);	
  
  save_historyshards(&HistoryRevShard);	

return 0;
}
void get_sat_concen(int HistID)
{
	int ShardID,Nsnap,FlagRvir,FlagTidal,HostID;
	float concen;
	if(HistoryRevShard.NBirth[HistID])
	{
	concen=0;
	ShardID=0;
	FlagTidal=FlagRvir=0;
	for(Nsnap=ClusterHistory[HistID].SnapBirth;Nsnap<MaxSnap;Nsnap++)
	{
		if(0==ClusterHistory[HistID].Member[Nsnap].SubRank)//update concen every time being a halo
		{
		HostID=ClusterHistory[HistID].Member[Nsnap].HostID;//use original history to get real host
		concen=halocon[Nsnap][HostID];
		}
		if(Nsnap==HistoryRevShard.Par[HistID][ShardID].SnapTidal)
		{
			HistoryRevShard.Par[HistID][ShardID].Csat[0]=concen;
			FlagTidal=1;//finished
		}
		if(Nsnap==HistoryRevShard.Par[HistID][ShardID].SnapRvir)
		{
			HistoryRevShard.Par[HistID][ShardID].Csat[1]=concen;
			FlagRvir=1;//finished
		}
		if(FlagTidal&&FlagRvir)//both finished
		{
			ShardID++;//new shard
			FlagTidal=0;
			FlagRvir=0;
			if(ShardID==HistoryRevShard.NBirth[HistID])//all done
			break;
		}
	}
	}
}
void get_shard_param(int HistID)
{
	struct ShardParam *ShardPar;
	int i,ShardID,Nsnap,SnapEnd,FlagRvir,FlagTidal,SubID,MainID,HostID;
	float Rup,Rdown,Rtidal,Tup,Tdown,r,rup,v,Pos[3],Vel[3],Kt,Kup,Kdown,j2up,j2down,vc2,sint2;
	for(ShardID=0;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)
	{
		ShardPar=HistoryRevShard.Par[HistID]+ShardID;
		if(ShardID==HistoryRevShard.NBirth[HistID]-1)//the last piece
		SnapEnd=MaxSnap;
		else
		SnapEnd=HistoryRevShard.Par[HistID][ShardID+1].SnapBirth;
		Rdown=2;
		Tdown=2;
		r=0;
		Kdown=0;
		j2down=0;
		FlagRvir=0;FlagTidal=0;
		for(Nsnap=ShardPar->SnapBirth;Nsnap<SnapEnd;Nsnap++)
		{
			SubID=ClusterHistoryRev[HistID].Member[Nsnap].SubID;
			HostID=ClusterHistoryRev[HistID].Member[Nsnap].HostID;
			if(HostID<0||SubID<0)//host or sub died, for last shard
			break;
			MainID=read_mainsubid(Nsnap,HostID);
			if(MainID==SubID)//self-host
			break;
			rup=r;
			r=distance(SubCoM[Nsnap][SubID],SubCoM[Nsnap][MainID]);
			v=distance(SubVCoM[Nsnap][SubID],SubVCoM[Nsnap][MainID]);
			for(i=0;i<3;i++)
			{
				Pos[i]=SubCoM[Nsnap][SubID][i]-SubCoM[Nsnap][MainID][i];
				Vel[i]=SubVCoM[Nsnap][SubID][i]-SubVCoM[Nsnap][MainID][i];
			}
			Rup=Rdown;
			Rdown=r/halosize[Nsnap][HostID].Rvir[0];
			Kup=Kdown;
			j2up=j2down;
			vc2=G*(ClusterHistoryRev[HistID].Member[Nsnap].Mhost+ClusterHistoryRev[HistID].Member[Nsnap].Mdm)*PartMass/r/scaleF[Nsnap];
			Kdown=v*v/vc2;
			sint2=1.0-pow(f_prod(Pos,Vel,3)/r/v,2);
			Kt=Kdown*sint2;
			j2down=Kt*(2-Kdown);
			if(FlagRvir==0)
			{
				if(Rdown<1)
				{
				FlagRvir=1;	
				ShardPar->SnapRvir=Nsnap;
				if(Nsnap==ShardPar->SnapBirth)
				{
				ShardPar->Mhost[1]=ClusterHistoryRev[HistID].Member[Nsnap].Mhost;
				ShardPar->Chost[1]=ClusterHistoryRev[HistID].Member[Nsnap].Chost;
				ShardPar->Mrate[1]=ClusterHistoryRev[HistID].Member[Nsnap].Mdm/ShardPar->Mhost[1];									
				ShardPar->Rhost[1]=r;
				}
				else
				{
				ShardPar->Mhost[1]=(ClusterHistoryRev[HistID].Member[Nsnap-1].Mhost-ClusterHistoryRev[HistID].Member[Nsnap].Mhost)
									*(1.0-Rdown)/(Rup-Rdown)+ClusterHistoryRev[HistID].Member[Nsnap].Mhost;
				if(ClusterHistoryRev[HistID].Member[Nsnap-1].Chost<0)//bad previous
				ShardPar->Chost[1]=ClusterHistoryRev[HistID].Member[Nsnap].Chost;
				else if(ClusterHistoryRev[HistID].Member[Nsnap].Chost<0)//bad current
				ShardPar->Chost[1]=ClusterHistoryRev[HistID].Member[Nsnap-1].Chost;
				else//both good
				ShardPar->Chost[1]=(ClusterHistoryRev[HistID].Member[Nsnap-1].Chost-ClusterHistoryRev[HistID].Member[Nsnap].Chost)
									*(1.0-Rdown)/(Rup-Rdown)+ClusterHistoryRev[HistID].Member[Nsnap].Chost;									
				ShardPar->Mrate[1]=((ClusterHistoryRev[HistID].Member[Nsnap-1].Mdm-ClusterHistoryRev[HistID].Member[Nsnap].Mdm)
									*(1.0-Rdown)/(Rup-Rdown)+ClusterHistoryRev[HistID].Member[Nsnap].Mdm)/ShardPar->Mhost[1];									
				ShardPar->Rhost[1]=(rup-r)*(1.0-Rdown)/(Rup-Rdown)+r;
				//~ ShardPar->Kappa[1]=(Kup-Kdown)*(1-Rdown)/(Rup-Rdown)+Kdown;
				//~ ShardPar->j2[1]=(j2up-j2down)*(1-Rdown)/(Rup-Rdown)+j2down;
				//~ ShardPar->Kappa[1]=(Rup-1<1-Rdown)?Kup:Kdown;
				//~ ShardPar->j2[1]=(Rup-1<1-Rdown)?j2up:j2down;
				//~ ShardPar->ErrK[1]=sqrt((ShardPar->Kappa[1]-Kdown)*(ShardPar->Kappa[1]-Kdown)+(ShardPar->Kappa[1]-Kup)*(ShardPar->Kappa[1]-Kup));
				//~ ShardPar->Errj2[1]=sqrt((ShardPar->j2[1]-j2down)*(ShardPar->j2[1]-j2down)+(ShardPar->j2[1]-j2up)*(ShardPar->j2[1]-j2up));
				}
				ShardPar->Kappa[1]=Kdown;
				ShardPar->j2[1]=j2down;
				ShardPar->ErrK[1]=Kup;
				ShardPar->Errj2[1]=j2up;
				}
			}
			if(FlagTidal==0)
			{
				Tup=Tdown;
				Tdown=Rdown/pow(2+Kt,1./3)/1.25;
				if(Tdown<1)
				{
				FlagTidal=1;
				if(Nsnap==ShardPar->SnapBirth)
				{
				Rtidal=Rdown;
				ShardPar->Mhost[0]=ClusterHistoryRev[HistID].Member[Nsnap].Mhost;
				ShardPar->Chost[0]=ClusterHistoryRev[HistID].Member[Nsnap].Chost;
				ShardPar->Mrate[0]=ClusterHistoryRev[HistID].Member[Nsnap].Mdm/ShardPar->Mhost[0];									
				ShardPar->Rhost[0]=r;
				}
				else
				{
				Rtidal=(Rup-Rdown)*(1.0-Tdown)/(Tup-Tdown)+Rdown;//r/Rhost when T=1
				ShardPar->Mhost[0]=(ClusterHistoryRev[HistID].Member[Nsnap-1].Mhost-ClusterHistoryRev[HistID].Member[Nsnap].Mhost)
									*(Rtidal-Rdown)/(Rup-Rdown)+ClusterHistoryRev[HistID].Member[Nsnap].Mhost;
				if(ClusterHistoryRev[HistID].Member[Nsnap-1].Chost<0)//bad previous
				ShardPar->Chost[1]=ClusterHistoryRev[HistID].Member[Nsnap].Chost;
				else if(ClusterHistoryRev[HistID].Member[Nsnap].Chost<0)//bad current
				ShardPar->Chost[1]=ClusterHistoryRev[HistID].Member[Nsnap-1].Chost;
				else//both good
				ShardPar->Chost[0]=(ClusterHistoryRev[HistID].Member[Nsnap-1].Chost-ClusterHistoryRev[HistID].Member[Nsnap].Chost)
									*(Rtidal-Rdown)/(Rup-Rdown)+ClusterHistoryRev[HistID].Member[Nsnap].Chost;									
				ShardPar->Mrate[0]=((ClusterHistoryRev[HistID].Member[Nsnap-1].Mdm-ClusterHistoryRev[HistID].Member[Nsnap].Mdm)
									*(Rtidal-Rdown)/(Rup-Rdown)+ClusterHistoryRev[HistID].Member[Nsnap].Mdm)/ShardPar->Mhost[0];
				ShardPar->Rhost[0]=(rup-r)*(Rtidal-Rdown)/(Rup-Rdown)+r;
				//~ ShardPar->Kappa[0]=(Kup-Kdown)*(Rtidal-Rdown)/(Rup-Rdown)+Kdown;
				//~ ShardPar->j2[0]=(j2up-j2down)*(Rtidal-Rdown)/(Rup-Rdown)+j2down;
				//~ ShardPar->Kappa[0]=(Rup-Rtidal<Rtidal-Rdown)?Kup:Kdown;;
				//~ ShardPar->j2[0]=(Rup-Rtidal<Rtidal-Rdown)?j2up:j2down;
				}
				if((Nsnap>ShardPar->SnapBirth)&&
				(ClusterHistoryRev[HistID].Member[Nsnap-1].Mdm>ClusterHistoryRev[HistID].Member[Nsnap].Mdm))
				{
				ShardPar->SnapTidal=Nsnap-1;
				ShardPar->Kappa[0]=Kup;
				ShardPar->j2[0]=j2up;
				ShardPar->ErrK[0]=Kdown;
				ShardPar->Errj2[0]=j2down;
				}
				else
				{
				ShardPar->SnapTidal=Nsnap;
				ShardPar->Kappa[0]=Kdown;;
				ShardPar->j2[0]=j2down;
				ShardPar->ErrK[0]=Kup;
				ShardPar->Errj2[0]=j2up;
				}
				//~ ShardPar->ErrK[0]=sqrt((ShardPar->Kappa[0]-Kdown)*(ShardPar->Kappa[0]-Kdown)+(ShardPar->Kappa[0]-Kup)*(ShardPar->Kappa[0]-Kup));
				//~ ShardPar->Errj2[0]=sqrt((ShardPar->j2[0]-j2down)*(ShardPar->j2[0]-j2down)+(ShardPar->j2[0]-j2up)*(ShardPar->j2[0]-j2up));
				}
			}
			if(FlagRvir&&FlagTidal)
			break;
		}
			if(FlagRvir==0)
			{
				ShardPar->SnapRvir=-1;
				ShardPar->Mhost[1]=0;
				ShardPar->Chost[1]=0;
				ShardPar->Rhost[1]=0;
				ShardPar->Kappa[1]=0;
				ShardPar->j2[1]=0;
				ShardPar->ErrK[1]=0;
				ShardPar->Errj2[1]=0;
			}
			if(FlagTidal==0)
			{
				ShardPar->SnapTidal=-1;
				ShardPar->Mhost[0]=0;
				ShardPar->Chost[0]=0;
				ShardPar->Rhost[0]=0;
				ShardPar->Kappa[0]=0;
				ShardPar->j2[0]=0;
				ShardPar->ErrK[0]=0;
				ShardPar->Errj2[0]=0;
			}
	}
	get_sat_concen(HistID);
}
void split_history(int HistID)//this act on ClusterHistoryRev
{
	int NBirth_cache,SnapBirth_cache[MaxSnap];
	int i,Nsnap,HostID_Up,HHistID_Up,HostID,HHistID;
	
	Nsnap=ClusterHistoryRev[HistID].SnapEnter+1;
	if(Nsnap<MaxSnap)
	{
		NBirth_cache=1;
		for(i=ClusterHistoryRev[HistID].SnapBirth;i<Nsnap;i++)
			if(ClusterHistoryRev[HistID].Member[i].HostID>=0)
				break;
		SnapBirth_cache[0]=i;//The first snapshot when sub and its host both exist,
							//just in case the host could be born later than sub
		HostID_Up=ClusterHistoryRev[HistID].Member[Nsnap].HostID;
		HHistID_Up=Sub2Hist[Nsnap][read_mainsubid(Nsnap,HostID_Up)];
		for(Nsnap=Nsnap+1;Nsnap<MaxSnap;Nsnap++)
		{
			HostID=ClusterHistoryRev[HistID].Member[Nsnap].HostID;//get host from rev
			if(ClusterHistoryRev[HistID].Member[Nsnap].SubID<0)//history ends
			break;
			if(ClusterHistory[HHistID_Up].Member[Nsnap].HostID!=HostID)//different dest-branch,compare to original host history(not extended)
			{
				if(HostID<0) 
				//~ {printf("error: HostID<0 not expected here\n");exit(1);}
				break;//host history ends
				HHistID=Sub2Hist[Nsnap][read_mainsubid(Nsnap,HostID)];
				if(ClusterHistory[HHistID].Member[Nsnap-1].HostID!=HostID_Up)//also different pro-branch,branch-crossing is happenning,compare to original host history(not extended)
				{
					SnapBirth_cache[NBirth_cache]=Nsnap;
					NBirth_cache++;
				}
				HHistID_Up=HHistID;//need to update HistID	
			}
			HostID_Up=HostID;//always need to update HostID
		}
		HistoryRevShard.NBirth[HistID]=NBirth_cache;
		HistoryRevShard.Par[HistID]=mymalloc(sizeof(struct ShardParam)*NBirth_cache);
		for(i=0;i<NBirth_cache;i++)
			HistoryRevShard.Par[HistID][i].SnapBirth=SnapBirth_cache[i];
	}
	else
	{
		HistoryRevShard.NBirth[HistID]=0;
		HistoryRevShard.Par[HistID]=NULL;
	}
}
void extend_history(int HistID,int SnapInfall)
{
	int i,Nsnap,HHistID_Up,HHistID_Down,HostID_Up,HostID_Down;
	float Rup,Rdown;
	Nsnap=SnapInfall;
	while(Nsnap<MaxSnap)
	{
		if(ClusterHistory[HistID].Member[Nsnap].SubRank>0)//still staying in previous host
		{
			HostID_Up=ClusterHistory[HistID].Member[Nsnap].HostID;
			HHistID_Up=Sub2Hist[Nsnap][read_mainsubid(Nsnap,HostID_Up)];
			Nsnap++;
		}
		else
		{
			SnapInfall=get_newinfall(Nsnap,HistID);
			if(SnapInfall==MaxSnap)//no subsequent infall
			{
				while(Nsnap<MaxSnap)
				{
				history_copy_host(ClusterHistoryRev[HistID].Member,ClusterHistory[HHistID_Up].Member,Nsnap);
				Nsnap++;
				}
				break;
			}
			else//extend Host info to fill the gap
			{
				HostID_Down=ClusterHistory[HistID].Member[SnapInfall].HostID;
				HHistID_Down=Sub2Hist[SnapInfall][read_mainsubid(SnapInfall,HostID_Down)];
				while(Nsnap<SnapInfall)
				{
					HostID_Up=ClusterHistory[HHistID_Up].Member[Nsnap].HostID;
					HostID_Down=ClusterHistory[HHistID_Down].Member[Nsnap].HostID;		
					if(HostID_Up>=0&&HostID_Down>=0)
					{
						Rup=distance(SubCoM[Nsnap][read_mainsubid(Nsnap,HostID_Up)],SubCoM[Nsnap][ClusterHistory[HistID].Member[Nsnap].SubID]);
						Rdown=distance(SubCoM[Nsnap][read_mainsubid(Nsnap,HostID_Down)],SubCoM[Nsnap][ClusterHistory[HistID].Member[Nsnap].SubID]);
						if(!(halosize[Nsnap][HostID_Up].flag_badvir[0])&&!(halosize[Nsnap][HostID_Down].flag_badvir[0]))
						{
							Rup=Rup/halosize[Nsnap][HostID_Up].Rvir[0];
							Rdown=Rdown/halosize[Nsnap][HostID_Down].Rvir[0];
						}
						else
						{
							Rup=Rup/pow(halosize[Nsnap][HostID_Up].mass,1./3);
							Rdown=Rdown/pow(halosize[Nsnap][HostID_Down].mass,1./3);
						}
						if(Rup<Rdown)
						history_copy_host(ClusterHistoryRev[HistID].Member,ClusterHistory[HHistID_Up].Member,Nsnap);	
						else
						history_copy_host(ClusterHistoryRev[HistID].Member,ClusterHistory[HHistID_Down].Member,Nsnap);	
					}
					else if(HostID_Up<0&&HostID_Down>=0)
						history_copy_host(ClusterHistoryRev[HistID].Member,ClusterHistory[HHistID_Down].Member,Nsnap);	
					else if(HostID_Up>=0&&HostID_Down<0)
						history_copy_host(ClusterHistoryRev[HistID].Member,ClusterHistory[HHistID_Up].Member,Nsnap);
					//else,HostID_Up<0&&HostID_Down<0,do nothing,leave the host info as is,i.e.,the host is itself
					Nsnap++;
				}
			}
		}
	}
}
void create_historyshards(struct HistoryShards *HistoryShard,int NumHist)
{
HistoryShard->NumHist=NumHist;
HistoryShard->NBirth=mymalloc(sizeof(int)*NumHist);
HistoryShard->HistoryOffset=mymalloc(sizeof(int)*NumHist);
HistoryShard->Par=mymalloc(sizeof(struct ShardParam *)*NumHist);
}
void save_historyshards(struct HistoryShards *HistoryShard)
{
int NumHist,i,Offset;
char buf[1024];
FILE *fp;
  
  NumHist=HistoryShard->NumHist;
  Offset=0;	
  for(i=0;i<NumHist;i++)
  {
  HistoryShard->HistoryOffset[i]=Offset;
  Offset+=HistoryShard->NBirth[i];
  }
  HistoryShard->NumShards=Offset;
  	
  sprintf(buf, "%s/history/HistoryShards", SUBCAT_DIR);
  myfopen(fp,buf,"w");
  fwrite(&NumHist,sizeof(int),1,fp);
  fwrite(&Offset,sizeof(int),1,fp);
  fwrite(HistoryShard->NBirth,sizeof(int),NumHist,fp);
  fwrite(HistoryShard->HistoryOffset,sizeof(int),NumHist,fp);
  for(i=0;i<NumHist;i++)
  fwrite(HistoryShard->Par[i],sizeof(struct ShardParam),HistoryShard->NBirth[i],fp);
  fclose(fp);	
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
void history_copy_host(SubNode HistoryTo[MaxSnap],SubNode HistoryFrom[MaxSnap],int Nsnap)
{//copy hostid, Mhost and Chost from HistoryFrom to HistoryTo at time Nsnap
HistoryTo[Nsnap].HostID=HistoryFrom[Nsnap].HostID;
HistoryTo[Nsnap].Mhost=HistoryFrom[Nsnap].Mhost;
HistoryTo[Nsnap].Chost=HistoryFrom[Nsnap].Chost;
}
int read_mainsubid(int Nsnap,int hostid)
{
FILE *fd;
char buf[1024];
int subid,Ngroups;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, Nsnap);
  myfopen(fd,buf,"r");
  
  fread(&Ngroups, sizeof(int), 1, fd);
  if(hostid<0||hostid>=Ngroups)
  {
	  printf("error: wrong hostid to read for mainsub\n");
	  exit(1);
  }
  fseek(fd,sizeof(int)*(2+Ngroups+hostid),SEEK_CUR);
  fread(&subid,sizeof(int),1,fd);
  fclose(fd);
  return subid;
}
int read_Ngroups(int Nsnap)
{
FILE *fd;
char buf[1024];
int Ngroups;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(int), 1, fd);
  
  fclose(fd);
  return Ngroups;
}
int read_subpos(int Nsnap,TYPE_XYZ *(*snappos))
{
FILE *fd;
char buf[1024];
int subid,Ngroups,Nsubs;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nsubs,sizeof(int),1,fd);
  *snappos=mymalloc(sizeof(float)*3*Nsubs);
  fseek(fd,sizeof(int)*(1+Ngroups*2+Nsubs*3),SEEK_CUR);
  fseek(fd,(sizeof(struct Chain_data)+sizeof(struct Hierarchy))*Nsubs,SEEK_CUR);
  for(subid=0;subid<Nsubs;subid++)
  {
  fread(*snappos+subid,sizeof(float),3,fd);
  fseek(fd,sizeof(float)*8,SEEK_CUR);
  }
  fclose(fd);
  return Nsubs;
}
int read_subvel(int Nsnap,TYPE_XYZ *(*snapvel))
{
FILE *fd;
char buf[1024];
int subid,Ngroups,Nsubs;

  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, Nsnap);
  myfopen(fd,buf,"r");

  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&Nsubs,sizeof(int),1,fd);
  *snapvel=mymalloc(sizeof(float)*3*Nsubs);
  fseek(fd,sizeof(int)*(1+Ngroups*2+Nsubs*3),SEEK_CUR);
  fseek(fd,(sizeof(struct Chain_data)+sizeof(struct Hierarchy))*Nsubs,SEEK_CUR);
  fseek(fd,sizeof(float)*3,SEEK_CUR);
  for(subid=0;subid<Nsubs;subid++)
  {
  fread(*snapvel+subid,sizeof(float),3,fd);
  fseek(fd,sizeof(float)*8,SEEK_CUR);
  }
  fclose(fd);
  return Nsubs;
}
int get_newinfall(int SnapBegin,int HistID)
{
	int i;
	for(i=SnapBegin;i<MaxSnap;i++)
	{
		if(ClusterHistory[HistID].Member[i].SubRank>0)
		break;
	}
	return i;
}
