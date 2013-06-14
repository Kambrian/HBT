/*extending the history to associate host halos all along each history;
 * for switch of host halos during cross, associate hosts according 
 * to distances scaled by host Rvir,choosing as the host to which the relative
 * distance is smaller
 * split history into pieces according to host-crossing (mergers)
 * get parameters for each merger*/
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
//~ #include "gas_proto.h"
#include "history_vars.h"
#include "history_proto.h"

HBTInt get_newinfall(HBTInt SnapBegin,HBTInt HistID);
void history_copy_host(SubNode *MemberTo,SubNode *MemberFrom);
void extend_history(HBTInt HistID,HBTInt SnapInfall);
void split_history(HBTInt HistID);
void get_shard_param(HBTInt HistID);

struct HistoryShards HistoryRevShard;
EVOLUTIONCAT EvoCat,EvoCatRev;
HALOSIZE *(halosize[MaxSnap]);
HBTReal *(halocon[MaxSnap]);
HBTxyz *(SubCoM[MaxSnap]),*(SubVCoM[MaxSnap]);
HBTInt * (Sub2Hist[MaxSnap]);
HBTReal scaleF[MaxSnap];
HBTReal PartMass;
void load_evocat_raw(EVOLUTIONCAT *EvoCat)
{
	HBTInt HistID,Nsnap;
	SubNode *Member;
	load_history_pre(EvoCat,SUBCAT_DIR);
	EvoCat->PartMass=header.mass[1];;
	for(HistID=0;HistID<EvoCat->NHist;HistID++)
	{
		for(Nsnap=EvoCat->History[HistID].SnapBirth;Nsnap<EvoCat->History[HistID].SnapDeath;Nsnap++)
		{
			Member=EvoCat->History[HistID].Member-EvoCat->History[HistID].SnapBirth+Nsnap;
			HBTInt HostID;
			HostID=Member->HostID;
			if(HostID<0)//quasi halo
			{
			Member->Mhost=0;
			Member->Chost=0;
			}
			else
			{
			//~ Member->Mhost=Cat.Len[HostID];
			Member->Mhost=halosize[Nsnap][HostID].Mvir[0];
			Member->Chost=halocon[Nsnap][HostID];
			}
		}
	}
}

int main(int argc,char **argv)
{
//~ SUBCATALOGUE SubCat;
//~ CATALOGUE Cat;
//~ GASSUBCAT GSubCat;

//~ char gasdir[1024]=GASCAT_DIR;
char buf[1024];FILE *fp;
HBTInt Nsnap,i,MemID,SnapInfall,Nsubs,Ngroups;
HBTInt NumHist,HostID,HistID,HHistID;

logfile=stdout;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
Ngroups=read_Ngroups(Nsnap);
halosize[Nsnap]=mymalloc(sizeof(HALOSIZE)*Ngroups);
load_halo_size(halosize[Nsnap],Ngroups,Nsnap);
halocon[Nsnap]=mymalloc(sizeof(HBTReal)*Ngroups);
load_halo_concentration(halocon[Nsnap],Nsnap);
read_subpos(Nsnap,SubCoM+Nsnap);
read_subvel(Nsnap,SubVCoM+Nsnap);
load_particle_header(Nsnap,SNAPSHOT_DIR);
scaleF[Nsnap]=header.time;
printf(HBTIFMT".",Nsnap);fflush(stdout);
}
printf("\n");
PartMass=header.mass[1];
load_evocat_raw(&EvoCat);
load_evocat_raw(&EvoCatRev);
NumHist=EvoCat.NHist;

create_historyshards(&HistoryRevShard,NumHist);

#pragma omp parallel for private(HistID,SnapInfall,HostID,HHistID,Nsnap)
for(HistID=0;HistID<NumHist;HistID++)
{
	//~ static HBTInt progress=0;
	//~ if(HistID>progress/100.0*NumHist)
	//~ {
		//~ printf("#");fflush(stdout);
		//~ progress++;
	//~ }
	HISTORY *History;
	History=EvoCat.History+HistID;
	SnapInfall=History->SnapEnter+1;
	if(SnapInfall<History->SnapDeath)
	{
	HostID=History->Member[SnapInfall-History->SnapBirth].HostID;
	if(HostID<0)
	{
		printf("error: hostid<0,quasi-halo as infall?\n");
		exit(1);
	}
	else
	{
	HHistID=Sub2Hist[SnapInfall][read_mainsubid(SnapInfall,HostID)];
	for(Nsnap=History->SnapBirth;Nsnap<SnapInfall;Nsnap++)//extend pre-infall part
		history_copy_host(GetMember(&EvoCatRev,HistID,Nsnap),GetMember(&EvoCat,HHistID,Nsnap));
	extend_history(HistID,SnapInfall);//cover gaps during branch-crossing
	}
	}
	split_history(HistID);	
	get_shard_param(HistID);
}
printf("\n");

  save_evocat_rev(&EvoCatRev,SUBCAT_DIR);  
  save_historyshards(&HistoryRevShard);	

return 0;
}

void get_sat_concen(HBTInt HistID)
{
	HBTInt ShardID,Nsnap,FlagRvir,FlagTidal,HostID;
	HBTReal concen;
	HBTInt SnapBirth,SnapDeath;
	SnapBirth=EvoCat.History[HistID].SnapBirth;
	SnapDeath=EvoCat.History[HistID].SnapDeath;
	if(HistoryRevShard.NBirth[HistID])
	{
	concen=0;
	ShardID=0;
	FlagTidal=FlagRvir=0;
	for(Nsnap=SnapBirth;Nsnap<SnapDeath;Nsnap++)
	{
		if(0==EvoCat.History[HistID].Member[Nsnap-SnapBirth].SubRank)//update concen every time being a halo
		{
		HostID=EvoCat.History[HistID].Member[Nsnap-SnapBirth].HostID;//use original history to get real host
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
void get_shard_param(HBTInt HistID)
{
	struct ShardParam *ShardPar;
	HBTInt i,ShardID,Nsnap,SnapEnd,FlagRvir,FlagTidal,SubID,MainID,HostID;
	HBTReal Rup,Rdown,Rtidal,Tup,Tdown,r,rup,v,Pos[3],Vel[3],Kt,Kup,Kdown,j2up,j2down,vc2,sint2;
	HBTInt SnapBirth,SnapDeath;
	SnapBirth=EvoCatRev.History[HistID].SnapBirth;
	SnapDeath=EvoCatRev.History[HistID].SnapDeath;
	for(ShardID=0;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)
	{
		ShardPar=HistoryRevShard.Par[HistID]+ShardID;
		if(ShardID==HistoryRevShard.NBirth[HistID]-1)//the last piece
		SnapEnd=SnapDeath;
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
			SubNode *MemberNow;
			MemberNow=GetMember(&EvoCatRev,HistID,Nsnap);
			SubID=MemberNow->SubID;
			HostID=MemberNow->HostID;
			if(HostID<0||SubID<0)//host or sub died, for last shard
			//~ {printf("error: HostID<0 or SubID<0 not expected here\n");exit(1);}
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
			vc2=G*(MemberNow->Mhost+MemberNow->Mdm)*PartMass/r/scaleF[Nsnap];
			Kdown=v*v/vc2;
			sint2=1.0-pow(vec_prod(Pos,Vel,3)/r/v,2);
			Kt=Kdown*sint2;
			j2down=Kt*(2-Kdown);
			if(FlagRvir==0)
			{
				if(Rdown<1)
				{
				FlagRvir=1;	
				ShardPar->SnapRvir=Nsnap;
				if(Nsnap==ShardPar->SnapBirth)//no previous snapshot
				{
				ShardPar->Mhost[1]=MemberNow->Mhost;
				ShardPar->Chost[1]=MemberNow->Chost;
				ShardPar->Mrate[1]=MemberNow->Mdm/ShardPar->Mhost[1];									
				ShardPar->Rhost[1]=r;
				}
				else
				{
				ShardPar->Mhost[1]=((MemberNow-1)->Mhost-MemberNow->Mhost)
									*(1.0-Rdown)/(Rup-Rdown)+MemberNow->Mhost;
				if((MemberNow-1)->Chost<0)//bad previous
				ShardPar->Chost[1]=MemberNow->Chost;
				else if(MemberNow->Chost<0)//bad current
				ShardPar->Chost[1]=(MemberNow-1)->Chost;
				else//both good
				ShardPar->Chost[1]=((MemberNow-1)->Chost-MemberNow->Chost)
									*(1.0-Rdown)/(Rup-Rdown)+MemberNow->Chost;									
				ShardPar->Mrate[1]=(((MemberNow-1)->Mdm-MemberNow->Mdm)
									*(1.0-Rdown)/(Rup-Rdown)+MemberNow->Mdm)/ShardPar->Mhost[1];									
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
				ShardPar->Mhost[0]=MemberNow->Mhost;
				ShardPar->Chost[0]=MemberNow->Chost;
				ShardPar->Mrate[0]=MemberNow->Mdm/ShardPar->Mhost[0];									
				ShardPar->Rhost[0]=r;
				}
				else
				{
				Rtidal=(Rup-Rdown)*(1.0-Tdown)/(Tup-Tdown)+Rdown;//r/Rhost when T=1
				ShardPar->Mhost[0]=((MemberNow-1)->Mhost-MemberNow->Mhost)
									*(Rtidal-Rdown)/(Rup-Rdown)+MemberNow->Mhost;
				if((MemberNow-1)->Chost<0)//bad previous
				ShardPar->Chost[0]=MemberNow->Chost;
				else if(MemberNow->Chost<0)//bad current
				ShardPar->Chost[0]=(MemberNow-1)->Chost;
				else//both good					
				ShardPar->Chost[0]=((MemberNow-1)->Chost-MemberNow->Chost)
									*(Rtidal-Rdown)/(Rup-Rdown)+MemberNow->Chost;									
				ShardPar->Mrate[0]=(((MemberNow-1)->Mdm-MemberNow->Mdm)
									*(Rtidal-Rdown)/(Rup-Rdown)+MemberNow->Mdm)/ShardPar->Mhost[0];
				ShardPar->Rhost[0]=(rup-r)*(Rtidal-Rdown)/(Rup-Rdown)+r;
				//~ ShardPar->Kappa[0]=(Kup-Kdown)*(Rtidal-Rdown)/(Rup-Rdown)+Kdown;
				//~ ShardPar->j2[0]=(j2up-j2down)*(Rtidal-Rdown)/(Rup-Rdown)+j2down;
				//~ ShardPar->Kappa[0]=(Rup-Rtidal<Rtidal-Rdown)?Kup:Kdown;;
				//~ ShardPar->j2[0]=(Rup-Rtidal<Rtidal-Rdown)?j2up:j2down;
				}
				if((Nsnap>ShardPar->SnapBirth)&&((MemberNow-1)->Mdm>MemberNow->Mdm))
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

void split_history(HBTInt HistID)//this act on EvoCatRev.History, split according to mergers
{
	HBTInt NBirth_cache,SnapBirth_cache[MaxSnap];
	HBTInt i,Nsnap,HostID_Up,HHistID_Up,HostID,HHistID;
	HBTInt SnapBirth,SnapDeath;
	SnapBirth=EvoCat.History[HistID].SnapBirth;
	SnapDeath=EvoCat.History[HistID].SnapDeath;
	Nsnap=EvoCatRev.History[HistID].SnapEnter+1;
	if(Nsnap<SnapDeath)
	{
		NBirth_cache=1;
		for(i=SnapBirth;i<Nsnap;i++)
			if(GetMember(&EvoCatRev,HistID,i)->HostID>=0)
				break;
		SnapBirth_cache[0]=i;//The first snapshot when sub and its host both exist,
							//just in case the host could be born later than sub
		HostID_Up=GetMember(&EvoCatRev,HistID,Nsnap)->HostID;
		HHistID_Up=Sub2Hist[Nsnap][read_mainsubid(Nsnap,HostID_Up)];
		for(Nsnap=Nsnap+1;Nsnap<SnapDeath;Nsnap++)
		{
			HostID=GetMember(&EvoCatRev,HistID,Nsnap)->HostID;//get host from rev
			if(GetMember(&EvoCatRev,HistID,Nsnap)->SubID<0)//history ends
			{printf("error: history ends before death:HistID="HBTIFMT",Nsnap="HBTIFMT"\n",HistID,Nsnap);exit(1);}
			if(GetHostID(GetMember(&EvoCat,HHistID_Up,Nsnap))!=HostID)//different dest-branch,compare to original host history(not extended)
			{
				if(HostID<0) 
				//~ {printf("error: HostID<0 not expected here:"HBTIFMT"\n",GetMember(&EvoCat,HistID,Nsnap)->HostID);exit(1);}
				break;//host history ends
				HHistID=Sub2Hist[Nsnap][read_mainsubid(Nsnap,HostID)];
				if(GetHostID(GetMember(&EvoCat,HHistID,Nsnap-1))!=HostID_Up)//also different pro-branch,branch-crossing is happenning,compare to original host history(not extended)
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
void extend_history(HBTInt HistID,HBTInt SnapInfall)
{
	HBTInt i,Nsnap,HHistID_Up,HHistID_Down,HostID_Up,HostID_Down;
	HBTInt SnapBirth,SnapDeath;
	HBTReal Rup,Rdown;
	Nsnap=SnapInfall;
	SnapBirth=EvoCat.History[HistID].SnapBirth;
	SnapDeath=EvoCat.History[HistID].SnapDeath;
	while(Nsnap<SnapDeath)
	{
		if(GetMember(&EvoCat,HistID,Nsnap)->SubRank>0)//still staying in previous host
		{
			HostID_Up=GetMember(&EvoCat,HistID,Nsnap)->HostID;
			HHistID_Up=Sub2Hist[Nsnap][read_mainsubid(Nsnap,HostID_Up)];
			Nsnap++;
		}
		else
		{
			SnapInfall=get_newinfall(Nsnap,HistID);
			if(SnapInfall==SnapDeath)//no subsequent infall
			{
				while(Nsnap<SnapDeath)
				{
				history_copy_host(GetMember(&EvoCatRev,HistID,Nsnap),GetMember(&EvoCat,HHistID_Up,Nsnap));
				Nsnap++;
				}
				break;
			}
			else//extend Host info to fill the gap
			{
				HostID_Down=GetMember(&EvoCat,HistID,SnapInfall)->HostID;
				HHistID_Down=Sub2Hist[SnapInfall][read_mainsubid(SnapInfall,HostID_Down)];
				while(Nsnap<SnapInfall)
				{
					HostID_Up=GetHostID(GetMember(&EvoCat,HHistID_Up,Nsnap));
					HostID_Down=GetHostID(GetMember(&EvoCat,HHistID_Down,Nsnap));		
					if(HostID_Up>=0&&HostID_Down>=0)
					{
						Rup=distance(SubCoM[Nsnap][read_mainsubid(Nsnap,HostID_Up)],SubCoM[Nsnap][GetMember(&EvoCat,HistID,Nsnap)->SubID]);
						Rdown=distance(SubCoM[Nsnap][read_mainsubid(Nsnap,HostID_Down)],SubCoM[Nsnap][GetMember(&EvoCat,HistID,Nsnap)->SubID]);
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
						history_copy_host(GetMember(&EvoCatRev,HistID,Nsnap),GetMember(&EvoCat,HHistID_Up,Nsnap));	
						else
						history_copy_host(GetMember(&EvoCatRev,HistID,Nsnap),GetMember(&EvoCat,HHistID_Down,Nsnap));	
					}
					else if(HostID_Up<0&&HostID_Down>=0)
						history_copy_host(GetMember(&EvoCatRev,HistID,Nsnap),GetMember(&EvoCat,HHistID_Down,Nsnap));	
					else if(HostID_Up>=0&&HostID_Down<0)
						history_copy_host(GetMember(&EvoCatRev,HistID,Nsnap),GetMember(&EvoCat,HHistID_Up,Nsnap));
					//else,HostID_Up<0&&HostID_Down<0,do nothing,leave the host info as is,i.e.,the host is itself
					Nsnap++;
				}
			}
		}
	}
}

void history_copy_host(SubNode *MemberTo,SubNode *MemberFrom)
{//copy hostid, Mhost and Chost from HistoryFrom to HistoryTo at time Nsnap
if(MemberTo==NULL)
{
	printf("error: bad SubNode destination to copy");
	exit(1);
}
if(MemberFrom==NULL)
{
MemberTo->HostID=-2;
MemberTo->Mhost=0;
MemberTo->Chost=0;
}
else
{
MemberTo->HostID=MemberFrom->HostID;
MemberTo->Mhost=MemberFrom->Mhost;
MemberTo->Chost=MemberFrom->Chost;
}
}

HBTInt get_newinfall(HBTInt SnapBegin,HBTInt HistID)
{
	HBTInt i;
	SubNode *Member;
	Member=EvoCat.History[HistID].Member-EvoCat.History[HistID].SnapBirth;
	for(i=SnapBegin;i<EvoCat.History[HistID].SnapDeath;i++)
	{
		if(Member[i].SubRank>0)
		break;
	}
	return i;
}
