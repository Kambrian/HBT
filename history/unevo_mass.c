//unevoled halo mass fun, upon Yang's request
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
#include "proto.h"
#include "history_vars.h"
#include "history_proto.h"

#define NSAMPLE 400
#define Mass_Min_at_Accr 0
#define Snap_Min_at_Accr 0   //or 40 to count only z>1
//~ #define Mmin 0.01
//~ #define Mmax 0.1

struct HistoryShards HistoryRevShard;
HBTInt * (Sub2Hist[MaxSnap]);
HBTInt * (GrpOffset_Sub[MaxSnap]);
EVOLUTIONCAT EvoCatPre,EvoCat;
HALOSIZE *(halosize[MaxSnap]);
HBTReal scaleF[MaxSnap];
HBTReal PartMass;
HBTInt get_mainsubid(HBTInt Nsnap,HBTInt hostid)
{
	 if(hostid<0)
  {
	  printf("error: wrong hostid to read for mainsub\n");
	  exit(1);
  }
  return GrpOffset_Sub[Nsnap][hostid];
}
#define read_mainsubid get_mainsubid
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
			//~ Member->Chost=0;
			}
			else
			{
			//~ Member->Mhost=Cat.Len[HostID];
			Member->Mhost=halosize[Nsnap][HostID].Mvir[0];
			//~ Member->Chost=halocon[Nsnap][HostID];
			}
		}
	}
}

HBTInt is_grp_family(HBTInt grpmainid,HBTInt HistID,HBTInt ShardID,HBTInt level)//return 1 if this shard is connected to the central halo
{
	HBTInt SnapRvir,HostID;
	if(level>100) return -1;//for safety check of deadlock
	
	if(HistID==Sub2Hist[MaxSnap-1][grpmainid])	return 1;
	
	SnapRvir=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	if(SnapRvir<0)	return 0; //no host
	
	HostID=GetMember(&EvoCat,HistID,SnapRvir)->HostID;
	HistID=Sub2Hist[SnapRvir][read_mainsubid(SnapRvir,HostID)];
	if(HistoryRevShard.NBirth[HistID]==0)//host have no host
	{
		if(HistID==Sub2Hist[MaxSnap-1][grpmainid])
		return 1;
		else
		return 0;
	}
	for(ShardID=1;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)//find host shardID
	{
		if(HistoryRevShard.Par[HistID][ShardID].SnapBirth>SnapRvir)
		break;
	}
	ShardID--;
	return is_grp_family(grpmainid,HistID,ShardID,level+1);
}

HBTInt is_direct_grp_family(HBTInt grpmainid,HBTInt HistID,HBTInt ShardID)//return 1 if this shard falls directly to the central halo
{
	HBTInt SnapRvir,HostID;

	if(HistID==Sub2Hist[MaxSnap-1][grpmainid])	return 0;//do not count the main-branch itself

	SnapRvir=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	if(SnapRvir<0)	return 0; //no host
	
	HostID=GetMember(&EvoCat,HistID,SnapRvir)->HostID;
	HistID=Sub2Hist[SnapRvir][read_mainsubid(SnapRvir,HostID)];
	if(HistID==Sub2Hist[MaxSnap-1][grpmainid])
	return 1;
	else
	return 0;
}
void select_groups(HBTInt *grpid,HBTInt Ngroups)
{
	HBTInt i,j,k,l;
	HBTInt Mvir[4],mvir[4];
	Mvir[0]=3e5/PartMass;
	Mvir[1]=1.08e4/PartMass;
	Mvir[2]=1.05e3/PartMass;
	Mvir[3]=1.02e2/PartMass;
	mvir[0]=0.2e5/PartMass;
	mvir[1]=0.92e4/PartMass;
	mvir[2]=0.95e3/PartMass;
	mvir[3]=0.98e2/PartMass;
	//~ HBTInt Mvir[3],mvir[3];
	//~ Mvir[0]=3e4/PartMass;
	//~ Mvir[1]=1.2e3/PartMass;
	//~ Mvir[2]=1.1e2/PartMass;
	//~ mvir[0]=0.3e4/PartMass;
	//~ mvir[1]=0.8e3/PartMass;
	//~ mvir[2]=0.9e2/PartMass;
	j=0;k=0;l=0;
	for(i=0;i<Ngroups;i++)
	{
		if(halosize[MaxSnap-1][i].flag_badvir[0]||
		   halosize[MaxSnap-1][i].flag_fakehalo||
		   halosize[MaxSnap-1][i].mass>halosize[MaxSnap-1][i].Mvir[0]*1.5)
		continue;
		if(halosize[MaxSnap-1][i].Mvir[0]<Mvir[j]&&halosize[MaxSnap-1][i].Mvir[0]>mvir[j])
		{
			grpid[k]=i;
			printf("%.1e,",halosize[MaxSnap-1][i].Mvir[0]*PartMass);fflush(stdout);
			k++;
			l++;
			if(l==100)
			{
				printf("\n");
				l=0;
				j++;
				if(j==NSAMPLE/100) break;
			}
		}
	}
	if(j<NSAMPLE/100)
	{
		printf("halos not enough\n");
		exit(1);
	}
}


HBTInt get_snapinfall(HBTInt HistID,HBTInt ShardID)
{
	HBTInt i,SnapBirth,SnapDeath,SnapEnd;
	SubNode *Member;
	Member=EvoCatPre.History[HistID].Member-EvoCatPre.History[HistID].SnapBirth;
	SnapBirth=HistoryRevShard.Par[HistID][ShardID].SnapBirth;
	SnapDeath=EvoCatPre.History[HistID].SnapDeath;
	if(ShardID==HistoryRevShard.NBirth[HistID]-1)//the last piece
	SnapEnd=SnapDeath;
	else
	SnapEnd=HistoryRevShard.Par[HistID][ShardID+1].SnapBirth;
	for(i=SnapBirth;i<SnapEnd;i++)
	{
		if(Member[i].SubRank>0)
		break;
	}
	if(i==SnapEnd) 
	{
		printf("unexpected: Shard never enter? "HBTIFMT","HBTIFMT"\n",HistID,ShardID);
		//~ exit(1);
		return i-1;
	}
	if(i==SnapBirth)//infall since birth
	return i;
	
	return i-1;
}

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
char buf[1024];FILE *fp;
HBTInt Nsnap,i,SnapInfall,SnapTidal,SnapRvir,SnapBirth;
HBTInt grpid[NSAMPLE],mainid[NSAMPLE],GrpLen[NSAMPLE],*DirectInfall;

HBTInt NumHist,NumShards,HistID,ShardID,Nsubs,Ngroups;
HBTInt *HistList[NSAMPLE],*ShardList[NSAMPLE],NList;
HBTReal *Mhost,*Msub,*tinfall,*Mmax,*tmax,*MRvir,*tRvir;
SubNode *Node;

logfile=stdout;

//~ Nsnap=59;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
Ngroups=read_Ngroups(Nsnap);
halosize[Nsnap]=mymalloc(sizeof(HALOSIZE)*Ngroups);
load_halo_size(halosize[Nsnap],Ngroups,Nsnap);
//~ halocon[Nsnap]=mymalloc(sizeof(HBTReal)*Ngroups);
//~ load_halo_concentration(halocon[Nsnap],Nsnap);
load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
GrpOffset_Sub[Nsnap]=SubCat.GrpOffset_Sub;
SubCat.GrpOffset_Sub=NULL;
free_sub_table(&SubCat);
load_particle_header(Nsnap,SNAPSHOT_DIR);
scaleF[Nsnap]=header.time;
printf(""HBTIFMT".",Nsnap);fflush(stdout);
}
printf("\n");
PartMass=header.mass[1];
select_groups(grpid,Ngroups);
for(i=0;i<NSAMPLE;i++)
	mainid[i]=read_mainsubid(MaxSnap-1,grpid[i]);
load_evocat_raw(&EvoCatPre);
load_evocat_rev(&EvoCat,SUBCAT_DIR);
NumHist=EvoCat.NHist;
load_historyshards(&HistoryRevShard);
NumShards=HistoryRevShard.NumShards;

for(i=0;i<NSAMPLE;i++)
{
GrpLen[i]=0;
HistList[i]=mymalloc(sizeof(HBTInt)*(NumHist/10));
ShardList[i]=mymalloc(sizeof(HBTInt)*(NumHist/10));
}

//~ HBTInt progress=0;
#pragma omp parallel for private(HistID,ShardID,SnapInfall,SnapTidal,i)
for(HistID=0;HistID<NumHist;HistID++)
{
	//~ if(HistID>progress*NumHist/100) 
	//~ {
		//~ printf("#");fflush(stdout);
		//~ progress++;
	//~ }
	
	if(HistoryRevShard.NBirth[HistID]&&EvoCat.History[HistID].ProHistID<0)//exclude splitters
	{
	//~ ShardID=HistoryRevShard.NBirth[HistID]-1;
	//~ ShardID=0;
	for(ShardID=0;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)
	{
		SnapInfall=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
		SnapTidal=HistoryRevShard.Par[HistID][ShardID].SnapTidal;
		if(SnapInfall>=Snap_Min_at_Accr&&SnapInfall<MaxSnap&&SnapTidal>=0) //SnapInfall==-1 if no infall
		{
			if(GetMember(&EvoCatPre,HistID,SnapTidal)->Mdm>=Mass_Min_at_Accr)
			{
				for(i=0;i<NSAMPLE;i++)  //record this history to its corresponding list in the sample
				{
					if(is_grp_family(mainid[i],HistID,ShardID,0)==1)
					{
					#pragma omp critical (addhist)
					{
					HistList[i][GrpLen[i]]=HistID;
					ShardList[i][GrpLen[i]]=ShardID;
					GrpLen[i]++;
					}
					break;//do not look for other grp association
					}
				}
			}
		}
	}
	}
}
NList=0;
for(i=0;i<NSAMPLE;i++)
{
	NList+=GrpLen[i];
	printf("grpid="HBTIFMT",mass="HBTIFMT",GrpLen="HBTIFMT",Max="HBTIFMT"\n",grpid[i],halosize[MaxSnap-1][grpid[i]].Mvir[0],GrpLen[i],NumHist/10);fflush(stdout);
}
printf("NList="HBTIFMT",PartMass=%g\n",NList,PartMass);
Mhost=mymalloc(sizeof(HBTReal)*NList);
Msub=mymalloc(sizeof(HBTReal)*NList);
tinfall=mymalloc(sizeof(HBTReal)*NList);
Mmax=mymalloc(sizeof(HBTReal)*NList);
tmax=mymalloc(sizeof(HBTReal)*NList);
MRvir=mymalloc(sizeof(HBTReal)*NList);
tRvir=mymalloc(sizeof(HBTReal)*NList);
DirectInfall=mymalloc(sizeof(HBTInt)*NList);

HBTInt j,k=0;
for(i=0;i<NSAMPLE;i++)
{
	for(j=0;j<GrpLen[i];j++)
	{	
	HistID=HistList[i][j];
	ShardID=ShardList[i][j];
	SnapTidal=HistoryRevShard.Par[HistID][ShardID].SnapTidal;
	SnapRvir=HistoryRevShard.Par[HistID][ShardID].SnapRvir;
	SnapInfall=get_snapinfall(HistID,ShardID);//what about change this to SnapRvir
	Node=GetMember(&EvoCatPre,HistID,SnapInfall);
	if(Node==NULL)
	{
		Mhost[k]=0.;
		Msub[k]=0.;
		tinfall[k]=0.;
		Mmax[k]=0.;
		tmax[k]=0.;
		MRvir[k]=0.;
		tRvir[k]=0.;
		DirectInfall[k]=0;
		printf("NULL node: "HBTIFMT","HBTIFMT"\n",HistID,ShardID);
	}
	else
	{
		if(Node->SubRank>0)
		Mhost[k]=Node->Mdm*PartMass;
		else
		Mhost[k]=Node->Mhost*PartMass;
		Msub[k]=Node->Mdm*PartMass;
		tinfall[k]=scaleF[SnapInfall];
		Mmax[k]=GetMember(&EvoCatPre,HistID,SnapTidal)->Mdm*PartMass;
		tmax[k]=scaleF[SnapTidal];
		MRvir[k]=GetMember(&EvoCatPre,HistID,SnapRvir)->Mdm*PartMass;
		tRvir[k]=scaleF[SnapTidal];
		DirectInfall[k]=is_direct_grp_family(mainid[i],HistID,ShardID);
	}
	k++;
	}
}

sprintf(buf,"%s/anal/infalldata",SUBCAT_DIR);
myfopen(fp,buf,"w");
for(i=0;i<NSAMPLE;i++)
{
HBTReal Mvir;
Mvir=halosize[MaxSnap-1][grpid[i]].Mvir[0]*PartMass;
fwrite(&Mvir,sizeof(HBTReal),1,fp);
}
fwrite(GrpLen,sizeof(HBTInt),NSAMPLE,fp);
fwrite(Mhost,sizeof(HBTReal),NList,fp);
fwrite(Msub,sizeof(HBTReal),NList,fp);
fwrite(tinfall,sizeof(HBTReal),NList,fp);
fwrite(Mmax,sizeof(HBTReal),NList,fp);
fwrite(tmax,sizeof(HBTReal),NList,fp);
fwrite(MRvir,sizeof(HBTReal),NList,fp);
fwrite(tRvir,sizeof(HBTReal),NList,fp);
fwrite(DirectInfall,sizeof(HBTInt),NList,fp);
fclose(fp);


HBTInt Hid[NSAMPLE];
free(Mhost);
free(Msub);
free(tinfall);
NList=0;
for(i=0;i<NSAMPLE;i++)
{
Hid[i]=Sub2Hist[MaxSnap-1][mainid[i]];
GrpLen[i]=EvoCat.HistLen[Hid[i]];
NList+=GrpLen[i];
}
Mhost=mymalloc(sizeof(HBTReal)*NList);
Msub=mymalloc(sizeof(HBTReal)*NList);
tinfall=mymalloc(sizeof(HBTReal)*NList);

HISTORY *Hist;
for(i=0,k=0;i<NSAMPLE;i++)
{	
	Hist=&EvoCatPre.History[Hid[i]];
	for(Nsnap=Hist->SnapBirth,j=0;Nsnap<Hist->SnapDeath;Nsnap++,j++,k++)
	{
		Mhost[k]=Hist->Member[j].Mhost*PartMass;
		Msub[k]=Hist->Member[j].Mdm*PartMass;
		tinfall[k]=scaleF[Nsnap];
	}
}
if(k!=NList)
{
	printf("error: k="HBTIFMT",Nlist="HBTIFMT"\n",k,NList);
	exit(1);
}
sprintf(buf,"%s/anal/mainsample",SUBCAT_DIR);
myfopen(fp,buf,"w");
for(i=0;i<NSAMPLE;i++)
{
HBTReal Mvir;
Mvir=halosize[MaxSnap-1][grpid[i]].Mvir[0]*PartMass;
fwrite(&Mvir,sizeof(HBTReal),1,fp);
}
fwrite(GrpLen,sizeof(HBTInt),NSAMPLE,fp);
fwrite(Mhost,sizeof(HBTReal),NList,fp);
fwrite(Msub,sizeof(HBTReal),NList,fp);
fwrite(tinfall,sizeof(HBTReal),NList,fp);
fclose(fp);

return 0;
}

