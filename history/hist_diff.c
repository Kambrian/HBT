//test file; ignore this.
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
//~ #include "history_vars.h"
//~ #include "history_proto.h"

typedef struct 
{
	int Mdm;//sub DM mass
	int Mhost;//virial mass if possible,otherwise fof-mass
	float Chost;//-1 if not NFW-fittable
	int SubID;
	int HostID;//host haloid,now dynamically extended in EvoCatRev
	int SubRank;
} SubNode;
typedef struct
{
	int SnapBirth;
	int SnapDeath;//when subhalo becomes under resolution
	int SnapEnter;
	int ProHistID;// -1 means the same as itself; positive integer means this is a splitter, 
					//so its progenitor is in a different history,given by ProHistID
	SubNode *Member;
} HISTORY;
typedef struct
{
	float PartMass;
	int NNode;
	int NHist;
	int *HistLen;
	int *HistOffset;
	HISTORY *History;
}EVOLUTIONCAT;
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

HALOSIZE *(halosize[MaxSnap]);
float *(halocon[MaxSnap]);
float PartMass;
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
void load_historyshards_old(struct HistoryShards *HistoryShard)
{
int NumHist,i,Offset;
char buf[1024];
FILE *fp;
  	
  sprintf(buf, "%s/history_old/HistoryShards", SUBCAT_DIR);
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
void load_sub2hist(int Nsnap,int **P2sub2hist,int *P2Nsubs,char *subcatdir)
{
FILE *fp;
char buf[1024];
int Nsubs,*Sub2Hist;

  sprintf(buf, "%s/history/sub2hist_%03d", subcatdir, Nsnap);
  myfopen(fp,buf,"r");
  fread(&Nsubs,sizeof(int),1,fp);
  Sub2Hist=mymalloc(sizeof(int)*Nsubs);
  fread(Sub2Hist,sizeof(int),Nsubs,fp);
  fclose(fp);
  *P2Nsubs=Nsubs;
  *P2sub2hist=Sub2Hist;
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
void load_evocat_rev(EVOLUTIONCAT *EvoCat,char *subcatdir)
{
FILE *fp;
char buf[1024];
int i;

  sprintf(buf, "%s/history/EvoCat_rev", subcatdir);
  myfopen(fp,buf,"r");
  fread(&EvoCat->PartMass,sizeof(float),1,fp);
  fread(&EvoCat->NNode,sizeof(int),1,fp);
  fread(&EvoCat->NHist,sizeof(int),1,fp);
  EvoCat->HistLen=mymalloc(sizeof(int)*EvoCat->NHist);
  EvoCat->HistOffset=mymalloc(sizeof(int)*EvoCat->NHist);
  EvoCat->History=mymalloc(sizeof(HISTORY)*EvoCat->NHist);
  EvoCat->History[0].Member=mymalloc(sizeof(SubNode)*EvoCat->NNode);
  fread(EvoCat->HistLen,sizeof(int),EvoCat->NHist,fp);
  fread(EvoCat->HistOffset,sizeof(int),EvoCat->NHist,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapBirth,sizeof(int),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapDeath,sizeof(int),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapEnter,sizeof(int),1,fp); 
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].ProHistID,sizeof(int),1,fp); 	 
 	
 fread(EvoCat->History[0].Member,sizeof(SubNode),EvoCat->NNode,fp);	 
 for(i=0;i<EvoCat->NHist;i++)
 	EvoCat->History[i].Member=EvoCat->History[0].Member+EvoCat->HistOffset[i];
  fclose(fp);
}
SubNode *GetMember(EVOLUTIONCAT *EvoCat,int HistID,int Nsnap)
{
	int SnapBirth,SnapDeath;
	SnapBirth=EvoCat->History[HistID].SnapBirth;
	SnapDeath=EvoCat->History[HistID].SnapDeath;
	if(Nsnap<SnapBirth||Nsnap>=SnapDeath)
		return NULL;
	else
		return &(EvoCat->History[HistID].Member[Nsnap-SnapBirth]);
}
void load_history_brf(EVOLUTIONCAT *EvoCat,char *subcatdir)
{
FILE *fp;
char buf[1024];
int i;

  sprintf(buf, "%s/history/EvoCat_brf", subcatdir);
  myfopen(fp,buf,"r");
  fread(&EvoCat->NNode,sizeof(int),1,fp);
  fread(&EvoCat->NHist,sizeof(int),1,fp);
  EvoCat->HistLen=mymalloc(sizeof(int)*EvoCat->NHist);
  EvoCat->HistOffset=mymalloc(sizeof(int)*EvoCat->NHist);
  EvoCat->History=mymalloc(sizeof(HISTORY)*EvoCat->NHist);
  EvoCat->History[0].Member=mymalloc(sizeof(SubNode)*EvoCat->NNode);
  fread(EvoCat->HistLen,sizeof(int),EvoCat->NHist,fp);
  fread(EvoCat->HistOffset,sizeof(int),EvoCat->NHist,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapBirth,sizeof(int),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapDeath,sizeof(int),1,fp);
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].SnapEnter,sizeof(int),1,fp); 
  for(i=0;i<EvoCat->NHist;i++)
 	 fread(&EvoCat->History[i].ProHistID,sizeof(int),1,fp); 	 
 for(i=0;i<EvoCat->NNode;i++)
 {	
 	fread(&EvoCat->History[0].Member[i].Mdm,sizeof(int),1,fp);	 
	fread(&EvoCat->History[0].Member[i].SubID,sizeof(int),1,fp);	 
	fread(&EvoCat->History[0].Member[i].HostID,sizeof(int),1,fp);	 
	fread(&EvoCat->History[0].Member[i].SubRank,sizeof(int),1,fp);	 
 }
 for(i=0;i<EvoCat->NHist;i++)
 	EvoCat->History[i].Member=EvoCat->History[0].Member+EvoCat->HistOffset[i];
  fclose(fp);
}
void load_evocat_raw(EVOLUTIONCAT *EvoCat)
{
	int HistID,Nsnap;
	SubNode *Member;
	load_history_brf(EvoCat,SUBCAT_DIR);
	EvoCat->PartMass=header.mass[1];;
	for(HistID=0;HistID<EvoCat->NHist;HistID++)
	{
		for(Nsnap=EvoCat->History[HistID].SnapBirth;Nsnap<EvoCat->History[HistID].SnapDeath;Nsnap++)
		{
			Member=EvoCat->History[HistID].Member-EvoCat->History[HistID].SnapBirth+Nsnap;
			int HostID;
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

int main()
{
SUBCATALOGUE SubCat;
EVOLUTIONCAT EvoCat,EvoCatRev;
struct HistoryShards HistoryRevShard,HistoryOldShard;
int NumHist,HistID,mainid,mainid2,Nsnap;
int * (Sub2Hist[MaxSnap]),Nsubs;
TYPE_XYZ *(SubCoM[MaxSnap]);
int Ngroups;
int HIDdiff[100],i,j,k;

for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
Ngroups=read_Ngroups(Nsnap);
halosize[Nsnap]=mymalloc(sizeof(HALOSIZE)*Ngroups);
load_halo_size(halosize[Nsnap],Ngroups,Nsnap);
halocon[Nsnap]=mymalloc(sizeof(float)*Ngroups);
load_halo_concentration(halocon[Nsnap],Nsnap);
load_particle_header(Nsnap,SNAPSHOT_DIR);
}

load_evocat_raw(&EvoCat);
load_evocat_rev(&EvoCatRev,SUBCAT_DIR);
load_historyshards(&HistoryRevShard);
load_historyshards_old(&HistoryOldShard);
j=0;
for(i=0;i<HistoryRevShard.NumHist;i++)
{
	if(HistoryRevShard.NBirth[i]!=HistoryOldShard.NBirth[i])
	{
	 	HIDdiff[j]=i;
		printf("%d,%d:%d,%d\n",j,i,HistoryRevShard.NBirth[i],HistoryOldShard.NBirth[i]);
		j++;
		if(j>=100) break;
	}
	for(k=0;k<HistoryRevShard.NBirth[i];k++)
	{
	if(HistoryRevShard.Par[i][k].SnapBirth!=HistoryOldShard.Par[i][k].SnapBirth)
	{
	 	HIDdiff[j]=i;
		printf("%d,%d,%d:birth%d,%d\n",j,i,k,HistoryRevShard.Par[i][k].SnapBirth,HistoryOldShard.Par[i][k].SnapBirth);
		j++;
		if(j>=100) break;
	}
	if(HistoryRevShard.Par[i][k].SnapRvir!=HistoryOldShard.Par[i][k].SnapRvir)
	{
	 	HIDdiff[j]=i;
		printf("%d,%d,%d:vir%d,%d\n",j,i,k,HistoryRevShard.Par[i][k].SnapRvir,HistoryOldShard.Par[i][k].SnapRvir);
		j++;
		if(j>=100) break;
	}
	if(HistoryRevShard.Par[i][k].SnapTidal!=HistoryOldShard.Par[i][k].SnapTidal)
	{
	 	HIDdiff[j]=i;
		printf("%d,%d,%d:tidal%d,%d\n",j,i,k,HistoryRevShard.Par[i][k].SnapTidal,HistoryOldShard.Par[i][k].SnapTidal);
		j++;
		if(j>=100) break;
	}
	if(HistoryRevShard.Par[i][k].Mrate[0]!=HistoryOldShard.Par[i][k].Mrate[0])
	{
	 	HIDdiff[j]=i;
		printf("%d,%d,%d:Mrate0 %f,%f\n",j,i,k,HistoryRevShard.Par[i][k].Mrate[0],HistoryOldShard.Par[i][k].Mrate[0]);
		j++;
		if(j>=100) break;
	}
	if(HistoryRevShard.Par[i][k].Kappa[1]!=HistoryOldShard.Par[i][k].Kappa[1])
	{
	 	HIDdiff[j]=i;
		printf("%d,%d,%d:Mrate0 %f,%f\n",j,i,k,HistoryRevShard.Par[i][k].Kappa[0],HistoryOldShard.Par[i][k].Kappa[0]);
		j++;
		if(j>=100) break;
	}
	}
}
mainid=read_mainsubid(20,605);
load_sub_catalogue(20,&SubCat,SUBCAT_DIR);
mainid2=SubCat.GrpOffset_Sub[605];
read_subpos(20,&SubCoM[20]);

return 0;
}
