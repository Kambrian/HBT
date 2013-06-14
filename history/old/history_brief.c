//to construct light-weighted history with only identities
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

typedef struct 
{
	int Mdm;//sub DM mass
	int SubID;
	int HostID;//host haloid
	int SubRank;
} SubNode;

typedef struct
{
	int SnapBirth;
	int SnapEnter;
	int SnapDeath;//when subhalo becomes under resolution
	int ProHistID;// -1 means the same as itself; positive integer means this is a splitter, 
					//so its progenitor is in a different history,given by ProHistID
	SubNode Member[MaxSnap];
 }SubHist;
 
SUBCATALOGUE SubCat;
CATALOGUE Cat;

void get_snapinfall(int NumHist,SubHist *ClusterHistory);
void add_member(SubNode *Member,int subid);
void add_Hist(SubHist *History,int SnapBirth,int subid);
void save_history_brf(int NumHist,SubHist *ClusterHistory,char *subcatdir);
void load_history_brf(int *P2NumHist,SubHist **P2ClusterHistory,char *subcatdir);
void save_sub2hist(int Nsnap,int *Sub2Hist,int Nsubs,char *subcatdir);
void load_sub2hist(int Nsnap,int **P2sub2hist,int *P2Nsubs,char *subcatdir);

int main(int argc,char **argv)
{
char buf[1024];
int SnapshotNum,i,SnapRange[2];

int *pro2dest,*sp2pro;
SubHist *ClusterHistory;
int NumHist,HistID,ProSubID,subid,Npro,Nspl,*Sub2Hist;

logfile=stdout;
ClusterHistory=mymalloc(sizeof(SubHist));
Sub2Hist=mymalloc(sizeof(int));
Npro=0;
NumHist=0;
for(SnapshotNum=0;SnapshotNum<MaxSnap;SnapshotNum++)
{
	load_sub_table(SnapshotNum,&SubCat,SUBCAT_DIR);
	load_pro2dest(SnapshotNum-1,&pro2dest,&Npro,SUBCAT_DIR);
	load_sp2pro(SnapshotNum,&Npro,&Nspl,&sp2pro,SUBCAT_DIR);
	ClusterHistory=realloc(ClusterHistory,sizeof(SubHist)*(SubCat.Nbirth+SubCat.Nsplitter+NumHist));
										//this would be big enough since NumHistNew=NumHist+Nbirth+Nsplitter-Nsplitter_death
										//those death splitters dies before they occupy a position as new sub, so there's no history for them
	for(HistID=0;HistID<NumHist;HistID++)//existing history
	{
		ProSubID=ClusterHistory[HistID].Member[SnapshotNum-1].SubID;//ok even if prosubid=-1;
		add_member(&(ClusterHistory[HistID].Member[SnapshotNum]),pro2dest[ProSubID]);
	}
	for(i=0;i<SubCat.Nsubs;i++)//new history
	{
		ProSubID=SubCat.HaloChains[i].ProSubID;
		if(ProSubID<0)//birth
		{
			ClusterHistory[HistID].ProHistID=-1;
			add_Hist(ClusterHistory+HistID,SnapshotNum,i);
			HistID++;
		}
		else if(ProSubID>=Npro)//splitter,and alive
		{
			ProSubID=sp2pro[ProSubID];//get mother id
			ClusterHistory[HistID].ProHistID=Sub2Hist[ProSubID];
			add_Hist(ClusterHistory+HistID,SnapshotNum,i);
			memcpy(ClusterHistory[HistID].Member,ClusterHistory[Sub2Hist[ProSubID]].Member,sizeof(SubNode)*SnapshotNum);
			HistID++;
		}
	}	
	NumHist=HistID;
	
	free(Sub2Hist);
	Sub2Hist=mymalloc(sizeof(int)*SubCat.Nsubs);
	for(HistID=0;HistID<NumHist;HistID++)
	{
		subid=ClusterHistory[HistID].Member[SnapshotNum].SubID;
		if(subid>=0)
		Sub2Hist[subid]=HistID;
	}
	save_sub2hist(SnapshotNum,Sub2Hist,SubCat.Nsubs,SUBCAT_DIR);
	
	free_sp2pro(sp2pro,Npro,Nspl);
	free_pro2dest(pro2dest);
	Npro=SubCat.Nsubs;	
	free_sub_table(&SubCat);
}
get_snapinfall(NumHist,ClusterHistory);
int NNodes=0,*HistLen,*HistOffset;
HistOffset=mymalloc(sizeof(int)*NumHist);
HistLen=mymalloc(sizeof(int)*NumHist);
for(HistID=0;HistID<NumHist;HistID++)
{
	int Nsnap;
	for(Nsnap=ClusterHistory[HistID].SnapBirth;Nsnap<MaxSnap;Nsnap++)
	if(ClusterHistory[HistID].Member[Nsnap].SubID<0) break;
	ClusterHistory[HistID].SnapDeath=Nsnap;
	HistOffset[HistID]=NNodes;
	HistLen[HistID]=ClusterHistory[HistID].SnapDeath-ClusterHistory[HistID].SnapBirth;
	NNodes+=HistLen[HistID];
}
printf("NNodes=%d\n",NNodes);
save_history_brf(NumHist,ClusterHistory,SUBCAT_DIR);
free(ClusterHistory);
free(Sub2Hist);

return 0;
}

void get_snapinfall(int NumHist,SubHist *ClusterHistory)
{
	int i,NSnap;
	for(i=0;i<NumHist;i++)
	{
		for(NSnap=ClusterHistory[i].SnapBirth;NSnap<MaxSnap;NSnap++)
		{
			if(ClusterHistory[i].Member[NSnap].SubRank>0)
			break;
		}
		ClusterHistory[i].SnapEnter=NSnap-1;//the snap just before crossing
	}
}
void add_member(SubNode *Member,int subid)
{
	if(subid<0)
	{
	Member->Mdm=0;
	Member->SubID=-1;
	Member->HostID=-2;
	Member->SubRank=-1;
	}
	else
	{
	Member->Mdm=SubCat.SubLen[subid];
	Member->SubID=subid;
	Member->HostID=SubCat.HaloChains[subid].HostID;
	if(Member->HostID<0)//quasi halo
	Member->SubRank=0;
	else
	Member->SubRank=SubCat.SubRank[subid];
	}
}
void add_Hist(SubHist *History,int SnapBirth,int subid)
{
	int i;
	History->SnapBirth=SnapBirth;
	for(i=0;i<SnapBirth;i++)
		add_member(History->Member+i,-1);
	add_member(History->Member+SnapBirth,subid);
}

void save_history_brf(int NumHist,SubHist *ClusterHistory,char *subcatdir)
{
FILE *fp;
char buf[1024];

  sprintf(buf, "%s/history2/ClusterHistory_brf", subcatdir);
  myfopen(fp,buf,"w");
  fwrite(&NumHist,sizeof(int),1,fp);
  fwrite(ClusterHistory,sizeof(SubHist),NumHist,fp);
  fclose(fp);		
}
void load_history_brf(int *P2NumHist,SubHist **P2ClusterHistory,char *subcatdir)
{
FILE *fp;
char buf[1024];
int NumHist;
SubHist *ClusterHistory;

  sprintf(buf, "%s/history/ClusterHistory_brf", subcatdir);
  myfopen(fp,buf,"r");
  fread(&NumHist,sizeof(int),1,fp);
  ClusterHistory=mymalloc(sizeof(SubHist)*NumHist);
  fread(ClusterHistory,sizeof(SubHist),NumHist,fp);
  fclose(fp);
  *P2NumHist=NumHist;
  *P2ClusterHistory=ClusterHistory;		
}

void save_sub2hist(int Nsnap,int *Sub2Hist,int Nsubs,char *subcatdir)
{
FILE *fp;
char buf[1024];

  sprintf(buf, "%s/history2/sub2hist_%03d", subcatdir, Nsnap);
  myfopen(fp,buf,"w");
  fwrite(&Nsubs,sizeof(int),1,fp);
  fwrite(Sub2Hist,sizeof(int),Nsubs,fp);
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
