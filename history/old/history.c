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

typedef struct
{
	int mass;//fof mass
	int Mvir[3];
	float Rvir[3];//[tophat,c200,b200],comoving
	int flag_badvir[3];
	int flag_fakehalo;//set to 1 when halo is not self-bound,to 0 otherwise.
} HALOSIZE;
void load_halo_size(HALOSIZE *halosize,int Ngroups,int Nsnap);
int load_halo_concentration(float *halocon,int Nsnap);

SUBCATALOGUE SubCat;
CATALOGUE Cat;
GASSUBCAT GSubCat;
HALOSIZE *halosize;
float *halocon;

void get_snapinfall(int NumHist,SubHist *ClusterHistory);
void add_member(SubNode *Member,int subid);
void add_Hist(SubHist *History,int SnapBirth,int subid);

int main(int argc,char **argv)
{
char gasdir[1024]=GASCAT_DIR;
char buf[1024],filemode[4];
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
	load_sub_catalogue(SnapshotNum,&SubCat,SUBCAT_DIR);
	load_gassubcat(SnapshotNum,&GSubCat,GASCAT_DIR);
	//~ load_group_catalogue(SnapshotNum,&Cat,GRPCAT_DIR);
	halosize=mymalloc(sizeof(HALOSIZE)*SubCat.Ngroups);
	load_halo_size(halosize,SubCat.Ngroups,SnapshotNum);
	halocon=mymalloc(sizeof(float)*SubCat.Ngroups);
	load_halo_concentration(halocon,SnapshotNum);
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
	//~ free_catalogue(&Cat);
	free(halocon);
	free(halosize);
	erase_gassubcat(&GSubCat);
	Npro=SubCat.Nsubs;	
	erase_sub_catalogue(&SubCat);
}
get_snapinfall(NumHist,ClusterHistory);
save_history(NumHist,ClusterHistory,SUBCAT_DIR);
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
	Member->Mgas=0;
	Member->Mhost=0;
	Member->Chost=0;
	Member->SubID=-1;
	Member->HostID=-2;
	Member->SubRank=-1;
	}
	else
	{
	Member->Mdm=SubCat.SubLen[subid];
	Member->Mgas=GSubCat.SubLen[subid];
	Member->SubID=subid;
	Member->HostID=SubCat.HaloChains[subid].HostID;
	if(Member->HostID<0)//quasi halo
	{
	Member->SubRank=0;
	Member->Mhost=0;
	Member->Chost=0;
	}
	else
	{
	Member->SubRank=SubCat.SubRank[subid];
	//~ Member->Mhost=Cat.Len[Member->HostID];
	Member->Mhost=halosize[Member->HostID].Mvir[0];
	Member->Chost=halocon[Member->HostID];
	}
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
void load_halo_param(float halopars[][7],int halostatus[],int Nsnap,int Ngroups)
{
	int grpid,Ngroups2;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_param_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	fread(&Ngroups,sizeof(int),1,fp);
	for(grpid=0;grpid<Ngroups;grpid++)
	{
		fread(halostatus+grpid,sizeof(int),1,fp);
		fread(halopars[grpid],sizeof(float),7,fp);
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
