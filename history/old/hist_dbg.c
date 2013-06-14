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
int main()
{
SUBCATALOGUE SubCat;
SubHist *ClusterHistory,*ClusterHistoryRev;
struct HistoryShards HistoryRevShard;
int NumHist,HistID,mainid,mainid2,Nsnap;
int * (Sub2Hist[MaxSnap]),Nsubs;
TYPE_XYZ *(SubCoM[MaxSnap]);

for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);

load_history(&NumHist,&ClusterHistory,SUBCAT_DIR);
load_history_rev(&NumHist,&ClusterHistoryRev,SUBCAT_DIR);
load_historyshards(&HistoryRevShard);
mainid=read_mainsubid(20,605);
load_sub_catalogue(20,&SubCat,SUBCAT_DIR);
mainid2=SubCat.GrpOffset_Sub[605];
read_subpos(20,&SubCoM[20]);

return 0;
}
