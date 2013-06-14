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

void save_history(int NumHist,SubHist *ClusterHistory,char *subcatdir)
{
FILE *fp;
char buf[1024];

  sprintf(buf, "%s/history/ClusterHistory", subcatdir);
  myfopen(fp,buf,"w");
  fwrite(&NumHist,sizeof(int),1,fp);
  fwrite(ClusterHistory,sizeof(SubHist),NumHist,fp);
  fclose(fp);		
}
void load_history(int *P2NumHist,SubHist **P2ClusterHistory,char *subcatdir)
{
FILE *fp;
char buf[1024];
int NumHist;
SubHist *ClusterHistory;

  sprintf(buf, "%s/history/ClusterHistory", subcatdir);
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

  sprintf(buf, "%s/history/sub2hist_%03d", subcatdir, Nsnap);
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
