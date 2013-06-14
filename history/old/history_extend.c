//extend history to associate hosts all along history
//this has been incorporated into and reploaced by history_revised.c
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

int read_mainsubid(int Nsnap,int hostid);

int main(int argc,char **argv)
{
//~ SUBCATALOGUE SubCat;
//~ CATALOGUE Cat;
//~ GASSUBCAT GSubCat;

//~ char gasdir[1024]=GASCAT_DIR;
char buf[1024];FILE *fp;
int Nsnap,i,SnapInfall,hostid,hosthistid,Nsubs;


SubHist *ClusterHistory,*ClusterHistoryExt;
int NumHist,HistID,* (Sub2Hist[MaxSnap]);

logfile=stdout;

load_history(&NumHist,&ClusterHistory,SUBCAT_DIR);
load_history(&NumHist,&ClusterHistoryExt,SUBCAT_DIR);
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);

for(HistID=0;HistID<NumHist;HistID++)
{
	SnapInfall=ClusterHistory[HistID].SnapEnter+1;
	if(SnapInfall<MaxSnap)
	{
	hostid=ClusterHistory[HistID].Member[SnapInfall].HostID;
	if(hostid<0)
	{
		printf("error: hostid<0\n");
		//exit(1);
	}
	else
	{
	hosthistid=Sub2Hist[SnapInfall][read_mainsubid(SnapInfall,hostid)];
	for(Nsnap=ClusterHistory[HistID].SnapBirth;Nsnap<SnapInfall;Nsnap++)
	{
		ClusterHistoryExt[HistID].Member[Nsnap].HostID=ClusterHistory[hosthistid].Member[Nsnap].HostID;
		ClusterHistoryExt[HistID].Member[Nsnap].Mhost=ClusterHistory[hosthistid].Member[Nsnap].Mhost;
	}
	for(Nsnap=SnapInfall+1;Nsnap<MaxSnap;Nsnap++)
	{
		if(ClusterHistory[HistID].Member[Nsnap].SubRank)
		{
			hostid=ClusterHistory[HistID].Member[Nsnap].HostID;
			hosthistid=Sub2Hist[Nsnap][read_mainsubid(Nsnap,hostid)];
		}
		else
		{
			ClusterHistoryExt[HistID].Member[Nsnap].HostID=ClusterHistory[hosthistid].Member[Nsnap].HostID;
			ClusterHistoryExt[HistID].Member[Nsnap].Mhost=ClusterHistory[hosthistid].Member[Nsnap].Mhost;
		}
	}
	}
	}
}

  sprintf(buf, "%s/history/ClusterHistoryExt", SUBCAT_DIR);
  myfopen(fp,buf,"w");
  fwrite(&NumHist,sizeof(int),1,fp);
  fwrite(ClusterHistoryExt,sizeof(SubHist),NumHist,fp);
  fclose(fp);		

return 0;
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
