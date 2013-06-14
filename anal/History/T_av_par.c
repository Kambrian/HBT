/*extract mass-averaged internal energy for infalling sat
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



int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
//~ GASHALOCAT GCat;
GASSUBCAT GSubCat;
struct SubNodeExp (**History);

char buf[1024];FILE *fp;
int Nsnap,i,SnapInfall,SnapTidal,SnapBirth,SnapEnd;
float Mmin,Mmax,Mratio;//mass ratio at RTidal

int NumHist,NumShards,HistID,ShardID,subid,hostsubid,Nsubs;
int *HistList,*ShardList,NList,*NNode;
SubNode *OldNode;


logfile=stdout;


sprintf(buf,"%s/anal/history_%03d_%03d.dat",SUBCAT_DIR,(int)(100*Mmin),(int)(100*Mmax));
myfopen(fp,buf,"r");
fread(&NList,sizeof(int),1,fp);
NNode=mymalloc(sizeof(int)*NList);
History=mymalloc(sizeof(struct SubNodeExp *)*NList);
HistSnap=mymalloc(sizeof(int *)*NList);
for(i=0;i<NList;i++)
{
	fread(NNode+i,sizeof(int),1,fp);
	History[i]=mymalloc(sizeof(struct SubNodeExp)*NNode[i])
	HistSnap[i]=mymalloc(sizeof(int)*NNode[i]);
	for(j=0;j<NNode[i];j++)
	{
		fread(HistSnap[i]+j,sizeof(int),1,fp);
		fread(History[i]+j,sizeof(struct SubNodeExp),1,fp);
	}
}
fread(&NList2,sizeof(int),1,fp);
if(NList2!=NList)
{
	fprintf(logfile,"error reading file %s: Nlist= %d, %d\n",buf,NList,NList2);
	exit(1);
}
fclose(fp);

sprintf(buf,"%s/anal/historypar_%03d_%03d.dat",SUBCAT_DIR,(int)(100*Mmin),(int)(100*Mmax));
myfopen(fp,buf,"r");
fread(&NList,sizeof(int),1,fp);
if(NList2!=NList)
{
	fprintf(logfile,"error reading file header %s: Nlist= %d, %d\n",buf,NList,NList2);
	exit(1);
}
HistPar=mymalloc(sizeof(struct ShardParam)*NList);
for(i=0;i<NList;i++)
{
	fread(HistPar+i,sizeof(struct ShardParam),1,fp);
}
fread(&NList2,sizeof(int),1,fp);
if(NList2!=NList)
{
	fprintf(logfile,"error reading file %s: Nlist= %d, %d\n",buf,NList,NList2);
	exit(1);
}
fclose(fp);

for(i=0;i<NList;i++)
{
	SnapInf=HistPar[i].SnapRvir;
	load_ga
}

return 0;
}


