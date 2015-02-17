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

#define Mmin 0.1
#define Mmax 0.2
struct SubNodeExp
{
	int HostID;
	int subid;
	int Mdm;
	int Mgas;
	int Mhost;//host halo mass
	int Mcen;//host subhalo mass
	float CoM[3];//relative position,comoving
	float VCoM[3];//relative Vel,physical
	float Pot;
	float Kin;
	float AM[3];
	//~ float HPot;//host Pot
	//~ float HKin;//host Kin
	//~ float HAM[3];
};
void load_history_ext(int *P2NumHist,SubHist **P2ClusterHistory,char *subcatdir)
{
FILE *fp;
char buf[1024];
int NumHist;
SubHist *ClusterHistory;

  sprintf(buf, "%s/history/ClusterHistoryExt", subcatdir);
  myfopen(fp,buf,"r");
  fread(&NumHist,sizeof(int),1,fp);
  ClusterHistory=mymalloc(sizeof(SubHist)*NumHist);
  fread(ClusterHistory,sizeof(SubHist),NumHist,fp);
  fclose(fp);
  *P2NumHist=NumHist;
  *P2ClusterHistory=ClusterHistory;		
}
int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
//~ CATALOGUE Cat;
//~ GASSUBCAT GSubCat;
struct SubNodeExp (*History)[MaxSnap];

//~ char gasdir[1024]=GASCAT_DIR;
char buf[1024];FILE *fp;
int Nsnap,i,SnapInfall,SnapBirth;
float Mratio;

SubHist *ClusterHistory;
int NumHist,HistID,subid,hostsubid,* (Sub2Hist[MaxSnap]),Nsubs;
int *HistList,NList,*NNode;
SubNode *OldNode;

logfile=stdout;

load_history_ext(&NumHist,&ClusterHistory,SUBCAT_DIR);
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
HistList=mymalloc(sizeof(int)*NumHist);
NList=0;
for(HistID=0;HistID<NumHist;HistID++)
{
	SnapInfall=ClusterHistory[HistID].SnapEnter+1;
	if(SnapInfall<MaxSnap)
	{
	Mratio=(float)(ClusterHistory[HistID].Member[SnapInfall].Mdm)/(float)(ClusterHistory[HistID].Member[SnapInfall].Mhost);
	if(Mratio>Mmin&&Mratio<Mmax&&ClusterHistory[HistID].Member[SnapInfall].HostID>=0)//to exclude quasi
	{
		HistList[NList]=HistID;
		NList++;		
	}
	}
}
printf("NList=%d\n",NList);
printf("sizeof(SubNodeExp)=%d\n",sizeof(struct SubNodeExp));
History=mymalloc(sizeof(struct SubNodeExp)*MaxSnap*NList);
NNode=calloc(NList,sizeof(int));
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	for(i=0;i<NList;i++)
	{
		HistID=HistList[i];
		SnapInfall=ClusterHistory[HistID].SnapEnter+1;
		SnapBirth=ClusterHistory[HistID].SnapBirth;
		OldNode=ClusterHistory[HistID].Member+Nsnap;
		if((Nsnap>=SnapBirth)&&(OldNode->Mdm)&&(OldNode->HostID>=0))//also exclude those whose hosthalo haven't been bore
		{
			NNode[i]++;
			hostsubid=SubCat.GrpOffset_Sub[OldNode->HostID];
			subid=OldNode->SubID;
			History[i][Nsnap].HostID=OldNode->HostID;
			History[i][Nsnap].subid=subid;
			History[i][Nsnap].Mdm=OldNode->Mdm;
			History[i][Nsnap].Mgas=OldNode->Mgas;
			History[i][Nsnap].Mhost=OldNode->Mhost;
			History[i][Nsnap].Mcen=SubCat.SubLen[hostsubid];
			History[i][Nsnap].CoM[0]=SubCat.Property[subid].CoM[0]-SubCat.Property[hostsubid].CoM[0];
			History[i][Nsnap].CoM[1]=SubCat.Property[subid].CoM[1]-SubCat.Property[hostsubid].CoM[1];
			History[i][Nsnap].CoM[2]=SubCat.Property[subid].CoM[2]-SubCat.Property[hostsubid].CoM[2];
			History[i][Nsnap].VCoM[0]=SubCat.Property[subid].VCoM[0]-SubCat.Property[hostsubid].VCoM[0];
			History[i][Nsnap].VCoM[1]=SubCat.Property[subid].VCoM[1]-SubCat.Property[hostsubid].VCoM[1];
			History[i][Nsnap].VCoM[2]=SubCat.Property[subid].VCoM[2]-SubCat.Property[hostsubid].VCoM[2];
			History[i][Nsnap].AM[0]=SubCat.Property[subid].AM[0];
			History[i][Nsnap].AM[1]=SubCat.Property[subid].AM[1];
			History[i][Nsnap].AM[2]=SubCat.Property[subid].AM[2];
			//~ History[i][Nsnap].HAM[0]=SubCat.Property[hostsubid].AM[0];
			//~ History[i][Nsnap].HAM[1]=SubCat.Property[hostsubid].AM[1];
			//~ History[i][Nsnap].HAM[2]=SubCat.Property[hostsubid].AM[2];
			History[i][Nsnap].Pot=SubCat.Property[subid].Pot;
			//~ History[i][Nsnap].HPot=SubCat.Property[hostsubid].Pot;
			History[i][Nsnap].Kin=SubCat.Property[subid].Kin;
			//~ History[i][Nsnap].HKin=SubCat.Property[hostsubid].Kin;
		}
		else
		{
			History[i][Nsnap].Mdm=0;
		}
	}
}
sprintf(buf,"%s/anal/history_%02d_%02d.dat",SUBCAT_DIR,(int)(10*Mmin),(int)(10*Mmax));
myfopen(fp,buf,"w");
fwrite(&NList,sizeof(int),1,fp);
for(i=0;i<NList;i++)
{
	//subid_birth
	fwrite(&(ClusterHistory[HistList[i]].Member[ClusterHistory[HistList[i]].SnapBirth].SubID),sizeof(int),1,fp);
	fwrite(NNode+i,sizeof(int),1,fp);
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
		if(History[i][Nsnap].Mdm)
		{
		fwrite(&Nsnap,sizeof(int),1,fp);
		fwrite(History[i]+Nsnap,sizeof(struct SubNodeExp),1,fp);
		}
	}
}
fwrite(&NList,sizeof(int),1,fp);
fclose(fp);

return 0;
}
