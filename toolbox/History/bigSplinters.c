//to extract some big splinters
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


int * (Sub2Hist[MaxSnap]);
EVOLUTIONCAT_Pre EvoCat;

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
SubNodePre *Member;
short *ProHistMask;
int *Mpri,*Msec,*SnapPro,*SubIDPro,MemberID,ProHistID;
int *SubIDPri,*SubIDSec;
int Nsp;
char buf[1024];FILE *fp;
int Nsnap,Nsubs,i,subid;
logfile=stdout;

load_evocat_pre(&EvoCat,SUBCAT_DIR);
ProHistMask=mymalloc(sizeof(short)*EvoCat.NHist);	
Mpri=mymalloc(sizeof(int)*EvoCat.NHist);//real main descendent mass
Msec=mymalloc(sizeof(int)*EvoCat.NHist);//recorded main descendent mass
SnapPro=mymalloc(sizeof(int)*EvoCat.NHist);
SubIDPro=mymalloc(sizeof(int)*EvoCat.NHist);
SubIDPri=mymalloc(sizeof(int)*EvoCat.NHist);
SubIDSec=mymalloc(sizeof(int)*EvoCat.NHist);
for(i=0;i<EvoCat.NHist;i++)
	ProHistMask[i]=0;
	
for(i=0;i<EvoCat.NHist;i++)
{
	ProHistID=EvoCat.History[i].ProHistID;
	if(ProHistID<0) continue;
	if(0==ProHistMask[ProHistID])
	{
		ProHistMask[ProHistID]=1;//this is a known splinter's progenitor
		MemberID=EvoCat.History[i].SnapBirth-EvoCat.History[ProHistID].SnapBirth;
		if(EvoCat.History[i].SnapBirth>=EvoCat.History[ProHistID].SnapDeath)
		{
			Mpri[ProHistID]=0;
			SubIDPri[ProHistID]=-1;
		}
		else
		{
			Mpri[ProHistID]=EvoCat.History[ProHistID].Member[MemberID].Mdm;
			SubIDPri[ProHistID]=EvoCat.History[ProHistID].Member[MemberID].SubID;
		}
		Msec[ProHistID]=0;
		SnapPro[ProHistID]=EvoCat.History[i].SnapBirth-1;
		SubIDPro[ProHistID]=EvoCat.History[ProHistID].Member[MemberID-1].SubID;
		SubIDSec[ProHistID]=-1;
	}
	else if(SnapPro[ProHistID]!=EvoCat.History[i].SnapBirth-1) continue; //splinter at a different time, skip it
	if(EvoCat.History[i].Member[0].Mdm>Mpri[ProHistID]) 
	{
		Msec[ProHistID]=Mpri[ProHistID];
		SubIDSec[ProHistID]=SubIDPri[ProHistID];
		Mpri[ProHistID]=EvoCat.History[i].Member[0].Mdm;		
		SubIDPri[ProHistID]=EvoCat.History[i].Member[0].SubID;	
	}
	else if(EvoCat.History[i].Member[0].Mdm>Msec[ProHistID])
	{
		Msec[ProHistID]=EvoCat.History[i].Member[0].Mdm;
		SubIDSec[ProHistID]=EvoCat.History[i].Member[0].SubID;	
	} 
}

Nsp=0;
sprintf(buf,"%s/anal/sppair",SUBCAT_DIR);
myfopen(fp,buf,"w");
for(i=0;i<EvoCat.NHist;i++)
{
	if(ProHistMask[i]==0) continue;
	if(Msec[i]>100) 
	{
	Nsp++;
	fprintf(fp,"%d,%d,%d,%d,%d,%d\n",Mpri[i],Msec[i],SnapPro[i],SubIDPro[i],SubIDPri[i],SubIDSec[i]);
	}
}
printf("%d\n",Nsp);
fclose(fp);
//~ 
//~ qsort(Slist,Nsp,sizeof(struct snode),compare_Mmax);
//~ 
//~ sprintf(buf,"%s/anal/splist_1000",SUBCAT_DIR);
//~ 
//~ int Nsp0=Nsp;
//~ 
//~ if(Nsp>1000) Nsp=1000;
//~ myfopen(fp,buf,"w");
//~ for(i=0;i<Nsp;i++)
//~ fprintf(fp,"%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",Slist[i].SnapSp,Slist[i].Sid,Slist[i].Msp,
								//~ Slist[i].SnapMax,Slist[i].Mid,Slist[i].Mmax);
//~ fclose(fp);
//~ 
//~ Nsp=Nsp0;
//~ qsort(Slist,Nsp,sizeof(struct snode),compare_Msp);
//~ 
//~ sprintf(buf,"%s/anal/splist0_1000",SUBCAT_DIR);
//~ 
//~ if(Nsp>1000) Nsp=1000;
//~ myfopen(fp,buf,"w");
//~ for(i=0;i<Nsp;i++)
//~ fprintf(fp,"%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",Slist[i].SnapSp,Slist[i].Sid,Slist[i].Msp,
								//~ Slist[i].SnapMax,Slist[i].Mid,Slist[i].Mmax);
//~ fclose(fp);

return 0;
}

