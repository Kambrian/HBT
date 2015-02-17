//pick out failure to assign main descendent to splinters
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
struct snode
{
	int ProHistID;
	int SnapSp; //splinter's birth time
	int Sid;//splinter id at snapsp
	int Msp;
};
struct bnode
{
	int SnapSp;
	int IDdes;//recorded descendent id
	int IDmain;//biggest splinter (excluding IDdes)
	int IDpro;//progenitor
	int Mdes;
	int Mmain;
};
static int compare_sp(const void *a, const void *b)
{
	  if(((struct snode *) a)->ProHistID > ((struct snode *) b)->ProHistID )
		return +1;  //move right

	 if(((struct snode *) a)->ProHistID < ((struct snode *) b)->ProHistID )
		return -1;  //move left

	if(((struct snode *) a)->SnapSp > ((struct snode *) b)->SnapSp )
		return +1;
	
	if(((struct snode *) a)->SnapSp < ((struct snode *) b)->SnapSp )
		return -1;	
	
	if(((struct snode *) a)->Msp < ((struct snode *) b)->Msp )
		return +1;	
		
	if(((struct snode *) a)->Msp > ((struct snode *) b)->Msp )
		return -1;	
					
	  return 0;
}

static int compare_Mmain(const void *a, const void *b)//used to sort Mmax in descending order
{
	  if(((struct bnode *) a)->Mmain > ((struct bnode *) b)->Mmain )
		return -1;

	 if(((struct bnode *) a)->Mmain < ((struct bnode *) b)->Mmain )
		return +1;

	  return 0;
}
int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
short *ProHistMask;
struct snode *Slist;
struct bnode *blist;
int SnapSp,MemberID,ProHistID;
int Nwrong,Nsp;
char buf[1024];FILE *fp;
int Nsnap,Nsubs,i,subid;
logfile=stdout;

load_evocat_pre(&EvoCat,SUBCAT_DIR);
ProHistMask=mymalloc(sizeof(short)*EvoCat.NHist);	

Nsp=0;
for(i=0;i<EvoCat.NHist;i++)
	if(EvoCat.History[i].ProHistID>=0) Nsp++;

Slist=mymalloc(sizeof(struct snode)*Nsp);
Nsp=0;
for(i=0;i<EvoCat.NHist;i++)//collect splinters
{
	ProHistID=EvoCat.History[i].ProHistID;
	if(ProHistID<0) continue;
	Slist[Nsp].ProHistID=ProHistID;
	Slist[Nsp].SnapSp=EvoCat.History[i].SnapBirth;
	Slist[Nsp].Sid=EvoCat.History[i].Member[0].SubID;
	Slist[Nsp].Msp=EvoCat.History[i].Member[0].Mdm;
	Nsp++;
}
qsort(Slist,Nsp,sizeof(struct snode),compare_sp);//sort according to ProID, then Nsnap, then Msp
blist=mymalloc(sizeof(struct bnode)*Nsp);
int Nbr=0;//number of split events
ProHistID=-1;
SnapSp=-1;
for(i=0;i<Nsp;i++)
{
	if(ProHistID==Slist[i].ProHistID)
	{
		if(SnapSp==Slist[i].SnapSp) 
			continue; //skip other sp
		else
			SnapSp=Slist[i].SnapSp;
	}
	else
	{
		ProHistID=Slist[i].ProHistID;
		SnapSp=Slist[i].SnapSp;
	}
	blist[Nbr].SnapSp=SnapSp;
	HISTORY_Pre *Hist;
	Hist=&EvoCat.History[ProHistID];
	blist[Nbr].IDpro=Hist->Member[SnapSp-1-Hist->SnapBirth].SubID;
	if(SnapSp>=Hist->SnapDeath)
	{
		blist[Nbr].IDdes=-1;
		blist[Nbr].Mdes=0;
	}
	else
	{
	blist[Nbr].IDdes=Hist->Member[SnapSp-Hist->SnapBirth].SubID;
	blist[Nbr].Mdes=Hist->Member[SnapSp-Hist->SnapBirth].Mdm;	
	}
	blist[Nbr].IDmain=Slist[i].Sid;
	blist[Nbr].Mmain=Slist[i].Msp;
	Nbr++;		
}

qsort(blist,Nbr,sizeof(struct bnode),compare_Mmain);

Nsp=0;
Nwrong=0;
int N=0;
sprintf(buf,"%s/anal/spmass",SUBCAT_DIR);
myfopen(fp,buf,"w");
for(i=0;i<Nbr;i++)
{
	if(blist[i].Mmain>300||blist[i].Mdes>300) 
	{
	N++;
	fprintf(fp,"%d,%d,%d,%d,%d,%d\n",blist[i].Mdes,blist[i].Mmain,blist[i].SnapSp,blist[i].IDdes,blist[i].IDmain,blist[i].IDpro);
	if(blist[i].Mdes<blist[i].Mmain) Nwrong++;
	}
}
printf("%d,%d,%d,%f\n",Nbr,N,Nwrong,(float)Nwrong/N);
fclose(fp);

return 0;
}

