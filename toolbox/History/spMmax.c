//pick up maximum subhalo mass for splinters' history
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
	int SnapSp;
	int Sid;  //id at snapquasi
	int Msp;
	int SnapMax;
	int Mid;//id at snapmax
	int Mmax;
};

static int compare_Mmax(const void *a, const void *b)//used to sort Mmax in descending order
{
	  if(((struct snode *) a)->Mmax > ((struct snode *) b)->Mmax )
		return -1;

	 if(((struct snode *) a)->Mmax < ((struct snode *) b)->Mmax )
		return +1;

	  return 0;
}

static int compare_Msp(const void *a, const void *b)//used to sort Mmax in descending order
{
	  if(((struct snode *) a)->Msp > ((struct snode *) b)->Msp )
		return -1;

	 if(((struct snode *) a)->Msp < ((struct snode *) b)->Msp )
		return +1;

	  return 0;
}
int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
SubNodePre *Member;
char buf[1024];FILE *fp;
int Nsnap,Nsubs,Nsp,i,subid;
struct snode *Slist;
logfile=stdout;

load_evocat_pre(&EvoCat,SUBCAT_DIR);
Nsp=0;
for(i=0;i<EvoCat.NHist;i++)
	if(EvoCat.History[i].ProHistID>=0) Nsp++;
	
Slist=mymalloc(sizeof(struct snode)*Nsp);	
Nsp=0;

for(i=0;i<EvoCat.NHist;i++)
{
	if(EvoCat.History[i].ProHistID<0) continue;
	
	Slist[Nsp].SnapSp=EvoCat.History[i].SnapBirth;
	Slist[Nsp].Sid=EvoCat.History[i].Member[0].SubID;
	Slist[Nsp].Msp=EvoCat.History[i].Member[0].Mdm;
	Slist[Nsp].Mmax=0;
	for(Nsnap=Slist[Nsp].SnapSp;Nsnap<EvoCat.History[i].SnapDeath;Nsnap++)
	{
		Member=&EvoCat.History[i].Member[Nsnap-EvoCat.History[i].SnapBirth];
		if(Slist[Nsp].Mmax<Member->Mdm)
		{
			Slist[Nsp].Mmax=Member->Mdm;
			Slist[Nsp].SnapMax=Nsnap;
			Slist[Nsp].Mid=Member->SubID;
		}
	}	
	Nsp++;
}

qsort(Slist,Nsp,sizeof(struct snode),compare_Mmax);

sprintf(buf,"%s/anal/splist_1000",SUBCAT_DIR);

int Nsp0=Nsp;

if(Nsp>1000) Nsp=1000;
myfopen(fp,buf,"w");
for(i=0;i<Nsp;i++)
fprintf(fp,"%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",Slist[i].SnapSp,Slist[i].Sid,Slist[i].Msp,
								Slist[i].SnapMax,Slist[i].Mid,Slist[i].Mmax);
fclose(fp);

Nsp=Nsp0;
qsort(Slist,Nsp,sizeof(struct snode),compare_Msp);

sprintf(buf,"%s/anal/splist0_1000",SUBCAT_DIR);

if(Nsp>1000) Nsp=1000;
myfopen(fp,buf,"w");
for(i=0;i<Nsp;i++)
fprintf(fp,"%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",Slist[i].SnapSp,Slist[i].Sid,Slist[i].Msp,
								Slist[i].SnapMax,Slist[i].Mid,Slist[i].Mmax);
fclose(fp);

return 0;
}

