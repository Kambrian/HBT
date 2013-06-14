//pick up maximum subhalo mass for quasi-halos' history
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

//~ #define SPLINTER
#define EJECTION

int * (Sub2Hist[MaxSnap]);
EVOLUTIONCAT_Pre EvoCat;
struct qnode
{
	int SnapQuasi;
	int Qid;  //id at snapquasi
	int Mquasi;
	int SnapMax;
	int Mid;//id at snapmax
	int Mmax;
};
void get_quasi_max(struct qnode *node)
{
	int HistID,Nsnap;
	SubNodePre *Member;
	HistID=Sub2Hist[node->SnapQuasi][node->Qid];
	node->Mmax=0;
	#ifdef SPLINTER
	if(EvoCat.History[HistID].ProHistID<0) return;
	#endif
	#ifdef EJECTION
	if(EvoCat.History[HistID].SnapEnter>=node->SnapQuasi) return;
	#endif
	for(Nsnap=node->SnapQuasi;Nsnap<EvoCat.History[HistID].SnapDeath;Nsnap++)
	{
		Member=&EvoCat.History[HistID].Member[Nsnap-EvoCat.History[HistID].SnapBirth];
		if(node->Mmax<Member->Mdm)
		{
			node->Mmax=Member->Mdm;
			node->SnapMax=Nsnap;
			node->Mid=Member->SubID;
		}
	}
}
static int compare_Mmax(const void *a, const void *b)//used to sort Mmax in descending order
{
	  if(((struct qnode *) a)->Mmax > ((struct qnode *) b)->Mmax )
		return -1;

	 if(((struct qnode *) a)->Mmax < ((struct qnode *) b)->Mmax )
		return +1;

	  return 0;
}
static int compare_Mquasi(const void *a, const void *b)//used to sort Mmax in descending order
{
	  if(((struct qnode *) a)->Mquasi > ((struct qnode *) b)->Mquasi )
		return -1;

	 if(((struct qnode *) a)->Mquasi < ((struct qnode *) b)->Mquasi )
		return +1;

	  return 0;
}
int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
char buf[1024];FILE *fp;
int Nsnap,Nsubs,Nquasi,i,subid;
struct qnode *Qlist;
logfile=stdout;

Nquasi=0;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
Nquasi+=SubCat.NQuasi;
free_sub_table(&SubCat);
printf("%d.",Nsnap);fflush(stdout);
}
printf("\n");
Qlist=mymalloc(sizeof(struct qnode)*Nquasi);
Nquasi=0;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	for(subid=SubCat.Nsubs-SubCat.NQuasi;subid<SubCat.Nsubs;subid++,Nquasi++)
	{
		Qlist[Nquasi].SnapQuasi=Nsnap;
		Qlist[Nquasi].Qid=subid;
		Qlist[Nquasi].Mquasi=SubCat.SubLen[subid];
	}
	free_sub_table(&SubCat);
}
load_evocat_pre(&EvoCat,SUBCAT_DIR);
for(i=0;i<Nquasi;i++)
get_quasi_max(Qlist+i);

//~ qsort(Qlist,Nquasi,sizeof(struct qnode),compare_Mmax);
qsort(Qlist,Nquasi,sizeof(struct qnode),compare_Mquasi);

#ifdef SPLINTER
sprintf(buf,"%s/anal/quasisplist_1000",SUBCAT_DIR);
#else 
#ifdef EJECTION
sprintf(buf,"%s/anal/quasi_ejectlist_1000",SUBCAT_DIR);
#else
sprintf(buf,"%s/anal/quasilist0_1000",SUBCAT_DIR);
#endif
#endif
if(Nquasi>1000) Nquasi=1000;
myfopen(fp,buf,"w");
for(i=0;i<Nquasi;i++)
fprintf(fp,"%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",Qlist[i].SnapQuasi,Qlist[i].Qid,Qlist[i].Mquasi,
								Qlist[i].SnapMax,Qlist[i].Mid,Qlist[i].Mmax);
fclose(fp);

return 0;
}

