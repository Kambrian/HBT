/* to show 2 selected snapshot for the selected subhalo history*/

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

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
//~ CATALOGUE Cat;
//~ GASSUBCAT GSubCat;
struct SubNodeExp (*History)[MaxSnap];

//~ char gasdir[1024]=GASCAT_DIR;
char buf[1024];FILE *fp;
int Nsnap,i,j,pid;

SubHist *ClusterHistory;
int HistID,subid[2],* (Sub2Hist[MaxSnap]),Nsubs,NumHist;
int subsnap[2]={59,63};
logfile=stdout;

load_history(&NumHist,&ClusterHistory,SUBCAT_DIR);
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);

HistID=Sub2Hist[6][497];
for (i=0;i<2;i++)
{
subid[i]=ClusterHistory[HistID].Member[subsnap[i]].SubID;
load_sub_catalogue(subsnap[i],&SubCat,SUBCAT_DIR);
load_particle_data(subsnap[i],SNAPSHOT_DIR);
fresh_ID2Index(&SubCat,FRSH_SUBCAT);
sprintf(buf,"%s/anal/subsnap_%02d_%d.dat",SUBCAT_DIR,subsnap[i],subid[i]);
myfopen(fp,buf,"w");
fwrite(SubCat.SubLen+subid[i],sizeof(int),1,fp);
for(j=0;j<SubCat.SubLen[subid[i]];j++)
{
	pid=SubCat.PSubArr[subid[i]][j];
	fwrite(Pdat.Pos[pid],sizeof(float),3,fp);
}
fwrite(SubCat.SubLen+subid[i],sizeof(int),1,fp);
fclose(fp);
erase_sub_catalogue(&SubCat);
}
for (i=1;i<2;i+=2)
{
subid[i]=ClusterHistory[HistID].Member[subsnap[i]].SubID;
load_sub_catalogue(subsnap[i],&SubCat,SUBCAT_DIR);
load_particle_data(subsnap[i-1],SNAPSHOT_DIR);
fresh_ID2Index(&SubCat,FRSH_SUBCAT);
sprintf(buf,"%s/anal/subsnap_%02d_%d.dat",SUBCAT_DIR,subsnap[i-1],subid[i]);
myfopen(fp,buf,"w");
fwrite(SubCat.SubLen+subid[i],sizeof(int),1,fp);
for(j=0;j<SubCat.SubLen[subid[i]];j++)
{
	pid=SubCat.PSubArr[subid[i]][j];
	fwrite(Pdat.Pos[pid],sizeof(float),3,fp);
}
fwrite(SubCat.SubLen+subid[i],sizeof(int),1,fp);
fclose(fp);
erase_sub_catalogue(&SubCat);
}
return 0;
}
