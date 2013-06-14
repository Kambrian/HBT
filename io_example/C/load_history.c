/*to pick out the history for a specified subhalo*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#include "history_vars.h"
#include "history_proto.h"


struct HistoryShards HistoryRevShard;
EVOLUTIONCAT EvoCat;
int * (Sub2Hist[MaxSnap]);

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;

char buf[1024];
FILE *fp;
int Nsnap,i,SnapBirth,SnapEnd,SnapDeath;
int SnapLoad,subid;

int NumHist,NumShards,HistID,ShardID,Nsubs;

if(argc!=3)
{
printf("usage:%s [Nsnap] [subid]\n",argv[0]);
exit(1);
}
SnapLoad=atoi(argv[1]);
subid=atoi(argv[2]);
logfile=stdout;

load_evocat_rev(&EvoCat,SUBCAT_DIR);
NumHist=EvoCat.NHist;

load_historyshards(&HistoryRevShard);
NumShards=HistoryRevShard.NumShards;

printf("History data loading complete.\n"
       "Total %d subhalos divided into %d histories, "
	   "fragmenting into %d pieces\n", 
	   EvoCat.NNode, NumHist, NumShards);
fflush(stdout);
	   
/* convert subhalo id to history id */
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
load_sub2hist(Nsnap,Sub2Hist+Nsnap,&Nsubs,SUBCAT_DIR);
HistID=Sub2Hist[SnapLoad][subid];

SnapBirth=EvoCat.History[HistID].SnapBirth;
SnapDeath=EvoCat.History[HistID].SnapDeath;

printf("Subhalo %d at Snapshot %d belongs to history %d\n",
        subid, SnapLoad, HistID);

printf("It lives %d snapshots, from snapshot %d to %d\n",
        EvoCat.HistLen[HistID],SnapBirth,SnapDeath-1);

SubNode *Member;		
Member=GetMember(&EvoCat,HistID,SnapLoad);
printf("Subhalo Mass %d, Hosthalo virial mass %d, Hosthalo Concentration %f\n",
		Member->Mdm, Member->Mhost, Member->Chost);
printf("SubID %d, Host HaloID %d\n", Member->SubID, Member->HostID);

printf("This subhalo undergoes %d mergers all its life\n", HistoryRevShard.NBirth[HistID]);
if(0==HistoryRevShard.NBirth[HistID]) return 0;
printf("the first merger crosses virial radius at snapshot %d\n",HistoryRevShard.Par[HistID][0].SnapRvir);
printf("with circularity %f\n",sqrt(HistoryRevShard.Par[HistID][0].j2[1]));

/*get shardID*/
for(ShardID=0;ShardID<HistoryRevShard.NBirth[HistID];ShardID++)
	if(SnapLoad<HistoryRevShard.Par[HistID][ShardID].SnapBirth)
		break;
ShardID--;		

printf("Snapshot %d is during the #%d merger event\n",SnapLoad,ShardID);	
printf("it crosses virial radius at snapshot %d\n",HistoryRevShard.Par[HistID][ShardID].SnapRvir);
printf("with circularity %f\n",sqrt(HistoryRevShard.Par[HistID][ShardID].j2[1]));

return 0;
}


