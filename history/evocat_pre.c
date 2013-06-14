//to construct light-weighted history with only identities
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
 
SUBCATALOGUE SubCat;
//~ CATALOGUE Cat;
BRFCAT BrfCat[MaxSnap];
EVOLUTIONCAT_Pre EvoCat;
SubNodePre *AllNodes;
HBTInt  MaxHist;

HBTInt  get_total_subs();
void get_snapinfall(HBTInt NumHist,HISTORY_Pre *History);
HBTInt new_history(HBTInt subid,HBTInt SnapBirth,HBTInt HistID,HBTInt Offset,HBTInt ProHistID);

int main(int argc,char **argv)
{
HBTInt Nsnap,i;
HBTInt * sp2pro,Npro,Nspl;
HBTInt NumHist,HistID,*Sub2Hist;
HBTInt ProSubID,subid,SnapBirth,HistOffset;
size_t mem_loaded;
logfile=stdout;
mem_loaded=0;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	mem_loaded+=load_brfcat(Nsnap,BrfCat+Nsnap);
	printf(""HBTIFMT".",Nsnap);fflush(stdout);
}
printf("data loaded: %ld GB memory used\n", mem_loaded/(1024L*1024L*1024L));
EvoCat.NNode=get_total_subs();
MaxHist=EvoCat.NNode;
EvoCat.History=mymalloc(sizeof(HISTORY_Pre)*MaxHist);
EvoCat.HistLen=mymalloc(sizeof(HBTInt)*MaxHist);
EvoCat.HistOffset=mymalloc(sizeof(HBTInt)*MaxHist);
AllNodes=mymalloc(sizeof(SubNodePre)*EvoCat.NNode);
Sub2Hist=mymalloc(sizeof(HBTInt));
Npro=0;
NumHist=0;
HistOffset=0;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	printf("Nsnap="HBTIFMT"\n",Nsnap);
	load_sp2pro(Nsnap,&Npro,&Nspl,&sp2pro,SUBCAT_DIR);
	HistID=NumHist;
	for(i=0;i<BrfCat[Nsnap].Nsubs;i++)
	{ 
		ProSubID=BrfCat[Nsnap].SubHalo[i].ProID;
		if(ProSubID<0)//birth 
		{
		  HistOffset+=new_history(i,Nsnap,HistID,HistOffset,-1);
		  HistID++;
	  	}
		else if(ProSubID>=Npro)//splitter
		{
		  HistOffset+=new_history(i,Nsnap,HistID,HistOffset,Sub2Hist[sp2pro[ProSubID]]);
		  HistID++;
		}		  
	}
	NumHist=HistID;
	free(Sub2Hist);
	Sub2Hist=mymalloc(sizeof(HBTInt)*BrfCat[Nsnap].Nsubs);
	for(HistID=0;HistID<NumHist;HistID++)
	{
		if(Nsnap>=EvoCat.History[HistID].SnapDeath) continue;
		SnapBirth=EvoCat.History[HistID].SnapBirth;
		subid=EvoCat.History[HistID].Member[Nsnap-SnapBirth].SubID;
		if(subid>=0)
		Sub2Hist[subid]=HistID;
		else
		{
			printf("error: subid<0 before death");
			exit(1);
		}
	}
	save_sub2hist(Nsnap,Sub2Hist,BrfCat[Nsnap].Nsubs,SUBCAT_DIR);
	free_sp2pro(sp2pro,Npro,Nspl);
}
if(HistOffset!=EvoCat.NNode)
{printf("error: Offset!=NNode,"HBTIFMT","HBTIFMT"\n",HistOffset,EvoCat.NNode);
exit(1);
}
EvoCat.HistLen=realloc(EvoCat.HistLen,sizeof(HBTInt)*NumHist);
EvoCat.HistOffset=realloc(EvoCat.HistOffset,sizeof(HBTInt)*NumHist);
EvoCat.History=realloc(EvoCat.History,sizeof(HISTORY_Pre)*NumHist);
EvoCat.NHist=NumHist;
get_snapinfall(NumHist,EvoCat.History);
save_evocat_pre(&EvoCat,SUBCAT_DIR);
free_evocat_pre(&EvoCat);
free(Sub2Hist);

return 0;
}

HBTInt new_history(HBTInt subid,HBTInt SnapBirth,HBTInt HistID,HBTInt Offset,HBTInt ProHistID)
{
	HISTORY_Pre *History;
	SubNodePre *NodeNow;
	HBTInt Nsnap;
	
	if(HistID>=MaxHist)
	{
		MaxHist*=2;
		EvoCat.HistLen=realloc(EvoCat.HistLen,sizeof(HBTInt)*MaxHist);
		EvoCat.HistOffset=realloc(EvoCat.HistOffset,sizeof(HBTInt)*MaxHist);
		EvoCat.History=realloc(EvoCat.History,sizeof(HISTORY_Pre)*MaxHist);
	}
	EvoCat.HistOffset[HistID]=Offset;
	History=&EvoCat.History[HistID];
	History->SnapBirth=SnapBirth;
	History->Member=AllNodes+Offset;
	History->ProHistID=ProHistID;
	for(Nsnap=SnapBirth,NodeNow=AllNodes+Offset;Nsnap<MaxSnap;Nsnap++,NodeNow++)
	{
	 if(subid<0) break;	
	 NodeNow->SubID=subid;
	 NodeNow->Mdm=BrfCat[Nsnap].SubHalo[subid].Mdm;
	 NodeNow->HostID=BrfCat[Nsnap].SubHalo[subid].HostID;
	 if(NodeNow->HostID<0)//quasi
	 NodeNow->SubRank=0;
	 else
	 NodeNow->SubRank=BrfCat[Nsnap].SubHalo[subid].SubRank;
	 subid=BrfCat[Nsnap].SubHalo[subid].DesID;
	}
	History->SnapDeath=Nsnap;
	EvoCat.HistLen[HistID]=Nsnap-SnapBirth;
	return EvoCat.HistLen[HistID];
}
void get_snapinfall(HBTInt NumHist,HISTORY_Pre *History)
{
	HBTInt i,NSnap;
	SubNodePre *Members;
	for(i=0;i<NumHist;i++)
	{
		Members=History[i].Member-History[i].SnapBirth;
		for(NSnap=History[i].SnapBirth;NSnap<History[i].SnapDeath;NSnap++)
		{
			if(Members[NSnap].SubRank>0)
			break;
		}
		History[i].SnapEnter=NSnap-1;//the snap just before crossing
	}
}
HBTInt get_total_subs()
{
	/* get total number of subhalos from all snapshots */
	
	FILE *fd;
	char buf[1024];
	HBTInt Nsnap,Ngroups,Nsubs;
	HBTInt Nall=0;
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
		sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, (int)Nsnap);
  		myfopen(fd,buf,"r");
  		fread(&Ngroups, sizeof(HBTInt), 1, fd);
  		fread(&Nsubs, sizeof(HBTInt), 1, fd);
		fclose(fd);
		//~ printf("Nsnap="HBTIFMT",Ngroups="HBTIFMT",Nsubs="HBTIFMT"\n",Nsnap,Ngroups,Nsubs);
		Nall+=Nsubs;
	}
	//~ printf("Nall="HBTIFMT",Mem=%2.1gG\n",Nall,Nall*4.0/1024/1024/1024);

	return Nall;
}
