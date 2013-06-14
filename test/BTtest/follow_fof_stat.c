//to test how much of a halo's particles can be marked as bound at least once throughout the subhalo's life
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "binding.c"

#define CoreFrac 0.1

int main(int nargc,char **argv)
{
int SnapLoad,SnapBind;
CATALOGUE Cat;
int GrpLen, *GrpPIDs,Lmain_removed,*PInd_main_removed;
int fofid,i,MaxSubLen,IniSubLen;
short *CatMask;
float s[3];
  FILE *fp;
  char buf[1024];

logfile=stdout;	
SnapLoad=atoi(argv[1]);
fofid=atoi(argv[2]);

		load_group_catalogue(SnapLoad,&Cat,GRPCAT_DIR);
		CatMask=mymalloc(sizeof(short)*NP_SIM);
		for(i=0;i<NP_SIM;i++)
		CatMask[i]=0;
		CatMask--;//PID from 1 to NP_SIM
		/*==back up the fof==*/
		GrpLen=Cat.Len[fofid];
		GrpPIDs=mymalloc(sizeof(int)*GrpLen);
		memcpy(GrpPIDs,Cat.PIDorIndex+Cat.Offset[fofid],sizeof(int)*GrpLen);
		#ifdef GRPINPUT_INDEX  //convert index to id
		load_particle_data(SnapLoad,SNAPSHOT_DIR);
		for(i=0;i<GrpLen;i++)
		#ifdef PID_ORDERED
			GrpPIDs[i]++;
		#else
			GrpPIDs[i]=Pdat.PID[GrpPIDs[i]];
		#endif
		#endif
		/*==compact Cat to contain just that fof==*/
		free_catalogue(&Cat);
		Cat.Nids=GrpLen;
		Cat.PIDorIndex=mymalloc(sizeof(int)*GrpLen);
		memcpy(Cat.PIDorIndex,GrpPIDs,sizeof(int)*GrpLen);
		for(SnapBind=SnapLoad;SnapBind<MaxSnap;SnapBind++)
		{
			load_particle_data(SnapBind,SNAPSHOT_DIR);
			fresh_ID2Index(Cat.PIDorIndex,Cat.Nids); 	
			unbind_core(&(Cat.Nids),&(Cat.PIDorIndex),s,&Lmain_removed,&PInd_main_removed,CoreFrac);		
			free(PInd_main_removed);
			if(SnapBind==SnapLoad) IniSubLen=Cat.Nids;
			for(i=0;i<Cat.Nids;i++)
			#ifdef PID_ORDERED
				CatMask[Cat.PIDorIndex[i]+1]=1;
			#else
				CatMask[Pdat.PID[Cat.PIDorIndex[i]]]=1;
			#endif
			Cat.Nids=GrpLen;
			free(Cat.PIDorIndex);
			Cat.PIDorIndex=mymalloc(sizeof(int)*GrpLen);
			memcpy(Cat.PIDorIndex,GrpPIDs,sizeof(int)*GrpLen);
		}
		
		MaxSubLen=0;
		for(i=1;i<=NP_SIM;i++)
			if(CatMask[i]) MaxSubLen++;
		sprintf(buf, "%s/anal/follow/follow_%03d_%d.stat",SUBCAT_DIR,SnapLoad,fofid);
		myfopen(fp,buf,"w");
		fprintf(fp,"MaxSubLen %d, IniSubLen %d, fofLen %d\n",MaxSubLen,IniSubLen,GrpLen);
		fclose(fp);
		return 0;
}
