#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

//~ #define SUBID_INPUT  //define this to use mainsubid as input
//choose one from below to choose a unbinding algorithm
//~ #define FOLLOW_MEAN
//~ #define FOLLOW_MINPOTH
//~ #define FOLLOW_POTW
//~ #define FOLLOW_CORE

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "binding.c"

int main(int nargc,char **argv)
{
int SnapLoad,SnapBind;
CATALOGUE Cat;
SUBCATALOGUE SubCat;
int GrpLen, *GrpPIDs,Lmain_removed,*PInd_main_removed;
int fofid,i,subid;
float s[3];
FILE *fp;
char buf[1024];
  
float CoreFrac=0.25;

logfile=stdout;	
SnapLoad=atoi(argv[1]);
#ifndef SUBID_INPUT
fofid=atoi(argv[2]);
#else
subid=atoi(argv[2]);
#endif
#ifdef FOLLOW_CORE
CoreFrac=atof(argv[3]);//0.04,0.1111,0.25,0.51; (5,3,2,1.4)
#endif
//~ SnapLoad=51;fofid=78;

		load_group_catalogue(SnapLoad,&Cat,GRPCAT_DIR);
		#ifdef SUBID_INPUT
		load_sub_table(SnapLoad,&SubCat,SUBCAT_DIR);
		fofid=SubCat.HaloChains[subid].HostID;
		printf("subid %d, fofid %d\n",subid,fofid);fflush(stdout);
		free_sub_table(&SubCat);
		#endif
		#ifdef FOLLOW_MEAN
		sprintf(buf, "%s/anal/follow/follow_%03d_%d.mean",SUBCAT_DIR,SnapLoad,fofid);
		#endif
		#ifdef FOLLOW_MINPOTH
		sprintf(buf, "%s/anal/follow/follow_%03d_%d.minpotH",SUBCAT_DIR,SnapLoad,fofid);
		#endif
		#ifdef FOLLOW_POTW
		sprintf(buf, "%s/anal/follow/follow_%03d_%d.potW",SUBCAT_DIR,SnapLoad,fofid);
		#endif
		#ifdef FOLLOW_CORE
		sprintf(buf, "%s/anal/follow/follow_%03d_%d.%d",SUBCAT_DIR,SnapLoad,fofid,(int)(100*CoreFrac));
		#endif
		#ifndef OUTPUT_CENTERS_TAG
	  	if(!(fp= fopen(buf, "w")))
		{
		fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
		exit(1);
		}
		#endif

		/*==back up the fof==*/
		GrpLen=Cat.Len[fofid];printf("%d\n",GrpLen);
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
		/*==make sub_infall to be that fof==*/
		SubCat.Ngroups=1;SubCat.GrpOffset_Sub=NULL;SubCat.GrpLen_Sub=NULL;
		SubCat.Nsubs=1;	create_sub_cat(&SubCat);
		SubCat.SubLen[0]=GrpLen;
		SubCat.PSubArr[0]=mymalloc(sizeof(int)*GrpLen);
		memcpy(SubCat.PSubArr[0],GrpPIDs,sizeof(int)*GrpLen);
		for(SnapBind=SnapLoad;SnapBind<MaxSnap;SnapBind++)
		{
			load_particle_data(SnapBind,SNAPSHOT_DIR);	
			fresh_ID2Index(Cat.PIDorIndex,Cat.Nids); 	
			fresh_ID2Index(SubCat.PSubArr[0],SubCat.SubLen[0]);
			#ifdef FOLLOW_MEAN
			unbind_mean(&(Cat.Nids),&(Cat.PIDorIndex),s,&Lmain_removed,&PInd_main_removed);
			#endif
			#ifdef FOLLOW_MINPOTH
			unbind_minpotH(&(Cat.Nids),&(Cat.PIDorIndex),s,&Lmain_removed,&PInd_main_removed);	
			#endif
			#ifdef FOLLOW_POTW
			unbind_potW(&(Cat.Nids),&(Cat.PIDorIndex),s,&Lmain_removed,&PInd_main_removed);
			#endif
			#ifdef FOLLOW_CORE
			unbind_core(&(Cat.Nids),&(Cat.PIDorIndex),s,&Lmain_removed,&PInd_main_removed,CoreFrac);		
			#endif
			free(PInd_main_removed);
			if(SubCat.SubLen[0])
			{
			#ifdef FOLLOW_MEAN				
			unbind_mean(SubCat.SubLen+0,SubCat.PSubArr+0,s,&Lmain_removed,&PInd_main_removed);	
			#endif
			#ifdef FOLLOW_MINPOTH
			unbind_minpotH(SubCat.SubLen+0,SubCat.PSubArr+0,s,&Lmain_removed,&PInd_main_removed);	
			#endif
			#ifdef FOLLOW_POTW
			unbind_potW(SubCat.SubLen+0,SubCat.PSubArr+0,s,&Lmain_removed,&PInd_main_removed);	
			#endif
			#ifdef FOLLOW_CORE
			unbind_core(SubCat.SubLen+0,SubCat.PSubArr+0,s,&Lmain_removed,&PInd_main_removed,CoreFrac);
			#endif
			free(PInd_main_removed);
			}
			fprintf(fp,"%d,%d,%d\n",SnapBind,Cat.Nids,SubCat.SubLen[0]);fflush(fp);
			Cat.Nids=GrpLen;
			free(Cat.PIDorIndex);
			Cat.PIDorIndex=mymalloc(sizeof(int)*GrpLen);
			memcpy(Cat.PIDorIndex,GrpPIDs,sizeof(int)*GrpLen);
			for(i=0;i<SubCat.SubLen[0];i++)
			#ifdef PID_ORDERED
				SubCat.PSubArr[0][i]++;
			#else
				SubCat.PSubArr[0][i]=Pdat.PID[SubCat.PSubArr[0][i]];
			#endif
		}
		fclose(fp);	
		return 0;
}
