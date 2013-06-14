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

int main(int nargc,char **argv)
{
int SnapLoad,SnapBind;
CATALOGUE Cat;
//~ SUBCATALOGUE SubCat;
int GrpLen, *GrpPIDs,Lmain_removed,*PInd_main_removed,GrpLen2,*GrpPIDs2,flag_Grp2=0;
int fofid,i,snap_skip;
float s[3];
FILE *fp;
char buf[1024];
float CoreFrac, MassRelax_Factor;

//~ SnapLoad=51;fofid=78;
logfile=stdout;	
SnapLoad=atoi(argv[1]);
fofid=atoi(argv[2]);
MassRelax_Factor=atof(argv[3]);//5,3,2,1.4 ;(0.04,0.1111,0.25,0.51)
snap_skip=atoi(argv[4]);

CoreFrac=1.0/MassRelax_Factor/MassRelax_Factor;

if(snap_skip==1)
sprintf(buf, "%s/anal/follow/follow_%03d_%d.adpt%2.1f",SUBCAT_DIR,SnapLoad,fofid,MassRelax_Factor);
else
sprintf(buf, "%s/anal/follow/follow_%03d_%d.adpt%2.1f_skip%d",SUBCAT_DIR,SnapLoad,fofid,MassRelax_Factor,snap_skip);
if(!(fp= fopen(buf, "w")))
{
fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
exit(1);
}

		load_group_catalogue(SnapLoad,&Cat,GRPCAT_DIR);
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
		/*==initialize GrpPIDs2 ==*/
		GrpLen2=GrpLen;
		GrpPIDs2=mymalloc(sizeof(int)*GrpLen2);
		memcpy(GrpPIDs2,GrpPIDs,sizeof(int)*GrpLen2);
		/*==compact Cat to contain just that fof==*/
		free_catalogue(&Cat);
		Cat.Nids=GrpLen;
		Cat.PIDorIndex=mymalloc(sizeof(int)*GrpLen);
		memcpy(Cat.PIDorIndex,GrpPIDs,sizeof(int)*GrpLen);
		
		for(SnapBind=SnapLoad;SnapBind<MaxSnap;SnapBind+=snap_skip)
		{
			load_particle_data(SnapBind,SNAPSHOT_DIR);	
			fresh_ID2Index(Cat.PIDorIndex,Cat.Nids); 	
			unbind_core(&(Cat.Nids),&(Cat.PIDorIndex),s,&Lmain_removed,&PInd_main_removed,CoreFrac);	
			//~ printf("CoreFrac: %f\n",CoreFrac);
			free(PInd_main_removed);
		
			if(flag_Grp2)  
			{
				if(Cat.Nids>GrpLen2) //update GrpPIDs2 when subhalo has grown big enough
				{
				GrpLen2=Cat.Nids;
				free(GrpPIDs2);
				GrpPIDs2=mymalloc(sizeof(int)*GrpLen2);
				for(i=0;i<GrpLen2;i++)
				#ifdef PID_ORDERED
					GrpPIDs2[i]=Cat.PIDorIndex[i]+1;
				#else
					GrpPIDs2[i]=Pdat.PID[Cat.PIDorIndex[i]];
				#endif
				CoreFrac=((float)GrpLen2)/((float)GrpLen)/MassRelax_Factor;//?? necessary to change this?
		 		}	
				else if(Cat.Nids<GrpLen*CoreFrac)//update GrpPIDs with GrpPIDs2, and update GrpPIDs2 with the current sub (or its nearest progenitor?)
				{
				GrpLen=GrpLen2;
				free(GrpPIDs);
				GrpPIDs=GrpPIDs2;
				GrpLen2=Cat.Nids;
				GrpPIDs2=mymalloc(sizeof(int)*GrpLen2);
				for(i=0;i<GrpLen2;i++)
				#ifdef PID_ORDERED
					GrpPIDs2[i]=Cat.PIDorIndex[i]+1;
				#else
					GrpPIDs2[i]=Pdat.PID[Cat.PIDorIndex[i]];
				#endif
				CoreFrac=((float)GrpLen2)/((float)GrpLen)/MassRelax_Factor;/*every time when GrpLen2 changes, 
																		* update CoreFrac to ensure we will be 
																		* operating in the same range ,
																		* that is, the lower limit sub will be unbind() with
																		* near real CoreFrac, rather than a possibly much higher
																		* CoreFrac if not updated*/
				}
			}
			else if(Cat.Nids<=GrpLen*CoreFrac*MassRelax_Factor)// register GrpPIDs2 when subhalo has reduced enough, and adapt CoreFrac
			{
				GrpLen2=Cat.Nids;
				free(GrpPIDs2);
				GrpPIDs2=mymalloc(sizeof(int)*GrpLen2);
				for(i=0;i<GrpLen2;i++)
				#ifdef PID_ORDERED
					GrpPIDs2[i]=Cat.PIDorIndex[i]+1;
				#else
					GrpPIDs2[i]=Pdat.PID[Cat.PIDorIndex[i]];
				#endif
				flag_Grp2=1;
				CoreFrac=((float)GrpLen2)/((float)GrpLen)/MassRelax_Factor;//?? necessary to change this?
			}
			fprintf(fp,"%d,%d,%d,%d,%.3f\n",SnapBind,Cat.Nids,GrpLen2,GrpLen,CoreFrac);fflush(fp);
			Cat.Nids=GrpLen;
			free(Cat.PIDorIndex);
			Cat.PIDorIndex=mymalloc(sizeof(int)*GrpLen);
			memcpy(Cat.PIDorIndex,GrpPIDs,sizeof(int)*GrpLen);
		}
		fclose(fp);	
		return 0;
}
