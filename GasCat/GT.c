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

int main(int argc,char **argv)
{
char gasdir[1024]=GASCAT_DIR;
char buf[1024],filemode[4];
HBTInt SnapshotNum,i,SnapRange[2];
GASHALOCAT GCat;
GASSUBCAT GSubCat;
GASSRCCAT GSrcCat;
SUBCATALOGUE SubCat;
time_t time_start,time_end; //for program timing

time_start=time(NULL);
/*=================program initialization==================*/
{
	if(argc==1)
	{
		SnapRange[0]=0;
		SnapRange[1]=MaxSnap-1;
	}
	 else if(argc==2)
	 {
		 SnapRange[0]=atoi(argv[1]);
		 SnapRange[1]=atoi(argv[1]);
	 }
	 else if(argc==3)
	 {
		 SnapRange[0]=atoi(argv[1]);
		 SnapRange[1]=atoi(argv[2]);
	 }
	 if((argc>3)||SnapRange[0]<0||SnapRange[0]>=MaxSnap||SnapRange[1]<0||SnapRange[1]>=MaxSnap||SnapRange[1]<SnapRange[0])
	  {
		 printf("Usage: %s [SnapShot_begin [SnapShot_end]]\n\
				\tSnapShot_begin and SnapShot_end being integers in the range 0~MaxSnap-1,specifying the snap range to be run\n\
				\t if no param is specified, run over all the snapshots\n\
				\t if only SnapShot_begin is specified, Snapshot_end is set to SnapShot_begin\n",argv[0]);
		 exit(1);
	 }	
	 
	 if(SnapRange[0]>0)//add logs to existing file
			sprintf(filemode,"%s","a");
	else//create new logfiles
			sprintf(filemode,"%s","w");
			
	if(0==(strcmp(LOGFILE_NAME,"stdout"))) 
		logfile=stdout;
	else
	{	
		sprintf(buf,"%s/logfile.GT",gasdir);
		if((logfile=fopen(buf,filemode))==NULL)	{fprintf(stderr,"Error opening file %s\n",buf);exit(1);}
	}

	if(SnapRange[0]>0)//pick up a previous subcat to continue,assign snapnum first
	{		
			fprintf(logfile,"\nrestarting program %s at %s from Snapshot %d\n",argv[0],ctime(&time_start),SnapRange[0]);
			SnapshotNum=SnapRange[0]-1;
			load_gassrccat(SnapshotNum,&GSrcCat,gasdir);
	}
	else//forge SnapshotNum=-1
	{	
			GSrcCat.Nsubs=0;create_gassrccat(&GSrcCat);
	}	
}

for(SnapshotNum=SnapRange[0];SnapshotNum<=SnapRange[1];SnapshotNum++)
{
	load_gashalocat(SnapshotNum,&GCat,gasdir);
	load_sub_catalogue(SnapshotNum,&SubCat,SUBCAT_DIR);
	load_particle_data(SnapshotNum,SNAPSHOT_DIR);
	load_gas_data(SnapshotNum,SNAPSHOT_DIR);
	fresh_ID2Index(&SubCat,FRSH_SUBCAT);	fresh_gasID2Index(&GCat,FRSH_GRPCAT);	fresh_gasID2Index(&GSrcCat,FRSH_SRCCAT);
	make_gassrccat(&GSrcCat,&GSubCat,&GCat,&SubCat,SnapshotNum);
	for(i=0;i<SubCat.Nsubs;i++)
	unbindgas(GSubCat.SubLen+i,GSubCat.PSubArr+i,GSubCat.Property+i,SubCat.SubLen[i],SubCat.PSubArr[i],SubCat.Property[i].CoM,SubCat.Property[i].VCoM);
	
	complete_N_save_gas(&GSrcCat,&GSubCat,SnapshotNum,gasdir);
	
	erase_gassubcat(&GSubCat);
	free_gashalocat(&GCat);
	erase_sub_catalogue(&SubCat);
	free_gas_data();
}

time_end=time(NULL);
fprintf(logfile,"Program %s Finished in %ld minutes (or %ld hours) (or %ld seconds).\n",argv[0], (time_end-time_start)/60,(time_end-time_start)/3600,time_end-time_start);
fclose(logfile);
return 0;
}

