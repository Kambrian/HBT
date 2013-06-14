/*to fix the bug of not saving subhalo using PIDs*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc,char **argv)
{
char outputdir[1024];
char inputdir[1024]=SUBCAT_DIR;
SUBCATALOGUE SubCat;
int SnapshotNum;
int snapbegin,snapend,subid,i;

if(argc!=3)
{
	printf("usage: %s [snap_begin] [snap_end]\n",argv[0]);
	return 1;
}
snapbegin=atoi(argv[1]);
snapend=atoi(argv[2]);

sprintf(outputdir,"%s/patched2",SUBCAT_DIR);
for(SnapshotNum=snapbegin;SnapshotNum<=snapend;SnapshotNum++)
{
	load_sub_catalogue(SnapshotNum,&SubCat,inputdir);
	load_particle_data(SnapshotNum,SNAPSHOT_DIR);
	for(subid=0;subid<SubCat.Nsubs;subid++)
		for(i=0;i<SubCat.SubLen[subid];i++)
			SubCat.PSubArr[subid][i]=Pdat.PID[SubCat.PSubArr[subid][i]];
	save_sub_catalogue(SnapshotNum,&SubCat,outputdir);
	erase_sub_catalogue(&SubCat);
}
return 0;
}
