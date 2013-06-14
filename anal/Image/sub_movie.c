#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc,char **argv)
{
	CATALOGUE Cat;
	int Nsnap,SnapMin=25,SnapMax=45,i,j,pid;
	int GrpLoad=6280,SnapLoad=30;
	//~ int GrpLoad=6508,SnapLoad=35;
	//~ int GrpLoad=2013,SnapLoad=56;
	int Nsubs,*pro2dest;
	FILE *fp;
	char buf[1024];

	int *PIndex,N;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/image",SUBCAT_DIR);
	mkdir(outputdir,0755);
	logfile=stdout;

//prepare particles
	load_group_catalogue(SnapLoad,&Cat,GRPCAT_DIR);
	
	PIndex=Cat.PIDorIndex+Cat.Offset[GrpLoad];
	N=Cat.Len[GrpLoad];
	#ifdef GRPINPUT_INDEX
	for(i=0;i<N;i++)
		#ifdef PID_ORDERED
		PIndex[i]++;
		#else
		PIndex[i]=Pdat.PID[PIndex[i]];
		#endif
	#endif
	
	for(Nsnap=SnapMin;Nsnap<=SnapMax;Nsnap++)
	{
		load_particle_data(Nsnap,SNAPSHOT_DIR);
		fresh_ID2Index(PIndex,N);

		sprintf(buf,"%s/posmap_S%dG%d.%d",outputdir,SnapLoad,GrpLoad,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<N;i++)
		{
			for(j=0;j<3;j++)
				fprintf(fp,"%g,",Pdat.Pos[PIndex[i]][j]);
			fprintf(fp,"\n");
		}
		fclose(fp);
		
		for(i=0;i<N;i++)
		#ifdef PID_ORDERED
		PIndex[i]++;
		#else
		PIndex[i]=Pdat.PID[PIndex[i]];
		#endif
	}
				
return 0;
}
