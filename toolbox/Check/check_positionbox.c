#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc, char** argv)
{	
	logfile=stdout;//redirect BT routines' log info to standard output	
	HBTInt Nsnap=0;
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
	load_particle_data_bypart(Nsnap,SNAPSHOT_DIR,FLAG_LOAD_POS);
	printf("Nsnap=%d\t",Nsnap);fflush(stdout);
	HBTReal box[2][3];
	HBTInt i,j;
	for(i=0;i<3;i++)
		box[0][i]=box[1][i]=Pdat.Pos[0][i];
	#pragma omp parallel for private(j,i)
	for(j=0;j<3;j++)
	for(i=0;i<NP_DM;i++)
	{
		if(Pdat.Pos[i][j]<box[0][j]) box[0][j]=Pdat.Pos[i][j];
		else if(Pdat.Pos[i][j]>box[1][j]) box[1][j]=Pdat.Pos[i][j];
	}
	for(i=0;i<2;i++)
	for(j=0;j<3;j++)
		box[i][j]/=BOXSIZE;
	printf("%.3f,%.3f,%.3f;%.3f,%.3f,%.3f\n",
		box[0][0],box[0][1],box[0][2],box[1][0],box[1][1],box[1][2]);
	fflush(stdout);
	free_particle_data();
	}
	
	return 0;
}
