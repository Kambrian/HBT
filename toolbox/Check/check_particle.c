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
	HBTInt Nsnap=MaxSnap-1;
	logfile=stdout;//redirect BT routines' log info to standard output
	
#define NP 3
	HBTInt pids[NP]={6464605614178645814L, 6464588470511166421L, 6464588470511166421L};
	HBTInt pids0[NP];
	memcpy(pids0, pids, sizeof(HBTInt)*NP);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&pids,NP);
	free_PIDHash();

	int i;
	for(i=0;i<NP;i++)
	{
	  HBTReal *x;
	  x=Pdat.Pos[pids[i]];
	  printf("%ld, %g, %g, %g\n", pids0[i], x[0], x[1], x[2]);
	  x=Pdat.Vel[pids[i]];
	  printf("%ld, %g, %g, %g\n", pids0[i], x[0], x[1], x[2]);
	}
	
	free_particle_data();
	return 0;
}
