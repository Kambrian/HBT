#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define NDIV 200
#include "linkedlist.c"

int main(int argc, char** argv)
{	
	logfile=stdout;//redirect BT routines' log info to standard output	
	HBTInt Nsnap=0;
	
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	printf("Nsnap=%d\t",Nsnap);fflush(stdout);
	makell(Pdat.Pos,NP_DM,NDIV);
	LINKLIST l;
	make_linklist(&l,Pdat.Pos,NP_DM,NDIV);
	printf("%d,%d\n",hoc[10][9][133],hoc[77][2][1]);
	printf("%d,%d\n",linklist_get_hoc(&l,10,9,133),linklist_get_hoc(&l,77,2,1));
	free_linklist(&l);
	free_particle_data();
	
	return 0;
}

