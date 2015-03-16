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
	HBTInt Nsnap=0;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	HBTInt i,j;
	load_particle_header(Nsnap, SNAPSHOT_DIR);
	IO_HEADER h;
	long np[6]={0}, ntot=0;
	for(i=0;i<header.num_files;i++)
	{
	  printf("%d:",i);
	  load_particle_header_i_into(Nsnap, SNAPSHOT_DIR, &h, i);
	  for(j=0;j<6;j++)
	  {
		np[j]+=h.npart[j];
		printf("%d ", h.npart[j]);
	  }
	  printf("\n");
	}
	printf("Num particles:\n");
	for(j=0;j<6;j++)
	{
	  printf("PartType %d: %ld (%u)\n", j, np[j], header.npartTotal[j]); 
	  ntot+=np[j];
	}
	
	printf("\nTotal=%ld\n", ntot);fflush(stdout);
// 	load_particle_data(Nsnap, SNAPSHOT_DIR);
	CATALOGUE Cat;
	load_group_catalogue(Nsnap, &Cat, GRPCAT_DIR);
	printf("groups loaded\n");
	return 0;
}
