#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

extern HBTInt read_part_pos_JING(int Nsnap,char *snapdir);
extern void test_open(int Nsnap, char *snapdir);
int main(int argc, char** argv)
{
	CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	SRCCATALOGUE SrcCat;
	HBTInt Nsubs,*pro2dest,Npro,Nsplitter,*sp2pro;
	
	HBTInt Nsnap=0;
	HBTInt grpid,subid,pid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
// 	test_open(snaplist[Nsnap], SNAPSHOT_DIR);
// 	read_part_pos_JING(snaplist[Nsnap],SNAPSHOT_DIR);
	
	load_particle_data_bypart(Nsnap, SNAPSHOT_DIR);
	
	pid=2;
	printf("%ld\n", Pdat.PID[pid]);
	int i;
	for(i=0;i<3;i++)
	  printf("%g,", Pdat.Pos[pid][i]);
	puts("\n");
	for(i=0;i<3;i++)
	  printf("%g,", Pdat.Vel[pid][i]);
	puts("\n");
	
	return 0;
}
