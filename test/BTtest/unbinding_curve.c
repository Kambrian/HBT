#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#define LOG_NBOUND
#include "binding.c"

int main(int argc, char** argv)
{
	CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	SRCCATALOGUE SrcCat;
	
	float CoreFrac=0.25;
	int Nsnap=0;
	int grpid,subid,pid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=3)
	{printf("usage: %s [Nsnap],[subid]\n",argv[0]);fflush(stdout);}
	else
	{
	Nsnap=atoi(argv[1]);
	subid=atoi(argv[2]);
	//~ grpid=atoi(argv[2]);
	}
	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_src_catalogue(Nsnap,&SrcCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fresh_ID2Index(&Cat,-1);
	fresh_ID2Index(&SrcCat,-3);
	int Len,*PIDs,L_removed,*PInd_removed;
	float s[3];
	Len=SrcCat.SubLen[subid];
	PIDs=SrcCat.PSubArr[subid];
	unbind_core(&(Len),&(PIDs),s,&L_removed,&PInd_removed,CoreFrac);		
	//~ Len=Cat.Len[grpid];
	//~ PIDs=mymalloc(sizeof(int)*Len);
	//~ memcpy(PIDs,Cat.PIDorIndex+Cat.Offset[grpid],sizeof(int)*Len);
	//~ unbind_core(&(Len),&(PIDs),s,&L_removed,&PInd_removed,CoreFrac);		
	//~ erase_src_catalogue(&SrcCat);
	return 0;
}
