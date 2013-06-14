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
	FILE *fd;
	char buf[1024];
	int Nsnap,Ngroups,Nsubs;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	int Nall=0;
	SUBCATALOGUE SubCat[MaxSnap];
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
		load_sub_table(Nsnap,SubCat+Nsnap,SUBCAT_DIR);
		printf("Snap %d\n",Nsnap);
		Nall+=SubCat[Nsnap].Nsubs;
	}
	printf("Nall=%d,Mem=%2.1gG\n",Nall,Nall*4.0/1024/1024/1024);

	return 0;
}
