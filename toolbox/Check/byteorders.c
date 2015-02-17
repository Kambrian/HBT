#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main()
{
char buf[1024];
int i,Nsnap,byteorder;

printf("Snapshots need byte-swap:\n");
for(i=0;i<MaxSnap;i++)
{
	Nsnap=i;
	#ifdef SNAPLIST
	Nsnap=snaplist[i];
	#endif
	sprintf(buf,"%s/snapdir_%03d/%s_%03d.0",SNAPSHOT_DIR,Nsnap,SNAPFILE_BASE,Nsnap);
    if(1==NFILES)
	 if(!try_readfile(buf))	sprintf(buf,"%s/%s_%03d",SNAPSHOT_DIR,SNAPFILE_BASE,Nsnap); //try the other convention
	
	byteorder=check_snapshot_byteorder(buf);	
	
	if(byteorder)
		printf("%d,%d\n",i, Nsnap);
}
printf("END\n");
	return 0;
}
