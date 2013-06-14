#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define SUBLEN(i) (subid[i]<0?0:SubCat.SubLen[subid[i]])
int sum_SubInSub_mass(int subid,SUBCATALOGUE *SubCat)
{
	int id,mass;
	mass=SubCat->SubLen[subid];
	id=SubCat->sub_hierarchy[subid].sub;
	while(id>=0) //add-up sub-in-sub mass
	{
	mass+=sum_SubInSub_mass(id,SubCat);
	id=SubCat->sub_hierarchy[id].next;
	}
	return mass;
}

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
int i,subid[5],grpid,Nsnap,sid;
int *pro2dest,Nsubs,flagstop;
char buf[1024];
FILE *fp;

logfile=stdout;

Nsnap=36;
grpid=33;
sprintf(buf,"%s/anal/follow/lostsub_%03d_%d",SUBCAT_DIR,Nsnap,grpid);
myfopen(fp,buf,"w");
for(Nsnap=36;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	if(Nsnap==36)
	{
		sid=SubCat.GrpOffset_Sub[grpid];
		for(i=0;i<5;i++)
		subid[i]=sid+i;
	}
	fprintf(fp,"%d,%d,%d,%d,%d,%d\n",Nsnap,SUBLEN(0),SUBLEN(1),SUBLEN(2),SUBLEN(3),SUBLEN(4));
	erase_sub_catalogue(&SubCat);
	if(Nsnap==99) break;
	load_pro2dest(Nsnap,&pro2dest,&Nsubs,SUBCAT_DIR);
	flagstop=1;
	for(i=0;i<5;i++)
	{
		subid[i]=pro2dest[subid[i]];								
		if(subid[i]>0) flagstop=0;
	}
	free_pro2dest(pro2dest);
	if(flagstop) break;
}
fclose(fp);
return 0;
}
	

