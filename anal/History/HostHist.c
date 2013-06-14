//print history of host haloids  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"


int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat;
	CATALOGUE Cat;
	int Nsnap,subid,HostID;
	
	logfile=stdout;
	subid=0;
	for(Nsnap=99;Nsnap>=0;Nsnap--)
	{
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	HostID=SubCat.HaloChains[subid].HostID;
	printf("%d:Host %d, HostLen %d,M0 %d\n",Nsnap,HostID,Cat.Len[HostID],Cat.Len[0]);
	subid=SubCat.HaloChains[subid].ProSubID;
	erase_sub_catalogue(&SubCat);
	free_catalogue(&Cat);
	if(subid<0) break;
	}
	
	return 0;
}
