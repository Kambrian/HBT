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
	int Nsnap,grpid,subid,HostID;
	
	grpid=10;//the target group
	
	logfile=stdout;
	Nsnap=MaxSnap-1;//final snapshot
	load_sub_table(Nsnap, &SubCat, SUBCAT_DIR);
	subid=SubCat.GrpOffset_Sub[grpid];//the main-sub of the target group
	free_sub_table(&SubCat);
	
	for(Nsnap=MaxSnap-1;Nsnap>=0;Nsnap--)
	{
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR); //load fof catalogue
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR); //load subhalo catalogue (only bulk properties loaded; particle lists not loaded. 
	//use load_sub_catalogue() if you need the particles
	HostID=SubCat.HaloChains[subid].HostID; //the host id of current subhalo
	printf("%d:Host %d, HostLen %d,M0 %d\n",Nsnap,HostID,Cat.Len[HostID],Cat.Len[0]);
	subid=SubCat.HaloChains[subid].ProSubID; //the subid of the progenitor of this subhalo
	free_sub_table(&SubCat);//use erase_sub_catalogue() if you used load_sub_catalogue() above.
	free_catalogue(&Cat);//free fof catalogue
	if(subid<0) break;//stop if no progenitor.
	}
	
	return 0;
}
