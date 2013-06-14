#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"


int main(int nargc,char **argv)
{
SUBCATALOGUE SubCat;
int Nsnap,subid,proid,Npro;
int *pro2dest;

logfile=stdout;
Npro=0;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	Npro+=SubCat.Nsplitter;
	pro2dest=mymalloc(sizeof(int)*Npro);
	for(subid=0;subid<Npro;subid++)
		pro2dest[subid]=-1;
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
		proid=SubCat.HaloChains[subid].ProSubID;
		if(proid>=0&&proid<Npro)
			pro2dest[proid]=subid;
	}
	save_pro2dest(Nsnap-1,pro2dest,Npro,SUBCAT_DIR);
	free(pro2dest);
	Npro=SubCat.Nsubs;
	erase_sub_catalogue(&SubCat);
}
return 0;
}
