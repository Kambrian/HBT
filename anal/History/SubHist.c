//print history of a given subhalo, for the mergermajors of subhalosGoNotts
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define LEN(x) ((x)<0?0:SubCat.SubLen[(x)])


int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat;
	CATALOGUE Cat;
	HBTInt Nsnap,subid[2],Nsubs, *pro2dest;
	HBTReal d;
	FILE *fp;
	char buf[1024];
	
	logfile=stdout;
	subid[0]=0;subid[1]=1;
	sprintf(buf,"%s/anal/notts_trace.dat",SUBCAT_DIR);
	myfopen(fp,buf,"w");
	
	for(Nsnap=0;Nsnap<MaxSnap-1;Nsnap++)
	{
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	if(subid[0]<0||subid[1]<0)
		d=0.;
	else
		d=distance(SubCat.Property[subid[0]].CoM,SubCat.Property[subid[1]].CoM);
	fprintf(fp,"%d %d %d %d %d %g\n",Nsnap,subid[0],subid[1], LEN(subid[0]), LEN(subid[1]), d);	
	free_sub_table(&SubCat);
	load_pro2dest(Nsnap,&pro2dest,&Nsubs,SUBCAT_DIR);	
	subid[0]=pro2dest[subid[0]];
	subid[1]=pro2dest[subid[1]];
	free_pro2dest(pro2dest);
	}
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	if(subid[0]<0||subid[1]<0)
		d=0.;
	else
		d=distance(SubCat.Property[subid[0]].CoM,SubCat.Property[subid[1]].CoM);
	fprintf(fp,"%d %d %d %d %d %g\n",Nsnap,subid[0],subid[1], LEN(subid[0]), LEN(subid[1]), d);	
	free_sub_table(&SubCat);
	return 0;
}
