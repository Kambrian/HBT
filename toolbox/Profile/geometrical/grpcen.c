#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
//~ void erase_sub_catalogue(SUBCATALOGUE *SubCat)
//~ {
	//~ int i;
	//~ for(i=0;i<SubCat->Nsubs;i++)
	//~ {
		//~ myfree(SubCat->PSubArr[i]);
	//~ }
	//~ free_sub_catalogue(SubCat);
//~ }
int main()
{
int Nsnap;
SUBCATALOGUE SubCat;
int grpid;
size_t recordsize;

char buf[1024]; FILE *fp; 
	
logfile=stdout;	

for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	sprintf(buf,"%s/anal/grpcen/grpcen.b20.%d.%04d",SUBCAT_DIR,RUN_NUM,snaplist[Nsnap]);
	myfopen(fp,buf,"w");	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR); 	
	
	recordsize=sizeof(float)*3*SubCat.Ngroups;
	fwrite(&recordsize,sizeof(int),1,fp);
	for(grpid=0;grpid<SubCat.Ngroups;grpid++)
		//~ fwrite(SubCat.CoM[grpid],sizeof(float),3,fp);
		fwrite(SubCat.Property[grpid].CoM,sizeof(float),3,fp);
	fwrite(&recordsize,sizeof(int),1,fp);
	fclose(fp);
	erase_sub_catalogue(&SubCat);
}

return 0;
}
