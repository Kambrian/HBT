//for tranforming sp2pro file in version <7.8 to v7.8 format
//i.e.,to add Nsplitter and Npro as header, save using Nsnap_dest, also move to sp2pro dir
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
	
int main()
{
SUBCATALOGUE SubCat;
int Nsnap,subid,proid,Npro,Nsplitter;
int *sp2pro;
char buf[1024];
FILE *fp_old,*fp_new;

logfile=stdout;
Npro=0;
Nsplitter=0;
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	sprintf(buf,"%s/sp2pro/sp2pro_%03d",SUBCAT_DIR,Nsnap);//should be moved to splitters later
	myfopen(fp_new,buf,"w");
	fwrite(&SubCat.Nsplitter,sizeof(int),1,fp_new);
	fwrite(&Npro,sizeof(int),1,fp_new);
	if(SubCat.Nsplitter)
	{
	sp2pro=mymalloc(sizeof(int)*SubCat.Nsplitter);
	Nsplitter=0;
	sprintf(buf,"%s/splitters/sp2pro_%03d",SUBCAT_DIR,Nsnap-1);
	myfopen(fp_old,buf,"r");
	while(1)
	{
		fread(sp2pro+Nsplitter,sizeof(int),1,fp_old);
		if(feof(fp_old)) break;
		Nsplitter++;
	}
	fclose(fp_old);
	if(Nsplitter!=SubCat.Nsplitter)
	{
		printf("error: splitter mismatch %d=%d\n",Nsplitter,SubCat.Nsplitter);
		exit(1);
	}	
	fwrite(sp2pro,sizeof(int),Nsplitter,fp_new);
	free(sp2pro);
	}
	fclose(fp_new);
	Npro=SubCat.Nsubs;
	erase_sub_catalogue(&SubCat);
}
return 0;
}
