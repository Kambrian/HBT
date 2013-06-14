/* to generate spatial sub-in-sub hierarchy */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define CEN_COM  //using CoM as center; otherwise using MostBnd (recommended).
//~ #define DYNAMICAL  //dynamical sub-in-sub

#define RHALF 0
#define R1SIG 1
#define R2SIG 2
#define R3SIG 3
#define RPOISSON 4

#define RSUB_TYPE R2SIG 

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
	
SUBCATALOGUE SubCat;
float *RSub;
struct Hierarchy *SubInSub;

void load_RSub(float *rsub,int Nsnap,int Nsub)
{
	char buf[1024];
	FILE *fp;
	int Nsubs,dummy;
	
	if(!rsub)
	{
		printf("error: allocate rsub first \n");
		exit(1);
	}
	#ifdef CEN_COM
	sprintf(buf,"%s/profile/RmaxVmax_%d.COM",SUBCAT_DIR,Nsnap);
	#else
	sprintf(buf,"%s/profile/RmaxVmax_%d.MBD",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(int),1,fp);
	if(Nsub!=Nsubs) 
	{
		printf("error loading %s:\n size not expected %d!=%d\n",buf,Nsub,Nsub);
		exit(1);
	}
	fseek(fp,sizeof(float)*Nsubs*(2+RSUB_TYPE),SEEK_CUR);
	fread(rsub,sizeof(float),Nsubs,fp);
	fseek(fp,sizeof(float)*Nsubs*(5-RSUB_TYPE-1),SEEK_CUR);
	fread(&dummy,sizeof(int),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent %d!=%d\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
}

void init_RSub(int Nsnap)
{
	RSub=mymalloc(sizeof(float)*SubCat.Nsubs);
	load_RSub(RSub,Nsnap,SubCat.Nsubs);
}


void load_SubInSub(int Nsnap)
{
	char buf[1024],outputdir[1024];
	FILE *fp;
	int Nsubs,dummy;

	#ifdef CEN_COM
	sprintf(buf,"%s/profile/SubInSub_%03d.COM.%d",SUBCAT_DIR,Nsnap,RSUB_TYPE);
	#else
	sprintf(buf,"%s/profile/SubInSub_%03d.MBD.%d",SUBCAT_DIR,Nsnap,RSUB_TYPE);
	#endif
	myfopen(fp,buf,"r");
	fread(&Nsubs,sizeof(int),1,fp);
	fread(SubInSub,sizeof(struct Hierarchy),Nsubs,fp);
	fread(&dummy,sizeof(int),1,fp);	
	if(Nsubs!=dummy)
	{
		printf("error read %s, check failed\n",buf);
		exit(1);
	}
	fclose(fp);
}

int main(int argc, char **argv)
{
	int Nsnap,subid;
	logfile=stdout;
	
	if(argc!=3)
	{
		printf("Usage: %s [Nsnap] [subid]\n",argv[0]);
		exit(1);
	}
	Nsnap=atoi(argv[1]);
	subid=atoi(argv[2]);
	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	init_RSub(Nsnap);
	#ifdef DYNAMICAL
	SubInSub=SubCat.sub_hierarchy;
	#else
	SubInSub=mymalloc(sizeof(struct Hierarchy)*SubCat.Nsubs);
	load_SubInSub(Nsnap);
	#endif
	
	FILE *fp;
	char buf[1024];
	sprintf(buf,"%s/anal/image",SUBCAT_DIR);
	mkdir(buf,0755);
	#ifdef DYNAMICAL
	sprintf(buf,"%s/anal/image/SubInSub_S%dB%d.DYN.%d",SUBCAT_DIR,Nsnap,subid,RSUB_TYPE);
	#else
	#ifdef CEN_COM
	sprintf(buf,"%s/anal/image/SubInSub_S%dB%d.COM.%d",SUBCAT_DIR,Nsnap,subid,RSUB_TYPE);
	#else
	sprintf(buf,"%s/anal/image/SubInSub_S%dB%d.MBD.%d",SUBCAT_DIR,Nsnap,subid,RSUB_TYPE);
	#endif
	#endif
	myfopen(fp,buf,"w");
	
	if(subid>0||SubCat.SubLen[subid])
	write_subinsub(subid,0,fp);

	fclose(fp);
}

int write_subinsub(int subid, int level, FILE *fp)
{//output the subhalo plus all its sub-in-subs
	int Nsubs;
	
	fprintf(fp,"%d,%d,%g,%g,%g,%g\n", subid, level, SubCat.Property[subid].CoM[0],
						SubCat.Property[subid].CoM[1],SubCat.Property[subid].CoM[2],
						RSub[subid]);
	Nsubs=1;
	
	subid=SubInSub[subid].sub;					
	while(subid>=0)
	{
		Nsubs+=write_subinsub(subid,level+1,fp);
		subid=SubInSub[subid].next;
	}
	return Nsubs;
}
