/* to save subhalo data for Subhalos Go Notts project */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define LUNIT (1000./header.HubbleParam)
#define VUNIT 1
#define MUNIT (1e4/header.HubbleParam)

#define RHALF 0
#define R1SIG 1
#define R2SIG 2
#define R3SIG 3
#define RPOISSON 4

#define RSUB_TYPE R2SIG 

#define CEN_COM  //using CoM as center; otherwise using MostBnd (recommended).
//~ #define DYNAMICAL  //dynamical sub-in-sub
	
SUBCATALOGUE SubCat;
float *RSub;
struct Hierarchy *SubInSub;

void init_RSub(int Nsnap);
void load_SubInSub(int Nsnap);
int sum_subinsubLen(int subid);
int write_subIDs(int subid, FILE *fp);
int write_subinsubIDs(int subid, FILE *fp);
int main(int argc, char** argv)
{
	int i;
	int Nsnap;
	int grpid,subid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	Nsnap=MaxSnap-1;
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fresh_ID2Index(&SubCat,-2);	
	init_RSub(Nsnap);
	SubInSub=mymalloc(sizeof(struct Hierarchy)*SubCat.Nsubs);
	load_SubInSub(Nsnap);
	
	printf("h=%f\n",header.HubbleParam);
	
	FILE *fptab,*fpid,*fpid2;
	char buf[1024];
	sprintf(buf,"%s/anal/notts_tab2_%d.txt",SUBCAT_DIR,Nsnap);
	myfopen(fptab,buf,"w");
	
	for(grpid=0;grpid<SubCat.Ngroups;grpid++)
	{
		//~ fprintf(fptab,"#FoFHalo #%d, %d subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		//~ fprintf(fpid,"#FoFHalo #%d, %d subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		//~ fprintf(fpid2,"#FoFHalo #%d, %d subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		for(subid=SubCat.GrpOffset_Sub[grpid];subid<SubCat.GrpOffset_Sub[grpid]+SubCat.GrpLen_Sub[grpid];subid++)
		{
			fprintf(fptab,"%d %g %g %g %g %g %g %g %g %d %d\n",subid,
			SubCat.Property[subid].CoM[0]*LUNIT,SubCat.Property[subid].CoM[1]*LUNIT,SubCat.Property[subid].CoM[2]*LUNIT,
			SubCat.Property[subid].VCoM[0]*VUNIT,SubCat.Property[subid].VCoM[1]*VUNIT,SubCat.Property[subid].VCoM[2]*VUNIT,
			SubCat.SubLen[subid]*header.mass[1]*MUNIT,
			RSub[subid]*LUNIT,SubCat.SubLen[subid],sum_subinsubLen(subid));
		}
	}
	//Quasi-halos
	//~ fprintf(fptab,"#FieldSubs, %d subhalos\n",SubCat.NQuasi);
	//~ fprintf(fpid,"#FieldSubs, %d subhalos\n",SubCat.NQuasi);
	//~ fprintf(fpid2,"#FieldSubs, %d subhalos\n",SubCat.NQuasi);
	for(subid=SubCat.Nsubs-SubCat.NQuasi;subid<SubCat.Nsubs;subid++)
	{
		fprintf(fptab,"%d %g %g %g %g %g %g %g %g %d %d\n",subid,
			SubCat.Property[subid].CoM[0]*LUNIT,SubCat.Property[subid].CoM[1]*LUNIT,SubCat.Property[subid].CoM[2]*LUNIT,
			SubCat.Property[subid].VCoM[0]*VUNIT,SubCat.Property[subid].VCoM[1]*VUNIT,SubCat.Property[subid].VCoM[2]*VUNIT,
			SubCat.SubLen[subid]*header.mass[1]*MUNIT,
			RSub[subid]*LUNIT,SubCat.SubLen[subid],sum_subinsubLen(subid));
	}
	fclose(fptab);

	erase_sub_catalogue(&SubCat);
	myfree(RSub);
	myfree(SubInSub);
	
	return 0;
}


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

int sum_subinsubLen(int subid)
{ //total SubLen including all sub-in-subs
	int Nids;
	Nids=SubCat.SubLen[subid];
	
	subid=SubInSub[subid].sub;					
	while(subid>=0)
	{
		Nids+=sum_subinsubLen(subid);
		subid=SubInSub[subid].next;
	}
	return Nids;
}
int write_subIDs(int subid, FILE *fp)
{
	int i;
	for(i=0;i<SubCat.SubLen[subid];i++)
	{
		fprintf(fp,"%d\n",Pdat.PID[SubCat.PSubArr[subid][i]]);
	}
	return SubCat.SubLen[subid];
}
int write_subinsubIDs(int subid, FILE *fp)
{//output the subhalo plus all its sub-in-subs
	int Nids;
	Nids=write_subIDs(subid,fp);
	
	subid=SubInSub[subid].sub;					
	while(subid>=0)
	{
		Nids+=write_subinsubIDs(subid,fp);
		subid=SubInSub[subid].next;
	}
	return Nids;
}

