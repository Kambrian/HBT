/* to save subhalo data for Subhalos Go Notts project */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#if defined CONVERT_LENGTH_MPC_KPC
#define LUNIT 1.
#else
#define LUNIT 1000.
#endif
#define VUNIT 1
#define MUNIT 1e4

#define RHALF 0
#define R1SIG 1
#define R2SIG 2
#define R3SIG 3
#define RPOISSON 4

#define RSUB_TYPE R2SIG 
#define DYNAMIC_SUBINSUB
#define CEN_COM  //using CoM as center; otherwise using MostBnd (recommended).

#define OUTPUT_ID_ONLY //ignore Rsub,CoM,VCoM,Rmax,Vmax...etc., only output particle ids
	
SUBCATALOGUE SubCat;
HBTReal *RSub, *Vmax, *Rmax, *Rhalf;
struct Hierarchy *SubInSub;
#ifdef HBTPID_RANKSTYLE
IDatInt *PIDs_Sorted;  
#endif

void init_RSub(HBTInt Nsnap);
void init_RmaxVmax(HBTInt Nsnap);
void load_SubInSub(HBTInt Nsnap);
HBTInt sum_subinsubLen(HBTInt subid);
HBTInt write_subIDs(HBTInt subid, FILE *fp);
HBTInt write_subinsubIDs(HBTInt subid, FILE *fp);
int main(int argc, char** argv)
{
	HBTInt i;
	HBTInt Nsnap;
	HBTInt grpid,subid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	Nsnap=MaxSnap-1;
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap="HBTIFMT"\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	//~ for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	//~ {
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	#ifdef HBTPID_RANKSTYLE
	PIDs_Sorted=load_PIDs_Sorted();
	#endif
	//~ fill_PIDHash();
	//~ fresh_ID2Index(&SubCat,-2);	
	//~ free_PIDHash();
	#ifndef OUTPUT_ID_ONLY
	init_RSub(Nsnap);
	init_RmaxVmax(Nsnap);
	#endif
	#ifdef DYNAMIC_SUBINSUB
	SubInSub=SubCat.sub_hierarchy;
	#else
	SubInSub=mymalloc(sizeof(struct Hierarchy)*SubCat.Nsubs);
	load_SubInSub(Nsnap);
	#endif
	
	printf("Nsnap=%d, h=%f\n",Nsnap, header.HubbleParam);
	
	FILE *fptab,*fpid,*fpid2;
	char buf[1024];
	sprintf(buf,"%s/anal/notts_tab_"HBTIFMT".txt",SUBCAT_DIR,Nsnap);
	myfopen(fptab,buf,"w");
	sprintf(buf,"%s/anal/notts_ids_"HBTIFMT".txt",SUBCAT_DIR,Nsnap);
	myfopen(fpid,buf,"w");
	#ifdef DYNAMIC_SUBINSUB
	sprintf(buf,"%s/anal/notts_iddup_"HBTIFMT".dynamic.txt",SUBCAT_DIR,Nsnap);
	#else
	sprintf(buf,"%s/anal/notts_iddup_"HBTIFMT".txt",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fpid2,buf,"w");
	
	for(grpid=0;grpid<SubCat.Ngroups;grpid++)
	{
		if(1==SubCat.GrpLen_Sub[grpid]&&0==SubCat.SubLen[SubCat.GrpOffset_Sub[grpid]]) continue;
		fprintf(fptab,"#FoFHalo #"HBTIFMT", "HBTIFMT" subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		fprintf(fpid,"#FoFHalo #"HBTIFMT", "HBTIFMT" subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		fprintf(fpid2,"#FoFHalo #"HBTIFMT", "HBTIFMT" subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		for(subid=SubCat.GrpOffset_Sub[grpid];subid<SubCat.GrpOffset_Sub[grpid]+SubCat.GrpLen_Sub[grpid];subid++)
		{
			#ifdef OUTPUT_ID_ONLY
			fprintf(fptab,""HBTIFMT" 0 0 0 0 0 0 %g 0 0 0\n",subid,
			//SubCat.Property[subid].CoM[0]*LUNIT,SubCat.Property[subid].CoM[1]*LUNIT,SubCat.Property[subid].CoM[2]*LUNIT,
			//SubCat.Property[subid].VCoM[0]*VUNIT,SubCat.Property[subid].VCoM[1]*VUNIT,SubCat.Property[subid].VCoM[2]*VUNIT,
			SubCat.SubLen[subid]*header.mass[1]*MUNIT);
			#else
			fprintf(fptab,""HBTIFMT" %g %g %g %g %g %g %g %g %g %g\n",subid,
			SubCat.Property[subid].CoM[0]*LUNIT,SubCat.Property[subid].CoM[1]*LUNIT,SubCat.Property[subid].CoM[2]*LUNIT,
			SubCat.Property[subid].VCoM[0]*VUNIT,SubCat.Property[subid].VCoM[1]*VUNIT,SubCat.Property[subid].VCoM[2]*VUNIT,
			SubCat.SubLen[subid]*header.mass[1]*MUNIT,
			RSub[subid]*LUNIT,Rmax[subid]*LUNIT,Vmax[subid]*VUNIT);
			#endif
			
			fprintf(fpid,""HBTIFMT" "HBTIFMT"\n", subid, SubCat.SubLen[subid]);
			write_subIDs(subid,fpid);
			
			fprintf(fpid2,""HBTIFMT" "HBTIFMT"\n", subid, sum_subinsubLen(subid));
			write_subinsubIDs(subid,fpid2);
		}
	}
	//Quasi-halos
	fprintf(fptab,"#FieldSubs, "HBTIFMT" subhalos\n",SubCat.NQuasi);
	fprintf(fpid,"#FieldSubs, "HBTIFMT" subhalos\n",SubCat.NQuasi);
	fprintf(fpid2,"#FieldSubs, "HBTIFMT" subhalos\n",SubCat.NQuasi);
	for(subid=SubCat.Nsubs-SubCat.NQuasi;subid<SubCat.Nsubs;subid++)
	{
		#ifdef OUTPUT_ID_ONLY
			fprintf(fptab,""HBTIFMT" 0 0 0 0 0 0 %g 0 0 0\n",subid,
			//SubCat.Property[subid].CoM[0]*LUNIT,SubCat.Property[subid].CoM[1]*LUNIT,SubCat.Property[subid].CoM[2]*LUNIT,
			//SubCat.Property[subid].VCoM[0]*VUNIT,SubCat.Property[subid].VCoM[1]*VUNIT,SubCat.Property[subid].VCoM[2]*VUNIT,
			SubCat.SubLen[subid]*header.mass[1]*MUNIT);
			#else
			fprintf(fptab,""HBTIFMT" %g %g %g %g %g %g %g %g %g %g\n",subid,
			SubCat.Property[subid].CoM[0]*LUNIT,SubCat.Property[subid].CoM[1]*LUNIT,SubCat.Property[subid].CoM[2]*LUNIT,
			SubCat.Property[subid].VCoM[0]*VUNIT,SubCat.Property[subid].VCoM[1]*VUNIT,SubCat.Property[subid].VCoM[2]*VUNIT,
			SubCat.SubLen[subid]*header.mass[1]*MUNIT,
			RSub[subid]*LUNIT,Rmax[subid]*LUNIT,Vmax[subid]*VUNIT);
			#endif
		
		fprintf(fpid,""HBTIFMT" "HBTIFMT"\n", subid, SubCat.SubLen[subid]);
		write_subIDs(subid,fpid);
		
		fprintf(fpid2,""HBTIFMT" "HBTIFMT"\n", subid, sum_subinsubLen(subid));
		write_subinsubIDs(subid,fpid2);
	}
	fclose(fptab);
	fclose(fpid);
	fclose(fpid2);

	erase_sub_catalogue(&SubCat);
	myfree(RSub);
	myfree(Rmax);
	myfree(Vmax);
	myfree(Rhalf);
	#ifndef DYNAMIC_SUBINSUB
	myfree(SubInSub);
	#endif
	#ifdef HBTPID_RANKSTYLE
	myfree(PIDs_Sorted);
	#endif
	free_particle_data();
	//~ }
	
	return 0;
}


void load_RSub(HBTReal *rsub,HBTInt Nsnap,HBTInt Nsub)
{
	char buf[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;
	
	if(!rsub)
	{
		printf("error: allocate rsub first \n");
		exit(1);
	}
	#ifdef CEN_COM
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
	#else
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	if(Nsub!=Nsubs) 
	{
		printf("error loading %s:\n size not expected "HBTIFMT"!="HBTIFMT"\n",buf,Nsub,Nsub);
		exit(1);
	}
	fseek(fp,sizeof(HBTReal)*Nsubs*(2+RSUB_TYPE),SEEK_CUR);
	fread(rsub,sizeof(HBTReal),Nsubs,fp);
	fseek(fp,sizeof(HBTReal)*Nsubs*(5-RSUB_TYPE-1),SEEK_CUR);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
}
int load_RmaxVmax(HBTReal *rmax,HBTReal *vmax, HBTReal *rhalf, HBTInt Nsnap)
{
	char buf[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;
	//~ HBTReal *rsig, *r2sig, *r3sig, *rpoisson; /* declare this as input var if you want them!!! */
	
	if(!rmax||!vmax||!rhalf)
	{
		printf("error: allocate rmax , vmax and rhalf first \n");
		exit(1);
	}
	#ifdef CEN_COM
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
	#else
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	fread(rmax,sizeof(HBTReal),Nsubs,fp);
	fread(vmax,sizeof(HBTReal),Nsubs,fp);
	fread(rhalf,sizeof(HBTReal),Nsubs,fp);
	fseek(fp,sizeof(HBTReal)*Nsubs*4,SEEK_CUR);
	//~ fread(rsig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r2sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r3sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(rpoisson,sizeof(HBTReal),Nsubs,fp);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
	return Nsubs;
}
void init_RSub(HBTInt Nsnap)
{
	RSub=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	load_RSub(RSub,Nsnap,SubCat.Nsubs);
}
void init_RmaxVmax(HBTInt Nsnap)
{
	Rmax=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	Vmax=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	Rhalf=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	load_RmaxVmax(Rmax,Vmax,Rhalf,Nsnap);
}
void load_SubInSub(HBTInt Nsnap)
{
	char buf[1024],outputdir[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;

	#ifdef CEN_COM
	sprintf(buf,"%s/profile/SubInSub_%03d.COM.%d",SUBCAT_DIR,Nsnap,RSUB_TYPE);
	#else
	sprintf(buf,"%s/profile/SubInSub_%03d.MBD.%d",SUBCAT_DIR,Nsnap,RSUB_TYPE);
	#endif
	myfopen(fp,buf,"r");
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	fread(SubInSub,sizeof(struct Hierarchy),Nsubs,fp);
	fread(&dummy,sizeof(HBTInt),1,fp);	
	if(Nsubs!=dummy)
	{
		printf("error read %s, check failed\n",buf);
		exit(1);
	}
	fclose(fp);
}

HBTInt sum_subinsubLen(HBTInt subid)
{ //total SubLen including all sub-in-subs
	HBTInt Nids;
	Nids=SubCat.SubLen[subid];
	
	subid=SubInSub[subid].sub;					
	while(subid>=0)
	{
		Nids+=sum_subinsubLen(subid);
		subid=SubInSub[subid].next;
	}
	return Nids;
}
HBTInt write_subIDs(HBTInt subid, FILE *fp)
{
	HBTInt i;
	for(i=0;i<SubCat.SubLen[subid];i++)
	{
		//~ fprintf(fp,""HBTIFMT"\n",Pdat.PID[SubCat.PSubArr[subid][i]]);
		#ifdef HBTPID_RANKSTYLE
		#ifdef INPUT_INT8
		fprintf(fp,"%lld\n",PIDs_Sorted[SubCat.PSubArr[subid][i]]);
		#else
		fprintf(fp,""HBTIFMT"\n",PIDs_Sorted[SubCat.PSubArr[subid][i]]);
		#endif
		#else
		#ifdef HBT_INT8
		fprintf(fp,"%lld\n",SubCat.PSubArr[subid][i]);
		#else
		fprintf(fp,""HBTIFMT"\n",SubCat.PSubArr[subid][i]);
		#endif
		#endif
	}
	return SubCat.SubLen[subid];
}
HBTInt write_subinsubIDs(HBTInt subid, FILE *fp)
{//output the subhalo plus all its sub-in-subs
	HBTInt Nids;
	Nids=write_subIDs(subid,fp);
	
	subid=SubInSub[subid].sub;					
	while(subid>=0)
	{
		Nids+=write_subinsubIDs(subid,fp);
		subid=SubInSub[subid].next;
	}
	return Nids;
}

