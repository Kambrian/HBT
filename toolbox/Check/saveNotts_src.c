/* to save srchalo data for Subhalos Go Notts project */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

//~ #define OutLen2  //whether to output SrcLen or SrcLen2

#ifndef OutLen2
#define SLEN SubLen
#define PARR PSubArr
#else
#define SLEN SubLen2
#define PARR PSubArr2
#endif

#ifdef CONVERT_LENGTH_MPC_KPC
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

#define CEN_COM  //using CoM as center; otherwise using MostBnd (recommended).
//~ #define DYNAMICAL  //dynamical sub-in-sub
	
SUBCATALOGUE SubCat;
SRCCATALOGUE SrcCat;
HBTReal *RSub;
struct Hierarchy *SubInSub;
#ifdef HBTPID_RANKSTYLE
IDatInt *PIDs_Sorted;  
#endif

void init_RSub(HBTInt Nsnap);
void load_SubInSub(HBTInt Nsnap);
HBTInt sum_srcinsrcLen(HBTInt subid);
HBTInt write_srcIDs(HBTInt subid, FILE *fp);
HBTInt write_srcinsrcIDs(HBTInt subid, FILE *fp);
int main(int argc, char** argv)
{
	HBTInt i;
	HBTInt Nsnap;
	HBTInt grpid,subid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	Nsnap=MaxSnap-1;
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_src_catalogue(Nsnap,&SrcCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	#ifdef HBTPID_RANKSTYLE
	PIDs_Sorted=load_PIDs_Sorted();
	#endif
	/* do not refresh!! because you need to save the IDs, so just work with the IDs*/
	//~ fill_PIDHash();
	//~ fresh_ID2Index(&SubCat,FRSH_SUBCAT);	
	//~ fresh_ID2Index(&SrcCat,FRSH_SRCCAT);	
	//~ free_PIDHash();
	init_RSub(Nsnap);
	SubInSub=mymalloc(sizeof(struct Hierarchy)*SubCat.Nsubs);
	load_SubInSub(Nsnap);
	
	printf("h=%f\n",header.HubbleParam);
	
	FILE *fptab,*fpid,*fpid2;
	char buf[1024];
	#ifndef OutLen2
	sprintf(buf,"%s/anal/notts_srctab_%d.txt",SUBCAT_DIR,Nsnap);
	myfopen(fptab,buf,"w");
	sprintf(buf,"%s/anal/notts_srcids_%d.txt",SUBCAT_DIR,Nsnap);
	myfopen(fpid,buf,"w");
	sprintf(buf,"%s/anal/notts_srciddup_%d.txt",SUBCAT_DIR,Nsnap);  //note this does not make much sense
	myfopen(fpid2,buf,"w");
	#else
	sprintf(buf,"%s/anal/notts_srctab2_%d.txt",SUBCAT_DIR,Nsnap);
	myfopen(fptab,buf,"w");
	sprintf(buf,"%s/anal/notts_srcids2_%d.txt",SUBCAT_DIR,Nsnap);
	myfopen(fpid,buf,"w");
	sprintf(buf,"%s/anal/notts_srciddup2_%d.txt",SUBCAT_DIR,Nsnap);
	myfopen(fpid2,buf,"w");
	#endif
	
	for(grpid=0;grpid<SubCat.Ngroups;grpid++)
	{
//		if(1==SubCat.GrpLen_Sub[grpid]&&0==SubCat.SubLen[SubCat.GrpOffset_Sub[grpid]]) continue;
		fprintf(fptab,"#FoFHalo #%d, %d subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		fprintf(fpid,"#FoFHalo #%d, %d subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		fprintf(fpid2,"#FoFHalo #%d, %d subhalos\n",grpid,SubCat.GrpLen_Sub[grpid]);
		for(subid=SubCat.GrpOffset_Sub[grpid];subid<SubCat.GrpOffset_Sub[grpid]+SubCat.GrpLen_Sub[grpid];subid++)
		{
			fprintf(fptab,"%d %g %g %g %g %g %g %g %g\n",subid,
			SubCat.Property[subid].CoM[0]*LUNIT,SubCat.Property[subid].CoM[1]*LUNIT,SubCat.Property[subid].CoM[2]*LUNIT,
			SubCat.Property[subid].VCoM[0]*VUNIT,SubCat.Property[subid].VCoM[1]*VUNIT,SubCat.Property[subid].VCoM[2]*VUNIT,
			SrcCat.SLEN[subid]*header.mass[1]*MUNIT,
			RSub[subid]*LUNIT);
			
			fprintf(fpid,"%d %d\n", subid, SrcCat.SLEN[subid]);
			write_srcIDs(subid,fpid);
			
			fprintf(fpid2,"%d %d\n", subid, sum_srcinsrcLen(subid));
			write_srcinsrcIDs(subid,fpid2);
		}
	}
	//Quasi-halos
	fprintf(fptab,"#FieldSubs, %d subhalos\n",SubCat.NQuasi);
	fprintf(fpid,"#FieldSubs, %d subhalos\n",SubCat.NQuasi);
	fprintf(fpid2,"#FieldSubs, %d subhalos\n",SubCat.NQuasi);
	for(subid=SubCat.Nsubs-SubCat.NQuasi;subid<SubCat.Nsubs;subid++)
	{
		fprintf(fptab,"%d %g %g %g %g %g %g %g %g\n",subid,
		SubCat.Property[subid].CoM[0]*LUNIT,SubCat.Property[subid].CoM[1]*LUNIT,SubCat.Property[subid].CoM[2]*LUNIT,
		SubCat.Property[subid].VCoM[0]*VUNIT,SubCat.Property[subid].VCoM[1]*VUNIT,SubCat.Property[subid].VCoM[2]*VUNIT,
		SrcCat.SLEN[subid]*header.mass[1]*MUNIT,
		RSub[subid]*LUNIT);
		
		fprintf(fpid,"%d %d\n", subid, SrcCat.SLEN[subid]);
		write_srcIDs(subid,fpid);
		
		fprintf(fpid2,"%d %d\n", subid, sum_srcinsrcLen(subid));
		write_srcinsrcIDs(subid,fpid2);
	}
	fclose(fptab);
	fclose(fpid);
	fclose(fpid2);

	erase_sub_catalogue(&SubCat);
	erase_src_catalogue(&SrcCat);
	myfree(RSub);
	myfree(SubInSub);
	#ifdef HBTPID_RANKSTYLE
	myfree(PIDs_Sorted);
	#endif
	free_particle_data();
	
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
	sprintf(buf,"%s/profile/RmaxVmax_%d.COM",SUBCAT_DIR,Nsnap);
	#else
	sprintf(buf,"%s/profile/RmaxVmax_%d.MBD",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	if(Nsub!=Nsubs) 
	{
		printf("error loading %s:\n size not expected %d!=%d\n",buf,Nsub,Nsub);
		exit(1);
	}
	fseek(fp,sizeof(HBTReal)*Nsubs*(2+RSUB_TYPE),SEEK_CUR);
	fread(rsub,sizeof(HBTReal),Nsubs,fp);
	fseek(fp,sizeof(HBTReal)*Nsubs*(5-RSUB_TYPE-1),SEEK_CUR);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent %d!=%d\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
}

void init_RSub(HBTInt Nsnap)
{
	RSub=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	load_RSub(RSub,Nsnap,SubCat.Nsubs);
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

HBTInt sum_srcinsrcLen(HBTInt subid)
{ //total SubLen including all sub-in-subs
	HBTInt Nids;
	Nids=SrcCat.SLEN[subid];
	
	subid=SubInSub[subid].sub;					
	while(subid>=0)
	{
		Nids+=sum_srcinsrcLen(subid);
		subid=SubInSub[subid].next;
	}
	return Nids;
}
HBTInt write_srcIDs(HBTInt subid, FILE *fp)
{
	HBTInt i;
	for(i=0;i<SrcCat.SLEN[subid];i++)
	{
		//~ fprintf(fp,"%d\n",Pdat.PID[SubCat.PSubArr[subid][i]]);
		#ifdef HBTPID_RANKSTYLE
		#ifdef INPUT_INT8
		fprintf(fp,"%lld\n",PIDs_Sorted[SrcCat.PARR[subid][i]]);
		#else
		fprintf(fp,"%d\n",PIDs_Sorted[SrcCat.PARR[subid][i]]);
		#endif
		#else
		#ifdef HBT_INT8
		fprintf(fp,"%lld\n",SrcCat.PARR[subid][i]);
		#else
		fprintf(fp,"%d\n",SrcCat.PARR[subid][i]);
		#endif
		#endif
	}
	return SrcCat.SLEN[subid];
}
HBTInt write_srcinsrcIDs(HBTInt subid, FILE *fp)
{//output the subhalo plus all its sub-in-subs
	HBTInt Nids;
	Nids=write_srcIDs(subid,fp);
	
	subid=SubInSub[subid].sub;					
	while(subid>=0)
	{
		Nids+=write_srcinsrcIDs(subid,fp);
		subid=SubInSub[subid].next;
	}
	return Nids;
}
