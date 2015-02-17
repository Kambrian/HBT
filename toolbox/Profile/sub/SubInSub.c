/* to generate spatial sub-in-sub hierarchy */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define CEN_COM  //using CoM as center; otherwise using MostBnd (recommended).

#define RHALF 0
#define R1SIG 1
#define R2SIG 2
#define R3SIG 3
#define RPOISSON 4

#define RSUB_TYPE R2SIG 

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
	
SUBCATALOGUE SubCat;
HBTReal *RSub;
struct Hierarchy *SubInSub;

void init_SubInSub()
{
	HBTInt subid;
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
		SubInSub[subid].pre=-1;
		SubInSub[subid].next=-1;
		SubInSub[subid].nibs=-1;
		SubInSub[subid].sub=-1;
	}
}

void insert_subinsub(HBTInt hostid,HBTInt subid)
{
	HBTInt oldsub;
	oldsub=SubInSub[hostid].sub;
	if(oldsub>=0)
	{
		SubInSub[oldsub].pre=subid;
		SubInSub[subid].next=oldsub;
	}
	SubInSub[hostid].sub=subid;	
	SubInSub[subid].nibs=hostid;
}

void find_hostsub(HBTInt subid,HBTInt idmin)
{
	HBTReal dr;
	HBTInt hostid;
	if(0==SubCat.SubLen[subid]) return;
	for(hostid=subid-1;hostid>=idmin;hostid--)
	{
	if(0==SubCat.SubLen[hostid]) continue;	
	#ifdef CEN_COM	
	dr=distance(SubCat.Property[subid].CoM,SubCat.Property[hostid].CoM);
	#else
	dr=distance(Pdat.Pos[SubCat.PSubArr[subid][0]],Pdat.Pos[SubCat.PSubArr[hostid][0]]);
	#endif
	if(dr<RSub[hostid])
	{
		insert_subinsub(hostid,subid);
		break;
	}
	}
}

void process_groups(HBTInt grpid)
{
	HBTInt subid,idmin,idmax;
	idmin=SubCat.GrpOffset_Sub[grpid];
	idmax=SubCat.GrpOffset_Sub[grpid]+SubCat.GrpLen_Sub[grpid]-1;
	for(subid=idmax;subid>=idmin;subid--)
	find_hostsub(subid,idmin);
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


void save_SubInSub(HBTInt Nsnap)
{
	char buf[1024],outputdir[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;
	
	sprintf(outputdir,"%s/profile/",SUBCAT_DIR);	
	mkdir(outputdir,0755);
	#ifdef CEN_COM
	sprintf(buf,"%s/profile/SubInSub_%03d.COM.%d",SUBCAT_DIR,Nsnap,RSUB_TYPE);
	#else
	sprintf(buf,"%s/profile/SubInSub_%03d.MBD.%d",SUBCAT_DIR,Nsnap,RSUB_TYPE);
	#endif
	myfopen(fp,buf,"w");
	fwrite(&SubCat.Nsubs,sizeof(HBTInt),1,fp);
	fwrite(SubInSub,sizeof(struct Hierarchy),SubCat.Nsubs,fp);
	fwrite(&SubCat.Nsubs,sizeof(HBTInt),1,fp);	
	fclose(fp);
}

int main()
{
	HBTInt Nsnap,grpid;
	logfile=stdout;
	
	Nsnap=MaxSnap-1;
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&SubCat,FRSH_SUBCAT);
	free_PIDHash();
	SubInSub=mymalloc(sizeof(struct Hierarchy)*SubCat.Nsubs);
	init_SubInSub();
	init_RSub(Nsnap);
	for(grpid=0;grpid<SubCat.Ngroups;grpid++)
	process_groups(grpid);
	
	//Quasi
	HBTInt subid,idmin,idmax;
	idmin=SubCat.Nsubs-SubCat.NQuasi;
	idmax=SubCat.Nsubs-1;
	for(subid=idmax;subid>=idmin;subid--)
	find_hostsub(subid,idmin);
	
	save_SubInSub(Nsnap);
	free_particle_data();
}
