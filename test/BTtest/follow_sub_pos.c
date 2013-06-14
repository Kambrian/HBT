//to follow the position of a subhalo
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
void load_halo_virial_size(int Mvir[][3],float Rvir[][3],int Ngroups,int Nsnap);
int main(int argc, char** argv)
{
	SUBCATALOGUE SubCat;
	int Nsubs,*pro2dest;
	int SnapLoad,Nsnap;
	int grpid,subid,hsubid;
	int (*Mvir)[3];
	float (*Rvir)[3],d;
	FILE *fp;
	char buf[1024];

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=3)
	{
	printf("usage: %s [SnapLoad], [fofid]\n",argv[0]);fflush(stdout);
	exit(1);
	}
	else
	{
	SnapLoad=atoi(argv[1]);
	grpid=atoi(argv[2]);
	}
	
	sprintf(buf, "%s/anal/follow/follow_%03d_%d.pos",SUBCAT_DIR,SnapLoad,grpid);
	myfopen(fp,buf,"w");

	for(Nsnap=SnapLoad;Nsnap<MaxSnap;Nsnap++)
	{
		load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
		if(Nsnap==SnapLoad)
			subid=SubCat.GrpOffset_Sub[grpid];
		Mvir=mymalloc(sizeof(int)*3*SubCat.Ngroups);
		Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);	
		load_halo_virial_size(Mvir,Rvir,SubCat.Ngroups,Nsnap);
		grpid=SubCat.HaloChains[subid].HostID;
		hsubid=SubCat.GrpOffset_Sub[grpid];
		d=distance(SubCat.Property[subid].CoM,SubCat.Property[hsubid].CoM);
		fprintf(fp,"%d,%d,%d,%f\n",Nsnap,SubCat.SubLen[subid],Mvir[grpid][0],d/Rvir[grpid][0]);fflush(fp);
		myfree(Mvir);
		myfree(Rvir);
		erase_sub_catalogue(&SubCat);
		if(Nsnap<MaxSnap-1)
		{
			load_pro2dest(Nsnap,&pro2dest,&Nsubs,SUBCAT_DIR);
			subid=pro2dest[subid];
			free_pro2dest(pro2dest);
			if(subid<0) break;	
		}
	}
	
	fclose(fp);
	return 0;
}

void load_halo_virial_size(int Mvir[][3],float Rvir[][3],int Ngroups,int Nsnap)
{
	char buf[1024];
	FILE *fp;
	int Nvir[3],i,j;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,14*4L,SEEK_CUR);
		fread(Mvir+i,sizeof(int),3,fp);
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}
