//this file produces compact binary output for halo and subhalo mass data
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
//~ #include "gas_vars.h"
//~ #include "gas_proto.h"

#define SUBFIND_DIR "/home/jvbq85/data/HBT/data/AqE5W/subfind"
#define OUTDIR "/home/jvbq85/data/HBT/data/AqE5W/subcat/anal/subfind"

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#define load_halo_virial_size load_subfind_halo_size
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

extern void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);

int main(int argc,char **argv)
{
	//~ CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	//~ GASSUBCAT GSubCat;
	FILE *fp;
	char buf[1024];
	int Nsnap=70;
	int i,j,pid;
	double partmass;
	float mass,(*Mvir)[3],(*Rvir)[3];
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	#ifdef SUBFIND_DIR
	sprintf(outputdir,"%s/massfun",OUTDIR);	
	#else
	sprintf(outputdir,"%s/anal/massfun",SUBCAT_DIR);
	#endif
	mkdir(outputdir,0755);
	
	logfile=stdout;
	if(argc!=2)
	{
	printf("usage:%s [Snap]\n",argv[0]);
	exit(1);
	}
	Nsnap=atoi(argv[1]);

	sprintf(buf,"%s/submass_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	//~ load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	//~ load_gassubcat(Nsnap,&GSubCat,GASCAT_DIR);
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];
	printf("dm %g,baryon %g, fraction %g,Nsplitter %lld, Nquasi %lld\n",
			partmass,header.mass[0],header.mass[0]/header.mass[1],(long long)SubCat.Nsplitter,(long long)SubCat.NQuasi);			
	
	fwrite(&SubCat.Nsubs,sizeof(int),1,fp);			
	for(i=0;i<SubCat.Nsubs;i++)
	{
		mass=SubCat.SubLen[i]*partmass;//+GSubCat.SubLen[i]*header.mass[0];
		fwrite(&mass,sizeof(float),1,fp);
	}
	fwrite(&SubCat.Nsubs,sizeof(int),1,fp);		
	fclose(fp);
	
	sprintf(buf,"%s/subcom_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	fwrite(&SubCat.Nsubs,sizeof(int),1,fp);		
	for(i=0;i<SubCat.Nsubs;i++)
	{
		float x;
		for(j=0;j<3;j++)
		{
			x=SubCat.Property[i].CoM[j];
			fwrite(&x,sizeof(float),1,fp);
		}
	}
	fwrite(&SubCat.Nsubs,sizeof(int),1,fp);		
	fclose(fp);
			
	sprintf(buf,"%s/cid_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
	fwrite(SubCat.GrpOffset_Sub,sizeof(int),SubCat.Ngroups,fp);
	fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
	fclose(fp);
	
	Mvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,(float)partmass,SubCat.Ngroups,Nsnap);
	sprintf(buf,"%s/grpsizeVIR_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
	for(i=0;i<SubCat.Ngroups;i++)
	{
	fwrite(&Mvir[i][0],sizeof(float),1,fp);
	fwrite(&Rvir[i][0],sizeof(float),1,fp);
	}
	fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
	fclose(fp);
	sprintf(buf,"%s/grpsizeC200_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
	for(i=0;i<SubCat.Ngroups;i++)
	{
	fwrite(&Mvir[i][1],sizeof(float),1,fp);
	fwrite(&Rvir[i][1],sizeof(float),1,fp);
	}
	fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
	fclose(fp);
	sprintf(buf,"%s/grpsizeB200_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
	for(i=0;i<SubCat.Ngroups;i++)
	{
	fwrite(&Mvir[i][2],sizeof(float),1,fp);
	fwrite(&Rvir[i][2],sizeof(float),1,fp);
	}
	fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
	fclose(fp);
	myfree(Mvir);
	myfree(Rvir);
	
	return 0;
}

#ifndef SUBFIND_DIR
void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
{
	char buf[1024];
	FILE *fp;
	int Nvir[3],i,j;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,14*4L,SEEK_CUR);
		fread(Nvir,sizeof(int),3,fp);
		for(j=0;j<3;j++)
		Mvir[i][j]=Nvir[j]*partmass;
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}
#endif
