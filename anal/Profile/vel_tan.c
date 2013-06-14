//tangential velocity distribution, in terms of Kt=vt^2/vc^2
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define GRPMAX 1

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
int load_halo_concentration(float *halocon,int Nsnap)
{
	int grpid,Ngroups,Ngroups2,halostatus;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_param_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	fread(&Ngroups,sizeof(int),1,fp);
	for(grpid=0;grpid<Ngroups;grpid++)
	{
		fread(&halostatus,sizeof(int),1,fp);
		fseek(fp,4*sizeof(float),SEEK_CUR);
		//~ fread(halostatus+grpid,sizeof(int),1,fp);
		fread(halocon+grpid,sizeof(float),1,fp);
		if(halostatus!=0) halocon[grpid]=-1; //set concentration to -1 if fitting not successful
		fseek(fp,sizeof(float)*2,SEEK_CUR);
	}	
	fread(&Ngroups2,sizeof(int),1,fp);
	if(Ngroups2!=Ngroups)
	{
		printf("error:Ngroups=%d,%d do not match when loading \n %s\n" 
			"probably file corruption or different file format\n",
			Ngroups,Ngroups2,buf);
		exit(-1);
	}
	fclose(fp);
	return Ngroups;
}
float NFWm(float x)
{ //x=r/rs;
	return log(1+x)-x/(1+x); //M(<x), normalized by Ms
}
int main(int argc, char** argv)
{
	SUBCATALOGUE SubCat;
	float (*Mvir)[3],(*Rvir)[3],*halocon;
	
	float *Kt,*R;
	int Nsnap=99;
	int i,grpid,subid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	Mvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	halocon=mymalloc(sizeof(float)*SubCat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,header.mass[1],SubCat.Ngroups,Nsnap);
	load_halo_concentration(halocon,Nsnap);
	Kt=mymalloc(sizeof(float)*SubCat.GrpOffset_Sub[GRPMAX]);
	R=mymalloc(sizeof(float)*SubCat.GrpOffset_Sub[GRPMAX]);
	float *cen,*vcen,Pos[3],Vel[3],r,v,vc2,sint2;
	for(grpid=0;grpid<GRPMAX;grpid++)
	{
		cen=SubCat.Property[grpid].CoM;
		vcen=SubCat.Property[grpid].VCoM;
		for(subid=SubCat.GrpOffset_Sub[grpid];subid<SubCat.GrpOffset_Sub[grpid+1];subid++)
		{
			r=distance(cen,SubCat.Property[subid].CoM);
			R[subid]=r/Rvir[grpid][0];
			v=distance(vcen,SubCat.Property[subid].VCoM);
			float mass;//mass within r
			if(R[subid]>1.0)
			mass=Mvir[grpid][0];
			else
			mass=Mvir[grpid][0]*NFWm(R[subid]*halocon[grpid])/NFWm(halocon[grpid]);
			vc2=G*mass/r/header.time;
			for(i=0;i<3;i++)
			{
				Pos[i]=SubCat.Property[subid].CoM[i]-cen[i];
				Vel[i]=SubCat.Property[subid].VCoM[i]-vcen[i];
			}
			sint2=1.0-pow(f_prod(Pos,Vel,3)/r/v,2);
			Kt[subid]=v*v/vc2*sint2;
		}
	}
	
	FILE *fp;
	char buf[1024];
	sprintf(buf,"%s/anal/Kt_%d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"w");
	for(subid=0;subid<SubCat.GrpOffset_Sub[GRPMAX];subid++)
	fprintf(fp,"%g,%g,%d\n",Kt[subid],R[subid],SubCat.SubLen[subid]);
	fclose(fp);
	
	free_sub_table(&SubCat);
	return 0;
}
