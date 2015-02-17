/* *
 * load halo virial size
 * Mvir[Ngroups][3]: [M_virial,M_critical200,M_mean200], in units 10^10Msun/h
 * Rvir[Ngroups][3]: [R_virial,R_critical200,R_mean200], comoving size,kpc/h
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

void load_halo_virial_size(int Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
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
		//~ Mvir[i][j]=Nvir[j]*partmass;
			Mvir[i][j]=Nvir[j];
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}

void load_halo_virial_size_old(int Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
{
	char buf[1024];
	FILE *fp;
	int Nvir[3],i,j;
	sprintf(buf,"%s/profile_old/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,14*4L,SEEK_CUR);
		fread(Nvir,sizeof(int),3,fp);
		for(j=0;j<3;j++)
		//~ Mvir[i][j]=Nvir[j]*partmass;
			Mvir[i][j]=Nvir[j];
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}

int main()
{
	double partmass;
	float (*Rvir)[3], (*Rvir_old)[3];
	int (*Mvir)[3], (*Mvir_old)[3];
	CATALOGUE Cat;
	int Nsnap=99;

	logfile=stdout;
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	#ifdef OMEGA0
	partmass=PMass;
	#else
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];
	#endif		
	Mvir=mymalloc(sizeof(int)*3*Cat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*Cat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,(float)partmass,Cat.Ngroups,Nsnap);
	Mvir_old=mymalloc(sizeof(int)*3*Cat.Ngroups);
	Rvir_old=mymalloc(sizeof(float)*3*Cat.Ngroups);
	load_halo_virial_size_old(Mvir_old,Rvir_old,(float)partmass,Cat.Ngroups,Nsnap);
	printf("Mvir[0]=%d,Rvir[0]=%f\n",Mvir[0][0],Rvir[0][0]);
	printf("Mc200[0]=%d,Rc200[0]=%f\n",Mvir[0][1],Rvir[0][1]);
	printf("Mb200[0]=%d,Rb200[0]=%f\n",Mvir[0][2],Rvir[0][2]);
	
	printf("Mvir[0]=%d,Rvir[0]=%f\n",Mvir_old[0][0],Rvir_old[0][0]);
	printf("Mc200[0]=%d,Rc200[0]=%f\n",Mvir_old[0][1],Rvir_old[0][1]);
	printf("Mb200[0]=%d,Rb200[0]=%f\n",Mvir_old[0][2],Rvir_old[0][2]);
	
	FILE *fp;
	char buf[1024];
	sprintf(buf,"%s/anal/virmass.%d", SUBCAT_DIR, (int)Nsnap);
	myfopen(fp,buf,"w");
	int i;
	for(i=0;i<Cat.Ngroups;i++)
	fprintf(fp,HBTIFMT","HBTIFMT","HBTIFMT"\n",Cat.Len[i],Mvir[i][0],Mvir_old[i][0]);
	fclose(fp);
	
	return 0;
}
