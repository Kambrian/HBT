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
// 		fread(Mvir,sizeof(int),3,fp);
		fread(Nvir,sizeof(int),3,fp);
		for(j=0;j<3;j++)
		  Mvir[i][j]=Nvir[j]*partmass;
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}

int main()
{
	double partmass;
	float (*Mvir)[3],(*Rvir)[3];
	CATALOGUE Cat;
	int Nsnap=59;

	logfile=stdout;
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	#ifdef OMEGA0
	partmass=PMass;
	#else
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];
	#endif		
	Mvir=mymalloc(sizeof(float)*3*Cat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*Cat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,(float)partmass,Cat.Ngroups,Nsnap);
	printf("Mvir[0]=%f,Rvir[0]=%f\n",Mvir[0][0],Rvir[0][0]);
	printf("Mc200[0]=%f,Rc200[0]=%f\n",Mvir[0][1],Rvir[0][1]);
	printf("Mb200[0]=%f,Rb200[0]=%f\n",Mvir[0][2],Rvir[0][2]);
	
	return 0;
}
