/*load virial radius and save, for majormerger project*/
//for host halos
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

void init_HostVir(HBTInt Nsnap);
void free_HostVir();


float (*Rvir)[3];
float (*Mvir)[3];
CATALOGUE Cat;
SUBCATALOGUE SubCat;
int main(int argc, char** argv)
{

	int Nsnap=99;

	logfile=stdout;
			
	FILE *fptab;
	char buf[1024];
	sprintf(buf,"%s/anal",SUBCAT_DIR);
	mkdir(buf,0755);

	sprintf(buf,"%s/hostsize",buf);
	myfopen(fptab,buf,"w");
	fprintf(fptab,"#Snapshot M0[1e10Msun/h] M1[1e10Msun/h] R0[Mpc/h] R1[Mpc/h]\n");
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
	  init_HostVir(Nsnap);
	  HBTInt host0=SubCat.HaloChains[0].HostID;
	  HBTInt host1=SubCat.HaloChains[1].HostID;
	  fprintf(fptab,"%d %g %g %g %g\n",(int)Nsnap,Mvir[host0][1],Mvir[host1][1],Rvir[host0][1],Rvir[host1][1]);fflush(fptab);
	  printf("%d %g %g %g %g\n",(int)Nsnap,Mvir[host0][1],Mvir[host1][1],Rvir[host0][1],Rvir[host1][1]);fflush(stdout);
	  free_HostVir();
	}
	fclose(fptab);
	
	return 0;
}
void free_HostVir()
{
  myfree(Mvir);
  myfree(Rvir);
  free_catalogue(&Cat);
  free_sub_table(&SubCat);
}
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
// 			Mvir[i][j]=Nvir[j];
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}
void init_HostVir(HBTInt Nsnap)
{
  double partmass;
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];	
	Mvir=mymalloc(sizeof(int)*3*Cat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*Cat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,(float)partmass,Cat.Ngroups,Nsnap);
}

