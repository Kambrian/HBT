/* *
 * load and dump halo virial size
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
struct HALO
{
  float cen[3];
  int Np;
  int Mvir[3];
  float Rvir[3];
  int flagbadvir[3];
};

struct HALO * load_halo_virial_size(int Ngroups,int Nsnap)
{
	char buf[1024];
	FILE *fp;
	int Nvir[3],i,j;
	struct HALO * h=mymalloc(sizeof(struct HALO)*Ngroups);
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,10*4L,SEEK_CUR);
		fread(h[i].cen,sizeof(float),3,fp);
		fread(&h[i].Np,sizeof(int),1,fp);
		fread(h[i].Mvir,sizeof(int),3,fp);
		fread(h[i].Rvir,sizeof(float),3,fp);
		fread(h[i].flagbadvir,sizeof(int),3,fp);
		fseek(fp,1*4L,SEEK_CUR);
	}
	fclose(fp);
	return h;
}

int main()
{
	double partmass;
	CATALOGUE Cat;
	int Nsnap=MaxSnap-1;

	logfile=stdout;
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	#ifdef OMEGA0
	partmass=PMass;
	#else
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];
	#endif		
	struct HALO *halo=load_halo_virial_size(Cat.Ngroups,Nsnap);
	if(halo[0].Np!=Cat.Len[0]||halo[1].Np!=Cat.Len[1])
	{
	  printf("error: fof len does not match %d=%d, %d=%d\n", (int)halo[0].Np, (int)Cat.Len[0], (int)halo[1].Np, (int)Cat.Len[1]);
	  exit(1);
	}
		
	FILE *fp;
	char buf[1024];
	sprintf(buf,"%s/anal/halos.%d", SUBCAT_DIR, (int)Nsnap);
	myfopen(fp,buf,"w");
	fprintf(fp,"#HaloID,X,Y,Z,NumPartFoF,MvirTH,MvirC200,MvirB200,RvirTH,RvirC200,RvirB200,flagBadTH,flagBadC200,flagBadB200\n");
	int i;
	for(i=0;i<Cat.Ngroups;i++)
	  fprintf(fp,"%d,%g,%g,%g,"HBTIFMT",%g,%g,%g,%g,%g,%g,%d,%d,%d\n",i,halo[i].cen[0],halo[i].cen[1],halo[i].cen[2],
		  halo[i].Np, halo[i].Mvir[0]*partmass, halo[i].Mvir[1]*partmass, halo[i].Mvir[2]*partmass,
		  halo[i].Rvir[0], halo[i].Rvir[1], halo[i].Rvir[2],
		  halo[i].flagbadvir[0],halo[i].flagbadvir[1],halo[i].flagbadvir[2]
 		);
	fclose(fp);
	return 0;
}
