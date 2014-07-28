#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

HBTReal *Rmax, *Vmax;
SUBCATALOGUE SubCat;
CATALOGUE Cat;
extern HBTReal group_VelDisp(HBTInt grpid);
extern void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);
float (*Mvir)[3], (*Rvir)[3];
int main(int argc,char **argv)
{

	
	HBTInt Nsnap;
	
	FILE *fp;
	char buf[1024],outputdir[1024];
	sprintf(outputdir,"/data/A4700r2d1/kambrain/6113/post/anal/");
	mkdir(outputdir,0755);		
	logfile=stdout;
	
	if(argc!=2)
	{
		printf(" %s [Nsnap]\n",argv[0]);
		exit(1);
	}
	Nsnap=atoi(argv[1]);
	HBTInt Nsnap0=Nsnap;
	
	sprintf(buf, "%s/Halo_%03d",outputdir,(int)Nsnap);
	myfopen(fp,buf,"w");

	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_particle_data_bypart(Nsnap,SNAPSHOT_DIR,FLAG_LOAD_VEL|FLAG_LOAD_ID);
	fill_PIDHash();
	fresh_ID2Index(&Cat,FRSH_GRPCAT); 
	free_PIDHash();
	
	float partmass=header.mass[1];	
	Mvir=mymalloc(sizeof(float)*3*Cat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*Cat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,partmass,Cat.Ngroups,Nsnap);
	HBTInt grpid;
	HBTReal *sigma;
	sigma=mymalloc(sizeof(HBTReal)*Cat.Ngroups);
#pragma omp parallel for
	for(grpid=0;grpid<Cat.Ngroups;grpid++)
	  sigma[grpid]=group_VelDisp(grpid);
	
	fprintf(fp,"HaloID,NpFoF,Nbound,Mvir,Rvir,x,y,z,Vx,Vy,Vz,SigmaV,SigmaVbnd\n");
	for(grpid=0;grpid<Cat.Ngroups;grpid++)
	{
	  HBTInt subid=SubCat.GrpOffset_Sub[grpid];
	  fprintf(fp,"%d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
		  (int)grpid,(int)Cat.Len[grpid],(int)SubCat.SubLen[subid],Mvir[grpid][0],Rvir[grpid][0],
	   SubCat.Property[subid].CoM[0],SubCat.Property[subid].CoM[1],SubCat.Property[subid].CoM[2],
	   SubCat.Property[subid].VCoM[0],SubCat.Property[subid].VCoM[1],SubCat.Property[subid].VCoM[2],
	   sigma[grpid],sqrt(2*SubCat.Property[subid].Kin));
	}
	fclose(fp);
	
	return 0;
}
HBTReal group_VelDisp(HBTInt grpid)
{
  HBTInt i,j,*PIndex,pid;
  HBTReal V2[3]={0.},V[3]={0.},sigma=0.;
  PIndex=Cat.PIDorIndex+Cat.Offset[grpid];
  for(i=0;i<Cat.Len[grpid];i++)
  {
    pid=PIndex[i];
    for(j=0;j<3;j++)
    {
      V2[j]+=Pdat.Vel[pid][j]*Pdat.Vel[pid][j];
      V[j]+=Pdat.Vel[pid][j];
    }
  }
  for(j=0;j<3;j++)
  {
    V2[j]/=Cat.Len[grpid];
    V[j]/=Cat.Len[grpid];
    sigma+=V2[j]-V[j]*V[j];
  }
  
  return sqrt(sigma);
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