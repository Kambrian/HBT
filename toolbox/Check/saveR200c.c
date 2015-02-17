/*load virial radius and save, for majormerger project*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define CEN_COM 0
#define CEN_MBD 1
#define CEN_MPT 2  //when using this, need to define HALO_PARA in the parameter file to make the tree-code thread-safe
#define CEN_TYPE CEN_COM

HBTReal  (*Mvir)[3], (*Rvir)[3];

void init_SubVir(HBTInt Nsnap);
void free_SubVir();

int main(int argc, char** argv)
{
	HBTInt Nsnap;
	
	logfile=stdout;//redirect BT routines' log info to standard output
	
	FILE *fptab;
	char buf[1024];
	sprintf(buf,"%s/anal",SUBCAT_DIR);
	mkdir(buf,0755);

	sprintf(buf,"%s/subsize",buf);
	myfopen(fptab,buf,"w");
	fprintf(fptab,"#Snapshot M0[1e10Msun/h] M1[1e10Msun/h] R0[Mpc/h] R1[Mpc/h]\n");
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
	  init_SubVir(Nsnap);
	  fprintf(fptab,"%d %g %g %g %g\n",(int)Nsnap,Mvir[0][1],Mvir[1][1],Rvir[0][1],Rvir[1][1]);fflush(fptab);
	  printf("%d %g %g %g %g\n",(int)Nsnap,Mvir[0][1],Mvir[1][1],Rvir[0][1],Rvir[1][1]);fflush(stdout);
	  free_SubVir();
	}
	fclose(fptab);
	
	return 0;
}
void init_SubVir(HBTInt Nsnap)
{
  char buf[1024];
  FILE *fp;
  HBTInt Nsubs,dummy;
   
  #if CEN_TYPE==CEN_COM
  sprintf(buf,"%s/profile/SubMvirRvir_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
  #elif CEN_TYPE==CEN_MBD
  sprintf(buf,"%s/profile/SubMvirRvir_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
  #elif CEN_TYPE==CEN_MPT
  sprintf(buf,"%s/profile/SubMvirRvir_"HBTIFMT".MPT",SUBCAT_DIR,Nsnap);
  #endif
  myfopen(fp,buf,"r");
  fread(&Nsubs,sizeof(HBTInt),1,fp);
  printf("Nsnap=%d,Nsubs=%d\n",(int)Nsnap,(int)Nsubs);fflush(stdout);
//   if(SubCat.Nsubs!=Nsubs) 
//   {
//     printf("error loading %s:\n size not consistent with existing catalogue "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,SubCat.Nsubs);
//     exit(1);
//   }
  Mvir=mymalloc(sizeof(HBTReal)*Nsubs*3);
  Rvir=mymalloc(sizeof(HBTReal)*Nsubs*3);
  fread(Mvir,sizeof(HBTReal),3*Nsubs,fp);
  fread(Rvir,sizeof(HBTReal),3*Nsubs,fp);
  fread(&dummy,sizeof(HBTInt),1,fp);
  if(dummy!=Nsubs) 
  {
    printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
    exit(1);
  }
  fclose(fp);
}

void free_SubVir()
{
  myfree(Mvir);
  myfree(Rvir);
}
