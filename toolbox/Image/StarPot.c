//To calculate Potential energy for each star particle, for Wenting
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "hdf_util.h"

HBTxyz *StarPos;
HBTInt Nstar;
#define GRPID 0
extern HBTInt load_star_pos(char *starfile);
extern void save_star_pot(char *outfile, float *pot);
int main(int argc,char **argv)
{
  CATALOGUE Cat;
  HBTInt Nsnap=MaxSnap-1;
  
  char starfile[1024];
  char outfile[1024];
  logfile=stdout;
  
  if(argc<2)
  {
    printf("Usage: %s [StarFile] <outfile>\n Outfile=StarFile if not explicitly specified.\n Now exit.\n", argv[0]);
    exit(1);
  }
  strcpy(starfile,argv[1]);
  if(argc<3)
    strcpy(outfile,argv[1]);
  else
    strcpy(outfile,argv[2]);
  
  //prepare particles
  load_star_pos(starfile);
  load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
  load_particle_data_bypart(Nsnap,SNAPSHOT_DIR,FLAG_LOAD_POS|FLAG_LOAD_ID);
  fill_PIDHash();
  fresh_ID2Index(&Cat,FRSH_GRPCAT); 
  free_PIDHash();
  
  HBTInt i, Np, *PIndex;
  Np=Cat.Len[GRPID];
  PIndex=Cat.PIDorIndex+Cat.Offset[GRPID];
  tree_tree_allocate(TREE_ALLOC_FACTOR*Np,Np);
  maketree(Np,PIndex,Pdat.Pos);
  float *pot=mymalloc(sizeof(float)*Nstar);
  #pragma omp parallel for
  for(i=0;i<Nstar;i++)
  {
    pot[i]=tree_treeevaluate_potential(StarPos[i],PIndex,Pdat.Pos);
    pot[i]*=header.mass[1]*G/header.time;
  }
  tree_tree_free();
 
  save_star_pot(outfile, pot);
  myfree(pot);
  free_catalogue(&Cat);
  free_particle_data();
  
  return 0;
}
void save_star_pot(char *outfile, float *pot)
{
    hid_t    file_id,group_id;
    herr_t      status;
    hsize_t n=Nstar;
    char dsetname[1024]="/Stars/Potential";
    
    if(try_readfile(outfile)) 
      file_id = H5Fopen(outfile, H5F_ACC_RDWR, H5P_DEFAULT);
    else
    {
      file_id = H5Fcreate (outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      group_id = HDFcreate_group(file_id, "/Stars");
      status = H5Gclose(group_id);
    }
    
    status = H5LTmake_dataset(file_id,dsetname,1,&n,H5T_NATIVE_FLOAT,pot);
    H5LTset_attribute_string(file_id,dsetname,"Unit","(km/s)^2, -GM/R_physical");
    status = H5Fclose(file_id); 
}

HBTInt load_star_pos(char *starfile)
{
  FloatMat A;
  strcpy(A.name,"/Stars/Coordinates");
  load_hdfmatrixF(starfile,&A,1);
  if(A.dim!=2||A.size[1]!=3) {printf("error: unexpected dimension\n");exit(1);};
  printf("Dataset:%ldx%ld\n",A.size[0],A.size[1]);
  Nstar=A.size[0];
  StarPos=mymalloc(sizeof(HBTxyz)*Nstar);
  HBTInt i,pid;
  for(pid=0;pid<Nstar;pid++)
    for(i=0;i<3;i++)
      StarPos[pid][i]=A.x[3*pid+i];
  printf("%d stars loaded\n",(int)Nstar);  fflush(stdout);  
  return Nstar;  
}

