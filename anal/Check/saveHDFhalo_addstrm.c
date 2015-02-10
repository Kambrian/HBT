#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "hdf_util.h"

#define h0 0.73
#define LUNIT ((1/h0)) //kpc, for A4
// #define LUNIT ((1/h0)*1000) //kpc, for B4
#define STRMCAT_DIR  "/gpfs/data/jvbq85/HBT/data/AqA4/strmcatL"
#define NsnapStrm 127

struct PList
{
  HBTInt np;
  HBTInt *PIndex;
};
CATALOGUE Cat;
SUBCATALOGUE SubCat, StrmCat;
SRCCATALOGUE SrcCat;
extern void collect_particles(struct PList * p);
extern HBTInt *prepare_ind2sub(SUBCATALOGUE *A);
extern void dump_particles_hdf(char *outfile, HBTInt *PIndex, HBTInt np, HBTInt *ID2Sub, HBTInt *ID2Halo,  HBTInt subid);
extern void dump_subhalos_hdf(char *outfile, HBTInt grpid);
int main(int argc, char** argv)
{
	char outfile[1024];

	
	HBTInt Nsnap=MaxSnap-1;
	HBTInt grpid=0,subid;
	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_sub_table(Nsnap, &SubCat, SUBCAT_DIR);
	load_sub_catalogue(NsnapStrm,&StrmCat, STRMCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&StrmCat,-2);	
	free_PIDHash();
	
	subid=SubCat.GrpOffset_Sub[grpid];

	sprintf(outfile,"%s/anal/allhalo.hdf5",SUBCAT_DIR);
	struct PList p;
	collect_particles(&p);
	
	HBTInt *ID2Strm=prepare_ind2sub(&StrmCat);
	hid_t file_id = H5Fopen(outfile, H5F_ACC_RDWR, H5P_DEFAULT);
	herr_t status=dump_haloid(p.PIndex, p.np, ID2Strm, file_id, "/StrmID");
	if(status>=0)
	{
	  int itmp=StrmCat.GrpLen_Sub[subid];
	  H5LTset_attribute_int(file_id, "/StrmID", "NSubInGrp", &itmp,1);
	}
    status = H5Fclose (file_id);	
	myfree(ID2Strm);
	
  free_sub_table(&SubCat);
  erase_sub_catalogue(&StrmCat);
  free_particle_data();
	return 0;
}

void collect_particles(struct PList * p)
{
  LINKLIST ll;
  make_linklist(&ll, NP_DM, 50, Pdat.Pos, GetArrPos, 0);
  p->np=SubCat.SubLen[0];
  p->PIndex=linklist_search_sphere(&ll, 500/LUNIT, SubCat.Property[0].CoM, &p->np);
  printf("Mv=%g\n", p->np*header.mass[1]/h0);
  free_linklist(&ll);
}

herr_t dump_haloid(HBTInt *PIndex, HBTInt np, HBTInt *ID2Halo, hid_t file_id, char *dsetname)
{//write haloid into hdf file
  int i;
  herr_t status;
  hsize_t dims;
  if(ID2Halo!=NULL)
  {
	int *sid=mymalloc(sizeof(int)*np);
	for(i=0;i<np;i++)
	  sid[i]=ID2Halo[PIndex[i]];
	dims=np;
	status = H5LTmake_dataset(file_id, dsetname, 1, &dims, H5T_NATIVE_INT, sid); 
	myfree(sid);
	return status;
  }
  else
	return -1;
}
	
void dump_subhalos_hdf(char *outfile, HBTInt grpid)
{
    hid_t    file_id,group_id;
    herr_t      status;
    size_t     i, j,nread=0;
    hsize_t dims[2];
    
    //if(try_readfile(outfile))
    file_id = H5Fcreate (outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file

    HBTInt np, subid, cenid;
    float (* pos)[3],(* vel)[3];
    int * m;
    np=SubCat.GrpLen_Sub[grpid];
    cenid=SubCat.GrpOffset_Sub[grpid];
    pos=mymalloc(sizeof(float)*3*np);
    vel=mymalloc(sizeof(float)*3*np);
    m=mymalloc(sizeof(int)*np);
    for(i=0;i<np;i++)
    {
      subid=cenid+i;
      m[i]=SubCat.SubLen[subid];
      for(j=0;j<3;j++)
      {
	pos[i][j]=(SubCat.Property[subid].CoM[j]-SubCat.Property[cenid].CoM[j])*LUNIT;
	vel[i][j]=SubCat.Property[subid].VCoM[j]-SubCat.Property[cenid].VCoM[j];
      }
    }
    dims[0]=np;
    dims[1]=3;
    status = H5LTmake_dataset(file_id,"/x",2,dims,H5T_NATIVE_FLOAT,pos);//relative to mainsub CoM
    status = H5LTmake_dataset(file_id,"/v",2,dims,H5T_NATIVE_FLOAT,vel);
    status = H5LTmake_dataset(file_id,"/m",1,dims,H5T_NATIVE_INT,m);//number of particles
    /* close file */
    status = H5Fclose (file_id);
    printf("Pos[1]: (%g,%g,%g)\n",pos[1][0],pos[1][1],pos[1][2]);
    myfree(pos);
    myfree(vel);
    myfree(m);
}

HBTInt *prepare_ind2sub(SUBCATALOGUE *A)
{
  HBTInt i,subid,pid;
  
  HBTInt *ID2Sub=mymalloc(sizeof(HBTInt)*NP_DM);
  
  #pragma omp parallel 
  {
	#pragma omp for
	for(i=0;i<NP_DM;i++)//initialization
	{
	  ID2Sub[i]=-1;/*return -1 if the PID does not belong to a Halo,
	  i.e,we consider the backgroud as a halo with haloid=-1; 
	  note that this only make sense when we try to find host for bound structures */
	}
	#pragma omp for private(i,pid,subid)
	for(subid=0;subid<A->Nsubs;subid++)
	{
	  for(i=0;i<A->SubLen[subid];i++)
	  {
		pid=A->PSubArr[subid][i];//Pindex ranges [0,NP_DM);
		ID2Sub[pid]=subid;//haloIDs begins from id=0
	  }
	}
  }
  return ID2Sub;
}
