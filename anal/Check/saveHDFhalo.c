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
// #define LUNIT ((1/h0)) //kpc, for A4
#define LUNIT ((1/h0)*1000) //kpc, for B4
struct PList
{
  HBTInt np;
  HBTInt *PIndex;
};
CATALOGUE Cat;
SUBCATALOGUE SubCat;
extern void collect_particles(struct PList * p);
extern void dump_particles_hdf(char *outfile, HBTInt *PIndex, HBTInt np, HBTInt subid);
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
	
/**/	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCat,-2);	//fresh_ID2Index(&SrcCat,-3);
	free_PIDHash();
	
	subid=SubCat.GrpOffset_Sub[grpid];
		
// 	sprintf(outfile,"%s/anal/halo.hdf5",SUBCAT_DIR);
// 	dump_particles_hdf(outfile,Cat.PIDorIndex+Cat.Offset[grpid],Cat.Len[grpid],subid);
	
	sprintf(outfile,"%s/anal/allhalo.hdf5",SUBCAT_DIR);
	struct PList p;
	collect_particles(&p);
	dump_particles_hdf(outfile,p.PIndex,p.np,subid);
	
// 	sprintf(outfile,"%s/anal/subhalo.hdf5",SUBCAT_DIR);
// 	dump_particles_hdf(outfile,SubCat.PSubArr[subid],SubCat.SubLen[subid],subid);
/**
	load_sub_table(Nsnap,&SubCat, SUBCAT_DIR);
	sprintf(outfile,"%s/anal/sublist.hdf5",SUBCAT_DIR);
	dump_subhalos_hdf(outfile, grpid);
	free_sub_table(&SubCat);
*/	
// 	free_catalogue(&Cat);
// 	erase_sub_catalogue(&SubCat);
	return 0;
}

void collect_particles(struct PList * p)
{
  LINKLIST ll;
  make_linklist(&ll, NP_DM, 50, Pdat.Pos, GetArrPos, 0);
  p->np=Cat.Len[0];
  p->PIndex=linklist_search_sphere(&ll, 500/LUNIT, SubCat.Property[0].CoM, &p->np);
  printf("Mv=%g\n", p->np*header.mass[1]/h0);
  free_linklist(&ll);
}

void dump_particles_hdf(char *outfile, HBTInt *PIndex, HBTInt np, HBTInt subid)
{
    hid_t    file_id,group_id;
    herr_t      status;
    size_t     i, j,nread=0;
    hsize_t dims[2];

    //if(try_readfile(outfile))
    file_id = H5Fcreate (outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file
//     group_id = HDFcreate_group(file_id, "/PartType0");
//     status = H5Gclose(group_id);

    float (* pos)[3],(* vel)[3];
    pos=mymalloc(sizeof(float)*3*np);
    vel=mymalloc(sizeof(float)*3*np);
    for(i=0;i<np;i++)
    {
      for(j=0;j<3;j++)
      {
	pos[i][j]=(Pdat.Pos[PIndex[i]][j]-SubCat.Property[subid].CoM[j])*LUNIT;
	vel[i][j]=Pdat.Vel[PIndex[i]][j]-SubCat.Property[subid].VCoM[j];
      }
    }
    dims[0]=1;
	float pmass=header.mass[1]/h0;
    status = H5LTmake_dataset(file_id,"/PartMass", 1, dims, H5T_NATIVE_FLOAT, &pmass); 
    dims[0]=np;
    dims[1]=3;
    status = H5LTmake_dataset(file_id,"/x",2,dims,H5T_NATIVE_FLOAT,pos);
    status = H5LTmake_dataset(file_id,"/v",2,dims,H5T_NATIVE_FLOAT,vel);
    /* close file */
    status = H5Fclose (file_id);
    printf("Pos[1]: (%g,%g,%g)\n",pos[1][0],pos[1][1],pos[1][2]);
    myfree(pos);
    myfree(vel);
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