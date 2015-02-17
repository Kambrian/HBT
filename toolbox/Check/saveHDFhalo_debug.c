#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "hdf_util.h"

#define SUBFIND_DIR GRPCAT_DIR
#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#define INPUT_DIR SUBFIND_DIR
#define EXTNAME ".subfind"
#else
#define EXTNAME ""
#define INPUT_DIR SUBCAT_DIR
#endif

#define h0 0.73
#ifdef CONVERT_LENGTH_MPC_KPC //already converted to kpc
  #define LUNIT ((1/h0)) //kpc, for A4
#else
  #define LUNIT ((1/h0)*1000) //kpc, for B4, A2, B2,....
#endif

struct PList
{
  HBTInt np;
  HBTInt *PIndex;
};
CATALOGUE Cat;
SUBCATALOGUE SubCat;
SRCCATALOGUE SrcCat;
extern void collect_particles(struct PList * p);
extern void collect_subhaloes(struct PList * p);
extern void update_CoMVCoM(SUBCATALOGUE *SubCat);
extern HBTInt *prepare_ind2sub(SUBCATALOGUE *A);
extern void dump_particles_hdf(char *outfile, HBTInt *PIndex, HBTInt np, HBTInt *ID2Sub, HBTInt *ID2Halo,  HBTInt subid);
extern void dump_subhalos_hdf(char *outfile, struct PList *p);
static int comp_IDatInt(const void *a, const void *b)//used to sort PID in ascending order; 
{
  if((*(IDatInt *)a) > (*(IDatInt *)b))
    return +1;

  if((*(IDatInt *)a) < (*(IDatInt *)b))
    return -1;

  return 0;
}
#define PrintWatch(x) {printf(x":\n");	printf("%g,%g,%g\n", Pdat.Pos[indwatch][0],Pdat.Pos[indwatch][1],Pdat.Pos[indwatch][2]);printf("%g,%g,%g\n", Pdat.Vel[indwatch][0],Pdat.Vel[indwatch][1],Pdat.Vel[indwatch][2]);}	
HBTInt indwatch;
int PIDFound=-1;
int main(int argc, char** argv)
{
	char outfile[1024];

	struct PList p;
	HBTInt Nsnap=MaxSnap-1;
	HBTInt grpid=0,subid;
	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
/**/
	long long idwatch=741164718857480175LL;
  #ifdef HBTPID_RANKSTYLE
  IDatInt *PIDs;  
  PIDs=load_PIDs_Sorted();
  indwatch=(IDatInt *)bsearch(&(idwatch),PIDs,NP_DM,sizeof(IDatInt),comp_IDatInt)-PIDs;
  printf("%lld, %lld\n", idwatch, PIDs[indwatch]);
  myfree(PIDs);
  #else
  indwatch=idwatch;
  #endif
  
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
// 	load_src_catalogue(Nsnap, &SrcCat, SUBCAT_DIR);
// 	load_sub_table(Nsnap, &SubCat, SUBCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,INPUT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fresh_ID2Index(&Cat,-1); 	
	fresh_ID2Index(&SubCat,-2);	
// 	fresh_ID2Index(&SrcCat,-3);
	indwatch=lookup_ID2Ind(indwatch);
	PrintWatch("Fresh")
	free_PIDHash();
	
#ifdef SUBFIND_DIR
	update_CoMVCoM(&SubCat);
#else
	printf("CoM=%g,%g,%g\n", SubCat.Property[0].CoM[0], SubCat.Property[0].CoM[1], SubCat.Property[0].CoM[2]);
	printf("Cen=%g,%g,%g\n", Pdat.Pos[SubCat.PSubArr[0][0]][0], Pdat.Pos[SubCat.PSubArr[0][0]][1], Pdat.Pos[SubCat.PSubArr[0][0]][2]);
	printf("VCoM=%g,%g,%g\n", SubCat.Property[0].VCoM[0], SubCat.Property[0].VCoM[1], SubCat.Property[0].VCoM[2]);
	printf("VCen=%g,%g,%g\n", Pdat.Vel[SubCat.PSubArr[0][0]][0], Pdat.Vel[SubCat.PSubArr[0][0]][1], Pdat.Vel[SubCat.PSubArr[0][0]][2]);
#endif
	PrintWatch("CoMVCoM")
	subid=SubCat.GrpOffset_Sub[grpid];
	
// 	sprintf(outfile,"%s/anal/halo.hdf5",SUBCAT_DIR);
// 	dump_particles_hdf(outfile,Cat.PIDorIndex+Cat.Offset[grpid],Cat.Len[grpid],subid);
	
	sprintf(outfile,"%s/anal/allpart"EXTNAME".hdf5.debug",SUBCAT_DIR);
	collect_particles(&p);
	PrintWatch("Collected")
	int i;
	for(i=0;i<p.np;i++)
	{
	  if(p.PIndex[i]==indwatch)
	  {
		PIDFound=i;
		break;
	  }
	}
	printf("PIDFound=%d\n", PIDFound);
	prepare_ind2halo(&Cat);
	PrintWatch("Ind2halo")
	HBTInt *ID2Sub=prepare_ind2sub(&SubCat);
	PrintWatch("Ind2sub")
	dump_particles_hdf(outfile,p.PIndex, p.np, ID2Sub, Cat.ID2Halo, subid);
	PrintWatch("Dumped")
	myfree(ID2Sub);
	myfree(p.PIndex);
	
	// 	load_sub_table(Nsnap,&SubCat, SUBCAT_DIR);
	sprintf(outfile,"%s/anal/sublist"EXTNAME".hdf5",SUBCAT_DIR);
	collect_subhaloes(&p);
	dump_subhalos_hdf(outfile, &p);
	myfree(p.PIndex);
// 	free_sub_table(&SubCat);
	
// 	sprintf(outfile,"%s/anal/subhalo.hdf5",SUBCAT_DIR);
// 	dump_particles_hdf(outfile,SubCat.PSubArr[subid],SubCat.SubLen[subid], NULL, NULL, subid);
	
// 	sprintf(outfile,"%s/anal/src.hdf5",SUBCAT_DIR);
// 	dump_particles_hdf(outfile,SrcCat.PSubArr[subid],SrcCat.SubLen[subid], NULL, NULL, subid);	

	
	free_catalogue(&Cat);
// free_sub_table(&SubCat);
// erase_src_catalogue(&SrcCat);
	erase_sub_catalogue(&SubCat);
	free_particle_data();
	return 0;
}

void collect_particles(struct PList * p)
{
  LINKLIST ll;
  make_linklist(&ll, NP_DM, 50, Pdat.Pos, GetArrPos, 0);
	PrintWatch("linklisted")  
  p->np=Cat.Len[0];
  p->PIndex=linklist_search_sphere(&ll, 500/LUNIT, SubCat.Property[0].CoM, &p->np);
  	PrintWatch("Searched")
  printf("Mv=%g\n", p->np*header.mass[1]/h0);
  free_linklist(&ll);
}

void collect_subhaloes(struct PList * p)
{
  LINKLIST ll;
  make_linklist(&ll, SubCat.Nsubs, 50, SubCat.Property, GetSubCoM, 0);
  p->np=SubCat.GrpLen_Sub[0];
  p->PIndex=linklist_search_sphere(&ll, 500/LUNIT, SubCat.Property[0].CoM, &p->np);
  printf("N=%d, %d\n", (int)p->np, (int)SubCat.GrpLen_Sub[0]);
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
	
void dump_particles_hdf(char *outfile, HBTInt *PIndex, HBTInt np, HBTInt *ID2Sub, HBTInt *ID2Halo,  HBTInt subid)
{//set ID2Sub=NULL and ID2Halo=NULL if you don't want to record the ids
    hid_t    file_id,group_id;
    herr_t      status;
    size_t     i, j,nread=0;
    hsize_t dims[2];

    //if(try_readfile(outfile))
    file_id = H5Fcreate (outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file
//     group_id = HDFcreate_group(file_id, "/PartType0");
//     status = H5Gclose(group_id);

	float sqa;
	#ifdef VEL_INPUT_PHYSICAL
	  sqa=1.0;
	#else
	  sqa = sqrt(header.time);
	#endif
    float (* pos)[3],(* vel)[3];
    pos=mymalloc(sizeof(float)*3*np);
    vel=mymalloc(sizeof(float)*3*np);
    for(i=0;i<np;i++)
    {
      for(j=0;j<3;j++)
      {
	pos[i][j]=(Pdat.Pos[PIndex[i]][j]-SubCat.Property[subid].CoM[j])*LUNIT;
	vel[i][j]=Pdat.Vel[PIndex[i]][j]*sqa-SubCat.Property[subid].VCoM[j];
      }
    }
    PrintWatch("Dumping...")
	if(PIDFound>=0)
	{
	printf("Pos=%g,%g,%g\n", pos[PIDFound][0], pos[PIDFound][1], pos[PIDFound][2]);
	printf("Vel=%g,%g,%g\n", vel[PIDFound][0], vel[PIDFound][1], vel[PIDFound][2]);	
	}
    dims[0]=1;
	float pmass=header.mass[1]/h0;
    status = H5LTmake_dataset(file_id,"/PartMass", 1, dims, H5T_NATIVE_FLOAT, &pmass); 
    dims[0]=np;
    dims[1]=3;
    status = H5LTmake_dataset(file_id,"/x",2,dims,H5T_NATIVE_FLOAT,pos);
    status = H5LTmake_dataset(file_id,"/v",2,dims,H5T_NATIVE_FLOAT,vel);
	float x0[3],v0[3];
	for(j=0;j<3;j++)
	{
	  x0[j]=SubCat.Property[subid].CoM[j];
	  v0[j]=SubCat.Property[subid].VCoM[j];
	}
	H5LTset_attribute_float(file_id, "/x", "x0", x0, 3);
	H5LTset_attribute_float(file_id, "/v", "v0", v0, 3);
	
	status=dump_haloid(PIndex, np, ID2Sub, file_id, "/SubID");
	if(status>=0)
	{
	  int itmp=SubCat.GrpLen_Sub[0];
	  H5LTset_attribute_int(file_id, "/SubID", "NSubInGrp", &itmp,1);
	}
	status=dump_haloid(PIndex, np, ID2Halo, file_id, "/HaloID");
		
    /* close file */
    status = H5Fclose (file_id);
    printf("Pos[1]: (%g,%g,%g)\n",pos[1][0],pos[1][1],pos[1][2]);
    myfree(pos);
    myfree(vel);
}

void dump_subhalos_hdf(char *outfile, struct PList *p)
{
    hid_t    file_id,group_id;
    herr_t      status;
    size_t     i, j,nread=0;
    hsize_t dims[2];
    
    //if(try_readfile(outfile))
    file_id = H5Fcreate (outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //always create a new file or overwite existing file

    HBTInt np, subid, cenid;
    float (* pos)[3],(* vel)[3];
    float * m;
	int *ID, *haloID;
    np=p->np;
    cenid=0;
    pos=mymalloc(sizeof(float)*3*np);
    vel=mymalloc(sizeof(float)*3*np);
    m=mymalloc(sizeof(float)*np);
	ID=mymalloc(sizeof(int)*np);
	haloID=mymalloc(sizeof(int)*np);
    for(i=0;i<np;i++)
    {
      subid=p->PIndex[i];
	  ID[i]=subid;
	  haloID[i]=SubCat.HaloChains[subid].HostID;
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
	float x0[3],v0[3];
	for(j=0;j<3;j++)
	{
	  x0[j]=SubCat.Property[subid].CoM[j];
	  v0[j]=SubCat.Property[subid].VCoM[j];
	}
	H5LTset_attribute_float(file_id, "/x", "x0", x0, 3);
	H5LTset_attribute_float(file_id, "/v", "v0", v0, 3);
    status = H5LTmake_dataset(file_id,"/PartMass",1,dims,H5T_NATIVE_FLOAT,m);//number of particles
	status = H5LTmake_dataset(file_id,"/SubID",1,dims,H5T_NATIVE_INT,ID);
	status = H5LTmake_dataset(file_id,"/HaloID",1,dims,H5T_NATIVE_INT,haloID);
    /* close file */
    status = H5Fclose (file_id);
    printf("Pos[1]: (%g,%g,%g)\n",pos[1][0],pos[1][1],pos[1][2]);
    myfree(pos);
    myfree(vel);
    myfree(m);
	myfree(ID);
	myfree(haloID);
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

void update_CoMVCoM(SUBCATALOGUE *SubCat)
{//update CoM and VCoM, using 10% most bound particles
#define CoreFracMbd 0.10 
#define CoreLenMinMbd 10
  HBTInt subid;
  HBTInt i,j,np;
  HBTReal sqa;
  #ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
  #else
	sqa = sqrt(header.time);
  #endif
  for(subid=0;subid<SubCat->Nsubs;subid++)
  {
	np=(int)(SubCat->SubLen[subid]*CoreFracMbd);
	if(np<CoreLenMinMbd) np=CoreLenMinMbd;
	if(np>SubCat->SubLen[subid]) np=SubCat->SubLen[subid];
	center_of_mass(SubCat->Property[subid].CoM, SubCat->PSubArr[subid], np, Pdat.Pos);
	for(j=0;j<3;j++) SubCat->Property[subid].VCoM[j]=0;
	for(i=0;i<np;i++)
	{
	  HBTInt pid=SubCat->PSubArr[subid][i];
	  for(j=0;j<3;j++)
		SubCat->Property[subid].VCoM[j]+=Pdat.Vel[pid][j];
	}
	for(j=0;j<3;j++)
	  SubCat->Property[subid].VCoM[j]*=sqa/np;
  }
  printf("CoM=%g,%g,%g\n", SubCat->Property[0].CoM[0], SubCat->Property[0].CoM[1], SubCat->Property[0].CoM[2]);
  printf("Cen=%g,%g,%g\n", Pdat.Pos[SubCat->PSubArr[0][0]][0], Pdat.Pos[SubCat->PSubArr[0][0]][1], Pdat.Pos[SubCat->PSubArr[0][0]][2]);
  printf("VCoM=%g,%g,%g\n", SubCat->Property[0].VCoM[0], SubCat->Property[0].VCoM[1], SubCat->Property[0].VCoM[2]);
  printf("VCen=%g,%g,%g\n", Pdat.Vel[SubCat->PSubArr[0][0]][0]*sqa, Pdat.Vel[SubCat->PSubArr[0][0]][1]*sqa, Pdat.Vel[SubCat->PSubArr[0][0]][2]*sqa);
}
