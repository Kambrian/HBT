//outputs final massfunction data for plots
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "hdf_util.h"

#define RMAX 1.5	 //statistics done in RMIN*rvi<r<RMAX*rvir


#define SUBFIND_DIR SNAPSHOT_DIR
#define OUTDIR "/gpfs/data/jvbq85/HBT/data/Millennium2/subcat/analysis/subfind/"

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#define load_halo_virial_size load_subfind_halo_size
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

typedef struct
{
  HBTInt MostBoundID;
  int IsCentral;
  int Nbound;
  int NboundPeak;
  float Vmax;
  float VmaxPeak;
  int BranchID;
  int HostID;
} DParticle_t;

typedef struct 
{
    float D;
    float MHost;
    float RHost;
    float VmaxHost;
    float Mbound;
    float Vmax;
    float MboundPeak;
    float VmaxPeak;
    int HostId;
    int HostFoFId;
    int IsFoFCentral;
    HBTInt MostBoundID;
} Satellite_t;

typedef struct
{
    Satellite_t *list;
    HBTInt len;
    HBTInt capacity;
} SatList_t;
void init_satlist(SatList_t *sats,HBTInt cap)
{
    sats->list=mymalloc(sizeof(Satellite_t)*cap);
    sats->capacity=cap;
    sats->len=0;
}
void expand_satlist(SatList_t *sats)
{
    sats->capacity*=2;
    sats->list=realloc(sats->list, sizeof(Satellite_t)*sats->capacity);
}
void append_satlist(SatList_t *Sats, SatList_t *sat)
{
    if(Sats->len+sat->len>Sats->capacity) expand_satlist(Sats);
    memcpy(Sats->list+Sats->len, sat->list, sizeof(Satellite_t)*sat->len);
    Sats->len+=sat->len;
}
void free_satlist(SatList_t *sats)
{
    free(sats->list);
    sats->capacity=0;
    sats->len=0;
}
void save(SatList_t *Satellites, int isnap);

double partmass;
float (*Mvir)[3],(*Rvir)[3];
SUBCATALOGUE SubCat;
float *SubVmax;
int VirType;

void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);
void collect_submass(int grpid,SatList_t *sats, LINKLIST *ll, DParticle_t *Branches);
extern float *export_Vmax();
DParticle_t * load_dparticles(HBTInt *np);
struct BranchPosData_t
{
  HBTxyz *Pos;
  HBTInt *ParticleIndex;
}BranchPos;
HBTReal *GetBranchPos(HBTInt pid, void *data)
{
  return BranchPos.Pos[BranchPos.ParticleIndex[pid]];
}
int main(int argc,char **argv)
{
	HBTInt Nsnap,i,j;
	
	logfile=stdout;
	if(argc!=3)
	{
	printf("usage:%s [Snap] [VirialType]\n",argv[0]);
	exit(1);
	}
	Nsnap=atoi(argv[1]);
	VirType=atoi(argv[2]);

	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
        SubVmax=export_Vmax();
	load_particle_data_bypart(Nsnap, SNAPSHOT_DIR, FLAG_LOAD_ID|FLAG_LOAD_POS);
	partmass=header.mass[1];
	Mvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,(float)partmass,SubCat.Ngroups,Nsnap);
	HBTInt NumBranches;
	DParticle_t *Branches=load_dparticles(&NumBranches);
	BranchPos.Pos=Pdat.Pos;
	BranchPos.ParticleIndex=malloc(sizeof(HBTInt)*NumBranches);
	for(i=0;i<NumBranches;i++)
	  BranchPos.ParticleIndex[i]=Branches[i].MostBoundID;
	fill_PIDHash();
	fresh_ID2Index(BranchPos.ParticleIndex, NumBranches);
	free_PIDHash();	
	
	for(i=0;i<NumBranches;i++)
	  if(Branches[i].IsCentral)
	  {
	    HBTInt subid=SubCat.GrpOffset_Sub[Branches[i].HostID];
	    memcpy(SubCat.Property[subid].CoM, GetBranchPos(i, NULL), sizeof(HBTxyz));
	  }
	
	LINKLIST ll;
	make_linklist(&ll, NumBranches, 128, NULL, GetBranchPos, 1);
	
        HBTInt grpid;
        SatList_t Satellites;
        init_satlist(&Satellites, SubCat.Nsubs);
    #pragma omp parallel for
    for(grpid=0; grpid<SubCat.Ngroups; grpid++)
    {
        if(SubCat.GrpLen_Sub[grpid]==0) continue;
        SatList_t satlist;
        collect_submass(grpid, &satlist, &ll, Branches);
        #pragma omp critical (insert_satellite)
        append_satlist(&Satellites, &satlist);
        free_satlist(&satlist);
    }
	
	save(&Satellites, Nsnap);
	free_linklist(&ll);
	free_particle_data();
	free(Branches);
	free(BranchPos.ParticleIndex);
	myfree(Mvir);
	myfree(Rvir);
	
	return 0;
}

void collect_submass( int grpid, SatList_t *satlist, LINKLIST *ll, DParticle_t *Branches)
{
	float rmax=RMAX*Rvir[grpid][VirType];
	int i,j,k,cenid,pid,subbox_grid[3][2];
	float rscale,dr;
	HBTReal *cen;
	
        init_satlist(satlist, SubCat.GrpLen_Sub[grpid]);
	cenid=SubCat.GrpOffset_Sub[grpid];
	if(SubCat.SubLen[cenid])
	{
	cen=SubCat.Property[cenid].CoM;
	rscale=rmax;
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rscale-ll->range[i][0])/ll->step[i]);
	subbox_grid[i][1]=floor((cen[i]+rscale-ll->range[i][0])/ll->step[i]);
	#ifndef PERIODIC_BDR //do not fix if periodic, since the search sphere is allowed to overflow the box in periodic case.
	subbox_grid[i][0]=linklist_fix_gridid(subbox_grid[i][0],ll);
	subbox_grid[i][1]=linklist_fix_gridid(subbox_grid[i][1],ll);
	#endif
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc_safe(ll,i,j,k); //in case the grid-id is out of box, in the periodic case
				while(pid>=0)
				{
					dr=distance(GetBranchPos(pid, NULL),cen);
					if(dr<rmax)
					{
                                            if(satlist->len==satlist->capacity) expand_satlist(satlist);
                                            Satellite_t *sat=satlist->list+satlist->len;
                                            sat->D=dr;
                                            sat->MHost=Mvir[grpid][VirType];
                                            sat->RHost=Rvir[grpid][VirType];
                                            sat->VmaxHost=SubVmax[cenid];
                                            sat->Mbound=Branches[pid].Nbound*partmass;
                                            sat->Vmax=Branches[pid].Vmax;
					    sat->MboundPeak=Branches[pid].NboundPeak*partmass;
					    sat->VmaxPeak=Branches[pid].VmaxPeak;
                                            sat->HostId=grpid;
					    sat->HostFoFId=Branches[pid].HostID;
					    sat->IsFoFCentral=Branches[pid].IsCentral;
					    sat->MostBoundID=Branches[pid].MostBoundID;
                                            satlist->len++;
					}
					pid=ll->list[pid];
				}
			}
	}
}

#ifndef SUBFIND_DIR
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
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}
#endif


hid_t BuildHDFSatellite()
{
    hid_t H5T_dtypeInMem=H5Tcreate(H5T_COMPOUND, sizeof (Satellite_t));
//     hsize_t dims[2]= {3,3};
//     hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);

#define InsertMember(x,t) H5Tinsert(H5T_dtypeInMem, #x, HOFFSET(Satellite_t, x), t)
    InsertMember(D, H5T_NATIVE_FLOAT);
    InsertMember(MHost, H5T_NATIVE_FLOAT);
    InsertMember(RHost, H5T_NATIVE_FLOAT);
    InsertMember(VmaxHost, H5T_NATIVE_FLOAT);
    InsertMember(Mbound, H5T_NATIVE_FLOAT);
    InsertMember(Vmax, H5T_NATIVE_FLOAT);
    InsertMember(MboundPeak, H5T_NATIVE_FLOAT);
    InsertMember(VmaxPeak, H5T_NATIVE_FLOAT);
    InsertMember(HostId, H5T_NATIVE_INT);
    InsertMember(HostFoFId, H5T_NATIVE_INT);
    InsertMember(IsFoFCentral, H5T_NATIVE_INT);
    InsertMember(MostBoundID, H5T_HBTInt);
#undef InsertMember

//     H5Tclose(H5T_HBTxyz);
    return H5T_dtypeInMem;
}

void save(SatList_t *Satellites, int isnap)
{
    char outdir[1024]=OUTDIR;
    mkdir(outdir,0755);
    char buf[1024];
    sprintf(buf,"%s/masslistDTree_%03d.%d.hdf5",outdir,isnap,VirType);
    hid_t file=H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dtype=BuildHDFSatellite();
    hid_t dtype2=H5Tcopy(dtype);
    H5Tpack(dtype2);
    hsize_t dims[]= {Satellites->len};
    writeHDFmatrix(file, Satellites->list, "Satellites", 1, dims, dtype, dtype2);
//     herr_t status = H5LTmake_dataset(file,"Satellites",1,dims,dtype,Satellites->list);
    H5Tclose(dtype2);
    H5Tclose(dtype);
    H5Fclose(file);
}

hid_t BuildHDFParticle()
{
    hid_t H5T_dtypeInMem=H5Tcreate(H5T_COMPOUND, sizeof (DParticle_t));
//     hsize_t dims[2]= {3,3};
//     hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);

#define InsertMember(x,t) H5Tinsert(H5T_dtypeInMem, #x, HOFFSET(DParticle_t, x), t)
    InsertMember(MostBoundID, H5T_HBTInt);
    InsertMember(IsCentral, H5T_NATIVE_INT);
    InsertMember(Nbound, H5T_NATIVE_INT);
    InsertMember(NboundPeak, H5T_NATIVE_INT);
    InsertMember(Vmax, H5T_NATIVE_FLOAT);
    InsertMember(VmaxPeak, H5T_NATIVE_FLOAT);
    InsertMember(BranchID, H5T_NATIVE_INT);
    InsertMember(HostID, H5T_NATIVE_INT);
#undef InsertMember

//     H5Tclose(H5T_HBTxyz);
    return H5T_dtypeInMem;
}

DParticle_t * load_dparticles(HBTInt *Np)
{
    hid_t dtype=BuildHDFParticle();
    char indir[1024]="/gpfs/data/jvbq85/HBT/data/Millennium2/subcat/analysis/DTrees";
    int ifile;
    HBTInt np=0, offset[128];
    for(ifile=0;ifile<128;ifile++)
    {
      offset[ifile]=np;
      char buf[1024];  
      sprintf(buf,"%s/DTreeMostBoundParticles.%d.hdf5",indir,ifile);
      hid_t file=H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t dset=H5Dopen2(file, "Branches", H5P_DEFAULT);
      hsize_t dims[1];
      GetDatasetDims(dset, dims);
      H5Dclose(dset);
      H5Fclose(file);
      np+=dims[0];
    }
    DParticle_t *DParticles=malloc(sizeof(DParticle_t)*np);
    *Np=np;
    
    for(ifile=0;ifile<128;ifile++)
    {
      char buf[1024];  
      sprintf(buf,"%s/DTreeMostBoundParticles.%d.hdf5",indir,ifile);
      hid_t file=H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t dset=H5Dopen2(file, "Branches", H5P_DEFAULT);
      H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, DParticles+offset[ifile]);
      H5Dclose(dset);
      H5Fclose(file);
    }
    H5Tclose(dtype);
    
    return DParticles;
}