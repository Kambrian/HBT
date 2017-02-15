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
    float D;
    float V;
    float MHost;
    float RHost;
    float VmaxHost;
    float Mbound;
    float Vmax;
    float padding;//to work around library bug
    HBTInt HostId;
    HBTInt SubId;
    HBTInt CenId;
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
void makell_sub();
void collect_submass(int grpid,SatList_t *sats);
void freell_sub();
extern float *export_Vmax();
int main(int argc,char **argv)
{
	HBTInt Nsnap,i,j;
	
	FILE *fp;
	char buf[1024];
	
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
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];
	Mvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,(float)partmass,SubCat.Ngroups,Nsnap);
	
	makell_sub();
        HBTInt grpid;
        SatList_t Satellites;
        init_satlist(&Satellites, SubCat.Nsubs);
    #pragma omp parallel for
    for(grpid=0; grpid<SubCat.Ngroups; grpid++)
    {
        if(SubCat.GrpLen_Sub[grpid]==0) continue;
        SatList_t satlist;
        collect_submass(grpid, &satlist);
        #pragma omp critical (insert_satellite)
        append_satlist(&Satellites, &satlist);
        free_satlist(&satlist);
    }
	
	freell_sub();
	
	save(&Satellites, Nsnap);
	myfree(Mvir);
	myfree(Rvir);
	
	return 0;
}

#define NDIV 200
static int hoc[NDIV][NDIV][NDIV],*ll;
static float range[3][2], step[3];
void makell_sub()
{
	int i,j,grid[3],np;
	
	#define POS(i,j) SubCat.Property[i].CoM[j]
	np=SubCat.Nsubs;
	printf("creating linked list..\n");
	ll=mymalloc(sizeof(int)*np);
	/*determining enclosing cube*/
	for(i=0;i<3;i++)
	{
		range[i][0]=0.;
		range[i][1]=BOXSIZE;
	}
	for(j=0;j<3;j++)
		step[j]=BOXSIZE/NDIV;
	//~ /*determining enclosing cube*/
	//~ for(i=0;i<3;i++)
		//~ for(j=0;j<2;j++)
			//~ range[i][j]=POS(0,i);
	//~ for(i=1;i<np;i++)
		//~ for(j=0;j<3;j++)
		//~ {
			//~ if(POS(i,j)<range[j][0])
				//~ range[j][0]=POS(i,j);
			//~ else if(POS(i,j)>range[j][1])
				//~ range[j][1]=POS(i,j);
		//~ }
	//~ for(j=0;j<3;j++)
		//~ step[j]=(range[j][1]-range[j][0])/NDIV;
	
	/*initialize hoc*/
	int *phoc=&(hoc[0][0][0]);
	for(i=0;i<NDIV*NDIV*NDIV;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((POS(i,j)-range[j][0])/step[j]);
			if(grid[j]<0) 
				grid[j]+=NDIV;
			else if(grid[j]>=NDIV)
				grid[j]-=NDIV;
		}
		ll[i]=hoc[grid[0]][grid[1]][grid[2]];
		hoc[grid[0]][grid[1]][grid[2]]=i; /*use hoc(floor(xsat)+1) as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
	#undef POS
}
void freell_sub()
{
	myfree(ll);
}
void collect_submass(int grpid, SatList_t *satlist)
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
	subbox_grid[i][0]=floor((cen[i]-rscale-range[i][0])/step[i]);
	//~ if(subbox_grid[i][0]<0)subbox_grid[i][0]=0;
	subbox_grid[i][1]=floor((cen[i]+rscale-range[i][0])/step[i]);
	//~ if(subbox_grid[i][1]>=NDIV)subbox_grid[i][1]=NDIV-1;
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	
#define GRID_IN_BOX(x) ((x)<0?((x)+NDIV):((x)>=NDIV?((x)-NDIV):(x)))
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=hoc[GRID_IN_BOX(i)][GRID_IN_BOX(j)][GRID_IN_BOX(k)];
				while(pid>=0)
				{
					dr=distance(SubCat.Property[pid].CoM,cen);
					if(dr<rmax)
					{
                                            if(satlist->len==satlist->capacity) expand_satlist(satlist);
                                            Satellite_t *sat=satlist->list+satlist->len;
                                            sat->D=dr;
                                            sat->V=raw_distance(SubCat.Property[pid].VCoM, SubCat.Property[cenid].VCoM);
                                            sat->MHost=Mvir[grpid][VirType];
                                            sat->RHost=Rvir[grpid][VirType];
                                            sat->VmaxHost=SubVmax[cenid];
                                            sat->Mbound=SubCat.SubLen[pid]*partmass;
                                            sat->Vmax=SubVmax[pid];
                                            sat->HostId=grpid;
                                            sat->SubId=pid;
                                            sat->CenId=cenid;
                                            satlist->len++;
					}
					pid=ll[pid];
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
    InsertMember(V, H5T_NATIVE_FLOAT);
    InsertMember(MHost, H5T_NATIVE_FLOAT);
    InsertMember(RHost, H5T_NATIVE_FLOAT);
    InsertMember(VmaxHost, H5T_NATIVE_FLOAT);
    InsertMember(Mbound, H5T_NATIVE_FLOAT);
    InsertMember(Vmax, H5T_NATIVE_FLOAT);
    InsertMember(HostId, H5T_HBTInt);
    InsertMember(SubId, H5T_HBTInt);
    InsertMember(CenId, H5T_HBTInt);
#undef InsertMember

//     H5Tclose(H5T_HBTxyz);
    return H5T_dtypeInMem;
}

void save(SatList_t *Satellites, int isnap)
{
    char outdir[1024]=OUTDIR;
    mkdir(outdir,0755);
    char buf[1024];
    sprintf(buf,"%s/masslist_%03d.%d.hdf5",outdir,isnap,VirType);
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
