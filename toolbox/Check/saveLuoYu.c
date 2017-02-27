/* to save subhalo data for Subhalos Go Notts project */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define CEN_COM 0
#define CEN_MBD 1
#define CEN_MPT 2  //when using this, need to define HALO_PARA in the parameter file to make the tree-code thread-safe

#define CEN_TYPE CEN_COM

#define VIR_TOPHAT 0
#define VIR_C200 1
#define VIR_B200 2

SUBCATALOGUE SubCat;
HBTReal *Vmax, *Rmax, *Rhalf;
float (*Mvir)[3], (*Rvir)[3];
HBTReal (*SubMvir)[3], (*SubRvir)[3];

typedef struct
{   //units: km/s, kpc/h, 1e10Msun/h
    HBTInt 	mostBoundID; //-1 if numberOfBoundParticles==0
    int 	subhaloId;//subhaloIndex at the current snapshot, 0~Nsub-1
    int 	numberOfBoundParticles; //can be 0 for the central of unbound FoF halos
    int 	descendantSubId;//-1 if none
    int 	progenitorSubId;//-1 if none
    int 	hostFoFId; //0~Nfof-1 for normal subhalos, can be -1 if the sub is not hosted by any fof (small subhalos in the field)
    int 	numberOfSubsInFoF;
    int 	firstSubInFoF;
    int 	isCentral;
    float 	comovingPosition[3];
    float 	physicalVelocity[3];
    float 	physicalSpecificAngularMomentum[3]; //<R_physical x V_physical>
    float 	physicalVelDisp;
    float 	physicalVmax;
    float 	comovingHalfMassRadius;
    float 	m_Mean200; //value for the host halo
    float 	m200;
    float 	m_TopHat;
} TREE;

void init_HaloSizes(HBTInt Nsnap);
void free_HaloSizes();
int main(int argc, char** argv)
{
    HBTInt Nsnap;
    HBTInt subid;

    logfile=stdout;//redirect BT routines' log info to standard output

// 	Nsnap=MaxSnap-1;
// 	if(argc!=2)
// 	{printf("usage: %s [Nsnap], otherwise Nsnap="HBTIFMT"\n",argv[0],Nsnap);fflush(stdout);}
// 	else
    int Nsnap0=atoi(argv[1]);
    int NsnapD=atoi(argv[2]);
    for(Nsnap=Nsnap0; Nsnap<MaxSnap; Nsnap+=NsnapD)
    {
        load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
        printf("Nsnap=%d, Nsubs=%ld\n",(int)Nsnap, SubCat.Nsubs);
        HBTInt *pro2dest;
        if(Nsnap<MaxSnap-1)
        {
            HBTInt nsubs;
            load_pro2dest(Nsnap,&pro2dest,&nsubs,SUBCAT_DIR);
            assert(nsubs>=SubCat.Nsubs);
        }
        else
        {
            pro2dest=mymalloc(sizeof(HBTInt)*SubCat.Nsubs);
            #pragma omp parallel for
            for(subid=0; subid<SubCat.Nsubs; subid++) pro2dest[subid]=-1;
        }
        HBTInt Npro, Nsplitter, *sp2pro;
        load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);
        load_particle_header(Nsnap,SNAPSHOT_DIR);
        float sqa = sqrt(header.time);
        init_HaloSizes(Nsnap);

        printf("filling tree...\n");
        fflush(stdout);
// #define FixProId(x) ((x)<Npro?(x):(sp2pro[x])) //or assign -1?
#define FixProId(x) ((x)<Npro?(x):-1)
#define copy_xyz(x,y) {x[0]=y[0];x[1]=y[1];x[2]=y[2];}
        TREE *subtree=malloc(sizeof(TREE)*SubCat.Nsubs);
        #pragma omp parallel for
        for(subid=0; subid<SubCat.Nsubs; subid++)
        {
            TREE *tree=subtree+subid;
            tree->subhaloId=subid;
            tree->mostBoundID=(SubCat.SubLen[subid]>0)?SubCat.PSubArr[subid][0]:-1;
            tree->descendantSubId=pro2dest[subid];
            tree->progenitorSubId=FixProId(SubCat.HaloChains[subid].ProSubID);
            HBTInt grpid=SubCat.HaloChains[subid].HostID;
            tree->hostFoFId=grpid;
            if(grpid<0)
            {
                tree->numberOfSubsInFoF=1;
                tree->firstSubInFoF=subid;
                tree->isCentral=1;
                tree->m_Mean200=SubMvir[subid][VIR_B200];
                tree->m200=SubMvir[subid][VIR_C200];
                tree->m_TopHat=SubMvir[subid][VIR_TOPHAT];
            }
            else
            {
                tree->numberOfSubsInFoF=SubCat.GrpLen_Sub[grpid];
                tree->firstSubInFoF=SubCat.GrpOffset_Sub[grpid];
                tree->isCentral=(SubCat.SubRank[subid]==0);
                tree->m_Mean200=Mvir[grpid][VIR_B200];
                tree->m200=Mvir[grpid][VIR_C200];
                tree->m_TopHat=Mvir[grpid][VIR_TOPHAT];
            }
            tree->numberOfSubsInFoF=(grpid<0?1:SubCat.GrpLen_Sub[grpid]);
            tree->firstSubInFoF=(grpid<0?subid:SubCat.GrpOffset_Sub[grpid]);
            tree->numberOfBoundParticles=SubCat.SubLen[subid];
            copy_xyz(tree->comovingPosition, SubCat.Property[subid].CoM);
            copy_xyz(tree->physicalVelocity, SubCat.Property[subid].VCoM);
            copy_xyz(tree->physicalSpecificAngularMomentum, SubCat.Property[subid].AM);
            tree->physicalVelDisp=sqrt(SubCat.Property[subid].Kin*2);//to improve: subtract (vcen-vmean)^2
            tree->physicalVmax=Vmax[subid]/sqa;
            tree->comovingHalfMassRadius=Rhalf[subid];
        }
        printf("saving ...\n");
        fflush(stdout);
        FILE *fp;
        char buf[1024];
        sprintf(buf,"%s/anal/subtree_%03d",SUBCAT_DIR,(int)Nsnap);
        myfopen(fp,buf,"w");
	{
	  int snapid=snaplist[Nsnap], ngroups=SubCat.Ngroups, nsubs=SubCat.Nsubs;
	  float scale=header.time;
	  fwrite(&snapid, sizeof(int), 1, fp);
	  fwrite(&scale, sizeof(float), 1, fp);
	  fwrite(&ngroups, sizeof(int), 1, fp);
	  fwrite(&nsubs, sizeof(int), 1, fp);
	  fwrite(subtree, sizeof(TREE), nsubs, fp);
	}
        fclose(fp);

        printf("saved.\n");
        fflush(stdout);
        erase_sub_catalogue(&SubCat);
        free_HaloSizes();
        free_pro2dest(pro2dest);
        free_sp2pro(sp2pro, Npro, Nsplitter);
        free(subtree);
    }

    return 0;
}


int load_RmaxVmax(HBTReal *rmax,HBTReal *vmax, HBTReal *rhalf, HBTInt Nsnap)
{
    char buf[1024];
    FILE *fp;
    HBTInt Nsubs,dummy;
    //~ HBTReal *rsig, *r2sig, *r3sig, *rpoisson; /* declare this as input var if you want them!!! */

    if(!rmax||!vmax||!rhalf)
    {
// 		printf("error: allocate rmax , vmax and rhalf first \n");
// 		exit(1);
    }
#if CEN_TYPE==CEN_COM
    sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
#elif CEN_TYPE==CEN_MPT
    sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
#else
    assert(0==1);
#endif
    myfopen(fp,buf,"r");
    fread(&Nsubs,sizeof(HBTInt),1,fp);
    fread(rmax,sizeof(HBTReal),Nsubs,fp);
    fread(vmax,sizeof(HBTReal),Nsubs,fp);
    fread(rhalf,sizeof(HBTReal),Nsubs,fp);
    fseek(fp,sizeof(HBTReal)*Nsubs*4,SEEK_CUR);
    //~ fread(rsig,sizeof(HBTReal),Nsubs,fp);
    //~ fread(r2sig,sizeof(HBTReal),Nsubs,fp);
    //~ fread(r3sig,sizeof(HBTReal),Nsubs,fp);
    //~ fread(rpoisson,sizeof(HBTReal),Nsubs,fp);
    fread(&dummy,sizeof(HBTInt),1,fp);
    if(dummy!=Nsubs)
    {
        printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
        exit(1);
    }
    fclose(fp);
    return Nsubs;
}


void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
{
    typedef struct
    {
        char pre[14*4];
        int Nvir[3];
        float Rvir[3];
        char tail[4*4];
    } halosize_t;
    halosize_t *halosize=mymalloc(sizeof(halosize_t)*Ngroups);

    char buf[1024];
    FILE *fp;
    int i;
    sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
    myfopen(fp,buf,"r");
    fread(halosize, sizeof(halosize_t), Ngroups, fp);
    fclose(fp);

    #pragma omp parallel for default(shared)
    for(i=0; i<Ngroups; i++)
    {
        int j;
        for(j=0; j<3; j++)
            Mvir[i][j]=halosize[i].Nvir[j]*partmass;
        for(j=0; j<3; j++)
            Rvir[i][j]=halosize[i].Rvir[j];
    }
    myfree(halosize);

}

void load_sub_virial_size(HBTReal Mvir[][3],HBTReal Rvir[][3], HBTInt Nsnap)
{
    FILE *fp;
    char buf[1024];
#if CEN_TYPE==CEN_COM
    sprintf(buf,"%s/profile/SubMvirRvir_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
#elif CEN_TYPE==CEN_MBD
    sprintf(buf,"%s/profile/SubMvirRvir_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
#elif CEN_TYPE==CEN_MPT
    sprintf(buf,"%s/profile/SubMvirRvir_"HBTIFMT".MPT",SUBCAT_DIR,Nsnap);
#endif
    HBTInt Nsubs, Nsubs2;
    myfopen(fp,buf,"r");
    fread(&Nsubs,sizeof(HBTInt),1,fp);
    fread(Mvir, sizeof(HBTReal), 3*Nsubs, fp);
    fread(Rvir, sizeof(HBTReal), 3*Nsubs, fp);
    fread(&Nsubs2,sizeof(HBTInt),1,fp);
    assert(Nsubs2==Nsubs);
    fclose(fp);
}

void init_HaloSizes(HBTInt Nsnap)
{
    Rmax=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
    Vmax=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
    Rhalf=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
    load_RmaxVmax(Rmax,Vmax,Rhalf,Nsnap);
    SubMvir=mymalloc(sizeof(HBTReal)*SubCat.Nsubs*3);
    SubRvir=mymalloc(sizeof(HBTReal)*SubCat.Nsubs*3);
    load_sub_virial_size(SubMvir, SubRvir, Nsnap);
    Mvir=malloc(sizeof(float)*3*SubCat.Ngroups);
    Rvir=malloc(sizeof(float)*3*SubCat.Ngroups);
    load_halo_virial_size(Mvir, Rvir, header.mass[1], SubCat.Ngroups, Nsnap);
}
void free_HaloSizes()
{
    myfree(Rmax);
    myfree(Vmax);
    myfree(Rhalf);
    myfree(SubMvir);
    myfree(SubRvir);
    myfree(Mvir);
    myfree(Rvir);
}
