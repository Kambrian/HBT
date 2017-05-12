#include <stdio.h>
#include <stdlib.h>

typedef long long HBTInt;

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

int main(int argc, char** argv)
{
    int Nsnap=10;
    FILE *fp;
    char buf[1024];
    sprintf(buf,"/mnt/ddnfs/jxhan/6620/subcat/anal/subtree_%03d",Nsnap);
    fp=fopen(buf,"r");
    int snapid, Ngroups, Nsubs;
    float scale_factor;
    fread(&snapid, sizeof(int), 1, fp);
    fread(&scale_factor, sizeof(float), 1, fp);
    fread(&Ngroups, sizeof(int), 1, fp);
    fread(&Nsubs, sizeof(int), 1, fp);
    TREE *subtree=malloc(sizeof(TREE)*Nsubs);
    fread(subtree, sizeof(TREE), Nsubs, fp);
    fclose(fp);

    /*read in the DescendantList*/
    int *descendantID=malloc(sizeof(int)*Nsubs);
    sprintf(buf,"/mnt/ddnfs/jxhan/6620/subcat/anal/DescendantList/snap_%03d",Nsnap);
    fp=fopen(buf,"r");
    fread(descendantID, sizeof(int), Nsubs, fp);
    fclose(fp);
    /*fix the descendantID in the tree*/
    int i;
    for(i=0;i<Nsubs;i++)
      subtree[i].descendantSubId=descendantID[i];
    free(descendantID);
    
    
    i=2;
    printf("snapid=%d, Ngroups=%d, Nsubs=%d, scale_factor=%g\n", snapid, Ngroups, Nsubs, scale_factor);
    printf("%d, %d, %ld, %d, %d, %g\n", subtree[i].subhaloId,  subtree[i].numberOfBoundParticles, subtree[i].mostBoundID, subtree[i].hostFoFId, subtree[i].isCentral, subtree[i].m200);

    free(subtree);
        
    return 0;
}