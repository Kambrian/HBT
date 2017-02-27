#include <stdio.h>
#include <stdlib.h>

typedef long long HBTInt;
typedef float HBTReal;

typedef struct 
{//units: km/s, kpc/h, 1e10Msun/h
  HBTInt 	subhaloId;//subhaloIndex at the current snapshot, 0~Nsub-1
  HBTInt	numberOfBoundParticles; //can be 0 for the central of unbound FoF halos
  HBTInt	mostBoundID; //-1 if numberOfBoundParticles==0
  HBTInt	descendantSubId;//-1 if none
  HBTInt        progenitorSubId;//-1 if none
  HBTInt        hostFoFId; //0~Nfof-1 for normal subhalos, can be -1 if the sub is not hosted by any fof (small subhalos in the field)
  HBTInt        numberOfSubsInFoF;
  HBTInt        firstSubInFoF;
// HBTInt		firstProgenitor;
// HBTInt		nextProgenitor;
// HBTInt		firstHaloInFOFgroup;
// HBTInt		nextHaloInFOFgroup;
  int 		isCentral;
  float 	comovingPosition[3];
  float 	physicalVelocity[3];
  float 	physicalSpecificAngularMomentum[3]; //<R_physical x V_physical>
  float 	physicalVelDisp;
  float 	physicalVmax;
  float 	comovingHalfMassRadius;
  float 	m_Mean200; //value for the host halo
  float 	m200;
  float 	m_TopHat;
//   int 		snapNum;
}TREE;

int main(int argc, char** argv)
{
int Nsnap=10;	
FILE *fp;
char buf[1024];
sprintf(buf,"/mnt/ddnfs/jxhan/6620/subcat/anal/subtree_%03d",Nsnap);
fp=fopen(buf,"r");
HBTInt snapid, Ngroups, Nsubs;
HBTReal scale_factor;
fread(&snapid, sizeof(HBTInt), 1, fp);
fread(&scale_factor, sizeof(HBTReal), 1, fp);
fread(&Ngroups, sizeof(HBTInt), 1, fp);
fread(&Nsubs, sizeof(HBTInt), 1, fp);
TREE *subtree=malloc(sizeof(TREE)*Nsubs);
fread(subtree, sizeof(TREE), Nsubs, fp);
fclose(fp);

int i=2;
printf("Ngroups=%ld, Nsubs=%ld, scale_factor=%g\n", Ngroups, Nsubs, scale_factor);
printf("%ld, %ld, %ld, %ld, %d, %g\n", subtree[i].subhaloId,  subtree[i].numberOfBoundParticles, subtree[i].mostBoundID, subtree[i].hostFoFId, subtree[i].isCentral, subtree[i].m200);

free(subtree);
return 0;
}