#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

typedef struct
{
  HBTInt BranchID; //id of branch, starting from 0.
  HBTInt SubID; 
  HBTInt HostID;
  HBTInt SubRank; //0: central; >0: satellite.
  HBTInt NpBnd; //number of bound particles
  HBTInt NpBndPeak; //maximum value of NBound in all previous snapshots
  HBTInt SnapNumPeak; //snapshot number at which NBoundPeak is found.
  HBTInt MstBndID; //ID of most bound particle
  HBTxyz MstBndPos; //comoving coordinate of MostBound particle
  HBTxyz MstBndVel; //physical velocity of MostBound particle
  HBTReal Vmax; //physical vmax
  HBTReal VmaxPeak; //peak of physical vmax
  HBTInt SnapNumVpeak; //snapshot number at VmaxPeak
} BranchNode;

typedef struct
{
  HBTInt BranchID; //id of branch, starting from 0.
  HBTInt SubID; 
  HBTInt HostID;
  HBTInt SubRank; //0: central; >0: satellite.
  HBTInt NpBnd; //number of bound particles
  HBTInt NpBndPeak; //maximum value of NBound in all previous snapshots
  HBTInt SnapNumPeak; //snapshot number at which NBoundPeak is found.
  HBTInt MstBndID; //ID of most bound particle
  HBTxyz MstBndPos; //comoving coordinate of MostBound particle
  HBTxyz MstBndVel; //physical velocity of MostBound particle
  HBTReal Vmax; //physical vmax
  HBTReal VmaxPeak; //peak of physical vmax
  HBTReal SnapNumVpeak; //snapshot number at VmaxPeak
} BranchNodeBug;

void copynode(BranchNodeBug *in, BranchNode *out)
{
    memcpy(out, in, sizeof(BranchNodeBug));
    out->SnapNumVpeak=in->SnapNumVpeak;
}

#define BLOCKSIZE 1024
void fixfile(int Nsnap, char *path)
{
    printf("Nsnap=%d\n", Nsnap);
  FILE *fp, *fpout;
  char buf[1024];
  sprintf(buf, "%s/CrossSection_%03d", path, (int)Nsnap);
  myfopen(fp,buf,"r");
  sprintf(buf, "%s/new/CrossSection_%03d", path, (int)Nsnap);
  myfopen(fpout,buf,"w");
  int i;
  HBTInt NumNode, NumNodeLeft;
  fread(&NumNode, sizeof(HBTInt), 1, fp);
  fwrite(&NumNode, sizeof(HBTInt), 1, fpout);
  BranchNodeBug* buffer=mymalloc(sizeof(BranchNodeBug)*BLOCKSIZE);
  BranchNode *bufferout=mymalloc(sizeof(BranchNode)*BLOCKSIZE);
  NumNodeLeft=NumNode;
  while(NumNodeLeft>BLOCKSIZE)
  {
      fread(buffer, sizeof(BranchNodeBug), BLOCKSIZE, fp);
      NumNodeLeft-=BLOCKSIZE;
      for(i=0;i<BLOCKSIZE;i++)
        copynode(buffer+i, bufferout+i);  
      fwrite(bufferout, sizeof(BranchNode), BLOCKSIZE, fpout);
  }
  if(NumNodeLeft)
  {
      fread(buffer, sizeof(BranchNodeBug), NumNodeLeft, fp);
      for(i=0;i<NumNodeLeft;i++)
        copynode(buffer+i, bufferout+i);  
      fwrite(bufferout, sizeof(BranchNode), NumNodeLeft, fpout);
  }
  
  HBTInt n;
  fread(&n, sizeof(HBTInt), 1, fp);
  fclose(fp);
  if(n!=NumNode)
  {
	fprintf(logfile, "Error: numbers of nodes do not match in %s\n"HBTIFMT","HBTIFMT"\nFile corruption?\n", buf, NumNode, n);
	exit(1);
  }
  fwrite(&NumNode, sizeof(HBTInt), 1, fpout);
  fclose(fpout);
}

int main(int argc,char **argv)
{
    printf("%zd, %zd, %zd, %zd\n", sizeof(BranchNodeBug), sizeof(BranchNode), sizeof(HBTReal), sizeof(HBTInt));
    int	Nsnap=atoi(argv[1]);
    int Nskip=atoi(argv[2]);
    int isnap;
    for(isnap=Nsnap;isnap<100;isnap+=Nskip)
        fixfile(isnap, "/mnt/ddnfs/jxhan/6610BranchTable");
   
	return 0;
}