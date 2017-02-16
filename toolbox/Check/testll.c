#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

// #define np 500
#define NDIV 256
#undef PERIODIC_BDR
// #include "linkedlist.c"

HBTInt print_cell(HBTInt id, LINKLIST *l);
// HBTInt print_cell2(HBTInt id);
int main(int argc, char** argv)
{	
	logfile=stdout;//redirect BT routines' log info to standard output	

	HBTInt np=NP_DM;// atoi(argv[1]);
// 	HBTInt ndiv=NDIV;
	HBTInt id=0;
	
	HBTxyz *Pos=calloc(np*3, sizeof(HBTReal));
	HBTInt i,j;
	/*
#pragma omp parallel for
	for(i=0;i<np;i++)
	{
	  for(j=0;j<3;j++)
	    Pos[i][j]=drand48()*2;
	}
	*/
// 	makell(Pos, np, NDIV); 
	
	LINKLIST l, l2;
	make_linklist(&l, np,NDIV, Pos, GetArrPos, 1);
// 	omp_set_num_threads(1);
// 	make_linklist(&l2, np,NDIV, Pos, GetArrPos, 0);
	
	print_cell(id, &l);
// 	print_cell(id, &l2);
// 	print_cell2(id);
	
	free_linklist(&l);
// 	free_linklist(&l2);
	myfree(Pos);
	
	return 0;
}

HBTInt print_cell(HBTInt id, LINKLIST *l)
{
  HBTInt n=0;
  HBTInt p=l->hoc[id];
  while(p>=0)
  {
//     printf("%d ", p);
    p=l->list[p];
    n++;
  }
  printf("\n%ld particles\n", n);
  return n;
}
/*
HBTInt print_cell2(HBTInt id)
{
  HBTInt n=0;
  HBTInt p=hoc[0][0][0];
  while(p>=0)
  {
//     printf("%d ", p);
    p=ll[p];
    n++;
  }
  printf("\n%ld particles\n", n);
  return n;
}*/