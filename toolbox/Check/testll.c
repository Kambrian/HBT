#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

// #define NDIV 40
// #define np 500
#define NDIV 3
#undef PERIODIC_BDR
#include "linkedlist.c"

int print_cell(int id, LINKLIST *l);
int print_cell2(int id);
int main(int argc, char** argv)
{	
	logfile=stdout;//redirect BT routines' log info to standard output	

	int np=atoi(argv[1]);
// 	int ndiv=NDIV;
	int id=0;
	
	HBTxyz Pos[np];
	int i,j;
	for(i=0;i<np;i++)
	{
	  for(j=0;j<3;j++)
	    Pos[i][j]=drand48()*2;
	}
	makell(Pos, np, NDIV); 
	
	LINKLIST l, l2;
	make_linklist(&l, np,NDIV, Pos, GetArrPos, 0);
	omp_set_num_threads(1);
	make_linklist(&l2, np,NDIV, Pos, GetArrPos, 0);
	
	print_cell(id, &l);
	print_cell(id, &l2);
	print_cell2(id);
	
	free_linklist(&l);
	free_linklist(&l2);
	
	return 0;
}

int print_cell(int id, LINKLIST *l)
{
  int n=0;
  int p=l->hoc[id];
  while(p>=0)
  {
    printf("%d ", p);
    p=l->list[p];
    n++;
  }
  printf("\n%d particles\n", n);
  return n;
}

int print_cell2(int id)
{
  int n=0;
  int p=hoc[0][0][0];
  while(p>=0)
  {
    printf("%d ", p);
    p=ll[p];
    n++;
  }
  printf("\n%d particles\n", n);
  return n;
}