//simple program to extract subhalo's distance with respect to host halo, for the MADHalos test data
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
char buf[1024];FILE *fp;
int Nsnap,i;
float r_mostbnd[MaxSnap],r_COM[MaxSnap];
	
logfile=stdout;
sprintf(buf,"%s/anal/dist",SUBCAT_DIR);
myfopen(fp,buf,"w");
fprintf(fp,"Nsnap,D_MostBound,D_CoM\n");	
for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fresh_ID2Index(&SubCat,FRSH_SUBCAT);
	r_COM[Nsnap]=distance(SubCat.Property[0].CoM,SubCat.Property[1].CoM);
	int pid0,pid1;
	pid0=SubCat.PSubArr[0][0];
	pid1=SubCat.PSubArr[1][0];
	r_mostbnd[Nsnap]=distance(Pdat.Pos[pid0],Pdat.Pos[pid1]);
	fprintf(fp,"%d,%g,%g\n",Nsnap,r_mostbnd[Nsnap],r_COM[Nsnap]);
	erase_sub_catalogue(&SubCat);
}
fclose(fp);

return 0;
}
