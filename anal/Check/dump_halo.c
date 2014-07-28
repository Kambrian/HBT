#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc, char** argv)
{
	CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	HBTInt Nsubs,Nhalo;
	
	HBTInt Nsnap=0;
	HBTInt grpid,subid,pid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
 load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
FILE *fp;
char buf[1024];
sprintf(buf,"%s/anal/halo_%d",SUBCAT_DIR,snaplist[Nsnap]);
myfopen(fp,buf,"w");
fprintf(fp,"ID\tNFoF\tNBound\tX\tY\tZ\tVx\tVy\tVz\n");
for(grpid=0;grpid<Cat.Ngroups;grpid++)
{
  subid=SubCat.GrpOffset_Sub[grpid];
  fprintf(fp,"%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n",grpid,Cat.Len[grpid],SubCat.SubLen[subid],
	  SubCat.Property[subid].CoM[0],SubCat.Property[subid].CoM[1],SubCat.Property[subid].CoM[2],
	  SubCat.Property[subid].VCoM[0],SubCat.Property[subid].VCoM[1],SubCat.Property[subid].VCoM[2]);
}
fclose(fp);
 
	free_catalogue(&Cat);
	free_sub_table(&SubCat);
	return 0;
}
