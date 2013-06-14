//to extract birth,death,quasi,splinter numbers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int read_NDeathSp(int Nsnap,char *dir)
{
FILE *fd;
char buf[1024];
int n;

sprintf(buf, "%s/srccat_%03d", dir, Nsnap);
myfopen(fd,buf,"r");
fseek(fd,-4L,SEEK_END);
fread(&n,sizeof(int),1,fd);
fclose(fd);
return n;
}

int main(int argc, char** argv)
{
	SUBCATALOGUE SubCat;
	SRCCATALOGUE SrcCat;
	
	int Nsnap=0;
	int N;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	FILE *fp;
	char buf[1024];
	sprintf(buf,"%s/anal/populations.dat",SUBCAT_DIR);
	myfopen(fp,buf,"w");
	fprintf(fp,"Nsnap\tNgroups\tNsubs\tNbirth\tNdeath\tNquasi\tNsp\tNdeathsp\n");
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
		Nsnap,SubCat.Ngroups,SubCat.Nsubs,SubCat.Nbirth,SubCat.Ndeath,SubCat.NQuasi,SubCat.Nsplitter);
	free_sub_table(&SubCat);
	N=read_NDeathSp(Nsnap,SUBCAT_DIR);
	fprintf(fp,"%d\n",N);
	}
	return 0;
}
