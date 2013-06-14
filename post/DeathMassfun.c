/* mass function for those subhalos dying
 * 
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
int *pro2dest,Nsubs,SnapPro,subid;
char buf[1024];
FILE *fp;

if(argc!=2)
{
printf("usage:%s [SnapPro]\n",argv[0]);
exit(1);
}
SnapPro=atoi(argv[1]);

sprintf(buf,"%s/anal/DeathMassfun_%03d",SUBCAT_DIR,SnapPro);
myfopen(fp,buf,"w");
load_pro2dest(SnapPro,&pro2dest,&Nsubs,SUBCAT_DIR);
load_sub_catalogue(SnapPro,&SubCat,SUBCAT_DIR);
for(subid=0;subid<SubCat.Nsubs;subid++)
{
	if(pro2dest[subid]<0)
	{
		fprintf(fp,"%d,%d,%d\n",SubCat.SubLen[subid],SubCat.HaloChains[subid].HostID,subid);
	}
}
fclose(fp);
erase_sub_catalogue(&SubCat);
free_pro2dest(pro2dest);

return 0;
}

