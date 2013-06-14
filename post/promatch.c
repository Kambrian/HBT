#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "globals.c"
#include "mymath.c"
#include "load_group.c"
#include "subfindread.c"

#define SUBBLEN(x) ((x)<0?-1:(SubCatB.SubLen[x]))
#define SUBBDEST(x) ((x)<0?-1:(DestSub[x]))

struct prolink
{
	int id;
	int count;
};
 int compare_count(const void *a, const void *b)//used to sort integer in descending order
{
  if(((struct prolink *)a)->count > ((struct prolink*)b)->count )
	return -1;

 if(((struct prolink *)a)->count < ((struct prolink *)b)->count)
    return +1;

  return 0;
};

int main(int nargc, char ** argv)
{
char subcatSdir[512]="/SANdisk1/wenting/mergetree/SIM6702/subcatS";
char subdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5";
char fofdir[512]="/raid1/hywang/ReSim/SIM6702/group_catalogue"; //"/home/kambrain/fof_hy";
char snapdir[512]="/SANdisk4/data/NewR/SIM6702";
char outputdir[1024]="/SANdisk5/kambrain/Sim6702/SubCat5/anal";

CATALOGUE CatA,CatB;
SUBCATALOGUE SubCatA,SubCatB,SubCatC;
int Nsnap,i,j,subid,prosubid;
struct prolink *sublinkcount;
int *PartHost,*DestSub,*sp2pro;
int sublist[15]={98,161,222,
								177,178,207,234,272,
								106,322,
								202,127,85,110,139};



char buf[1024];
FILE *fp;
logfile=stdout;

if(nargc!=2)
{
printf("usage: promatch <SnapNum>\n");
return 1;
}
Nsnap=atoi(argv[1]);

	//~ load_group_catalogue(Nsnap,&CatA,fofdir);
	//~ load_group_catalogue(Nsnap,&CatB,fofdir);
	subread(99,subcatSdir);
	subfindcat2BT(&SubCatA);
	load_sub_catalogue(Nsnap,&SubCatB,subdir);
	load_sub_catalogue(Nsnap+1,&SubCatC,subdir);
	if(SubCatC.Nsplitter)
		sp2pro=load_sp2pro(Nsnap,SubCatC.Nsplitter,subdir);
	//~ load_particle_data(Nsnap,&CatA,&SubCatA,snapdir);
	//~ load_particle_data(Nsnap,&CatB,&SubCatB,snapdir);

	PartHost=mymalloc(sizeof(int)*NP_SIM);
	PartHost--;
	for(i=1;i<=NP_SIM;i++)
		PartHost[i]=-1;	
	//~ PartHost--;//so that it's accessed through [1,Ndm]
for(subid=0;subid<SubCatB.Nsubs;subid++)
{
	for(i=0;i<SubCatB.SubLen[subid];i++)
		PartHost[SubCatB.PSubArr[subid][i]]=subid;
}
	DestSub=mymalloc(sizeof(int)*SubCatB.Nsubs);
	for(i=0;i<SubCatB.Nsubs;i++)
		DestSub[i]=-1;
	for(subid=0;subid<SubCatC.Nsubs;subid++)
	{
		if((prosubid=SubCatC.ProSubID[subid])>=0)
		{
			if(prosubid<SubCatB.Nsubs)
			{
				if(DestSub[prosubid]==-1)//not a splitter and no other dest
					DestSub[prosubid]=subid;
			}
			else
			{
				prosubid=sp2pro[prosubid-SubCatB.Nsubs];
				DestSub[prosubid]=-2;// a splitter with multiple dest
			}
		}
	}
		
	sprintf(buf,"%s/subfind_lostpro_%03d",outputdir,Nsnap);
	if((fp=fopen(buf,"w"))==NULL)
	{
		printf("error: file open failed for %s!\n",buf);
		exit(1);
	}

	sublinkcount=mymalloc(sizeof(struct prolink)*(SubCatB.Nsubs+1));
	sublinkcount++;
for(i=0;i<15;i++)
{
	subid=sublist[i]-1;
	for(prosubid=-1;prosubid<SubCatB.Nsubs;prosubid++)
	{
		sublinkcount[prosubid].id=prosubid;
		sublinkcount[prosubid].count=0;
	}
	for(j=0;j<SubCatA.SubLen[subid];j++)
	{
		prosubid=PartHost[SubCatA.PSubArr[subid][j]];
		(sublinkcount[prosubid].count)++;
	}
	qsort(sublinkcount-1,SubCatB.Nsubs+1,sizeof(struct prolink),compare_count);
	fprintf(fp,"%d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\n",sublinkcount[-1].id,sublinkcount[0].id,sublinkcount[1].id,
																																	sublinkcount[-1].count,sublinkcount[0].count,sublinkcount[1].count,
																																	SUBBLEN(sublinkcount[-1].id),SUBBLEN(sublinkcount[0].id),SUBBLEN(sublinkcount[1].id),
																																	SUBBDEST(sublinkcount[-1].id),SUBBDEST(sublinkcount[0].id),SUBBDEST(sublinkcount[1].id),
																																	subid,SubCatA.SubLen[subid]);																																
}		
fclose(fp);
return 0;
}
