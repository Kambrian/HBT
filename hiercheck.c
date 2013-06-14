#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

static short *subhier_mask;
void init_subhier_mask(HBTInt nsubs)
{ HBTInt i;
	subhier_mask=mymalloc(sizeof(short)*nsubs);
	for(i=0;i<nsubs;i++)
		subhier_mask[i]=0;
}
void check_hierarchy_recursive(HBTInt mainsubid,HBTInt mainnibs,HBTInt mainpre,HBTInt mainnext,SUBCATALOGUE *SubCat)
{ 
	HBTInt nibs,pre,next,subid;
	HBTInt hostid;
	HBTInt idl,idh;
		if(SubCat->sub_hierarchy[mainsubid].nibs!=mainnibs||SubCat->sub_hierarchy[mainsubid].pre!=mainpre||SubCat->sub_hierarchy[mainsubid].next!=mainnext)
		{
			printf("error hier wrong arc\n");
			exit(11);
		}
		if(subhier_mask[mainsubid])
		{
			printf("error hier dupmask\n");
			exit(12);
		}
		subhier_mask[mainsubid]=1;
		hostid=SubCat->HaloChains[mainsubid].HostID;
		if(hostid>=0)
		{
			idl=SubCat->GrpOffset_Sub[hostid];
			idh=SubCat->GrpOffset_Sub[hostid]+SubCat->GrpLen_Sub[hostid];
		}
		else
		{
			idl=SubCat->Nsubs-SubCat->NQuasi;
			idh=SubCat->Nsubs;
		}	
		if(mainsubid<idl||mainsubid>=idh)
		{
			printf("error hier out of range\n");
			exit(13);
		}
		nibs=mainsubid;
		subid=SubCat->sub_hierarchy[nibs].sub;
		pre=-1;
		next=SubCat->sub_hierarchy[subid].next;
		while(subid>0)
		{
			check_hierarchy_recursive(subid,nibs,pre,next,SubCat);
			pre=subid;
			subid=next;
			next=SubCat->sub_hierarchy[subid].next;
		}
}
void check_subhier_mask(HBTInt nsubs)
{
	HBTInt i;
	for(i=0;i<nsubs;i++)
		if(subhier_mask[i]==0)
		{
			printf("error hier mask0\n");
			exit(14);
		}
	if(nsubs!=0) free(subhier_mask);	
}


int main(HBTInt argc,char **argv)
{
HBTInt SnapshotNum=0,i;
CATALOGUE CatB;
SUBCATALOGUE SubCat;
SRCCATALOGUE SrcCat;
HBTInt Nsnap[2];
	
	logfile=stdout;	
if(argc==2)
{
Nsnap[0]=atoi(argv[1]);
Nsnap[1]=atoi(argv[1])+1;
}
else if(argc==3)
{
	Nsnap[0]=atoi(argv[1]);
	Nsnap[1]=atoi(argv[2])+1;
}
else
{
	printf("usage: %s <snapbegin> <snapend>\n",argv[0]);
	exit(1);
}



	for(SnapshotNum=Nsnap[0];SnapshotNum<Nsnap[1];SnapshotNum++)
	{
			load_sub_catalogue(SnapshotNum,&SubCat,SUBCAT_DIR);
			load_src_catalogue(SnapshotNum,&SrcCat,SUBCAT_DIR);
			init_subhier_mask(SubCat.Nsubs);
			for(i=0;i<SubCat.Ngroups;i++)
			{
				check_hierarchy_recursive(SubCat.GrpOffset_Sub[i],-1,-1,-1,&SubCat);
			}
			for(i=SubCat.Nsubs-SubCat.NQuasi;i<SubCat.Nsubs;i++)
			{
				if(SubCat.sub_hierarchy[i].nibs<0)
					check_hierarchy_recursive(i,-1,-1,-1,&SubCat);
			}
			check_subhier_mask(SubCat.Nsubs);
			for(i=0;i<SubCat.Nsubs;i++)
			{
				free(SubCat.PSubArr[i]);
				free(SrcCat.PSubArr[i]);
				myfree(SrcCat.PSubArr2[i]);
			}
			free_sub_catalogue(&SubCat);
			free_src_catalogue(&SrcCat);
	}
printf("all check passed\n");
return 0;
}
