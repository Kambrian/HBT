#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

int main(int argc,char **argv)
{
	FILE *fp;
	char buf[1024];
	int Nsnap,subid,fofid,pid,pid0,prosubid;
	int Lfofinf,flag_death;
	SUBCATALOGUE SubCat;
	SRCCATALOGUE SrcCat;
	CATALOGUE Cat;
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal",SUBCAT_DIR);	
	
	logfile=stdout;
	//~ Nsnap=19;
	//~ Lfofinf=154439;
	//~ prosubid=4590;
if(argc!=3){printf("usage: %s <Nsnap> <subid>\n",argv[0]);exit(1);}
Nsnap=atoi(argv[1]);
prosubid=atoi(argv[2]);

	load_group_catalogue(Nsnap,&Cat,fofdir);
	//~ for(fofid=0;fofid<Cat.Ngroups;fofid++)
	//~ {
		//~ if(Cat.Len[fofid]==Lfofinf) break;
	//~ }
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	load_src_catalogue(Nsnap,&SrcCat,inputdir);
	//~ prosubid=SubCat.GrpOffset_Sub[fofid];
	fofid=SubCat.HaloChains[prosubid].HostID;
	printf("Snap:\t(Host\tSub\tSubRank)\t(HostLen\tSubLen)\t(SrcLen\tSrcLen2\tCoreFrac)\n%d:\t(%d\t%d\t%d)\t(%d\t%d)\t(%d\t%d\t%f)\n",
			Nsnap,fofid,prosubid,SubCat.SubRank[prosubid],Cat.Len[fofid],SubCat.SubLen[prosubid],SrcCat.SubLen[prosubid],SrcCat.SubLen2[prosubid],SrcCat.CoreFrac[prosubid]);
	
	for(Nsnap++;Nsnap<100;Nsnap++)
	{
	flag_death=1;
	free_catalogue(&Cat);
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
		myfree(SubCat.PSubArr[subid]);
	}
	free_sub_catalogue(&SubCat);
	for(subid=0;subid<SrcCat.Nsubs;subid++)
	{
		myfree(SrcCat.PSubArr[subid]);
		myfree(SrcCat.PSubArr2[subid]);
	}
	free_src_catalogue(&SrcCat);
	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	load_src_catalogue(Nsnap,&SrcCat,inputdir);
		for(subid=0;subid<SubCat.Nsubs;subid++)
		{
			if(SubCat.HaloChains[subid].ProSubID==prosubid)
			{
				fofid=SubCat.HaloChains[subid].HostID;
				printf("%d:\tID(G%d\tS%d\tR%d)\tLen(G%d\tS%d)\tSrc(L%d\tLL%d\tF%f)\n",Nsnap,fofid,subid,SubCat.SubRank[subid],Cat.Len[fofid],SubCat.SubLen[subid],SrcCat.SubLen[subid],SrcCat.SubLen2[subid],SrcCat.CoreFrac[subid]);
				prosubid=subid;
				flag_death=0;
				break;
			}
		}
		if(flag_death)
		{
				printf("die\n");	
				break;
		}
	}
return 0;
}
