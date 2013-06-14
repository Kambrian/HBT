//last infall
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
	SUBCATALOGUE SubCat;
	CATALOGUE Cat;
	struct
	{
		int subid;
		int prosubid;
		int sublen;
		int foflen;
		int subrank;
		int SnapInfall;
		int DirectInfall;// 1 means direct infall to the final fof; 0 means hierachical infall
		int hostsubpro;
	} *prosub;
	int i,subid,Nsnap,Nsub,flag_continue;
	int *sp2pro,*lastlen,Npro,Nsplitter;
	
	char outputdir[1024];		
	FILE *fp;
	char buf[1024];	
	logfile=stdout;
	sprintf(outputdir,"%s/anal/steller",SUBCAT_DIR);	
	mkdir(outputdir,0755);
	
	if(argc!=2)
	{
		printf(" %s [Nsnap]\n",argv[0]);
		exit(1);
	}
	Nsnap=atoi(argv[1]);
	
	sprintf(buf, "%s/SnapInfall_%03d",outputdir,Nsnap);
	myfopen(fp,buf,"w");
	
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	Nsub=SubCat.Nsubs;
	prosub=mymalloc(sizeof(*prosub)*Nsub);
	lastlen=mymalloc(sizeof(int)*Nsub);
	for(i=0;i<Nsub;i++)
	{
		prosub[i].subid=i;
		prosub[i].prosubid=SubCat.HaloChains[i].ProSubID;
		prosub[i].DirectInfall=1;
		if(!(prosub[i].subrank=SubCat.SubRank[i]))// a central 
		{
		prosub[i].sublen=SubCat.SubLen[i];
		prosub[i].foflen=((SubCat.HaloChains[i].HostID<0)?(SubCat.SubLen[i]):(Cat.Len[SubCat.HaloChains[i].HostID]));
		prosub[i].SnapInfall=Nsnap;
		}
		else
		{
			prosub[i].hostsubpro=SubCat.HaloChains[SubCat.GrpOffset_Sub[SubCat.HaloChains[i].HostID]].ProSubID;//its last fof's main sub chain;
		}
		lastlen[i]=SubCat.SubLen[i];
	}
	free_sub_table(&SubCat);
	free_catalogue(&Cat);
	
	flag_continue=Nsub;
	while(flag_continue)
	{
		load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);
		Nsnap--;
		flag_continue=0;
		load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
		load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
		for(i=0;i<Nsub;i++)
		{
			if(prosub[i].subrank)//a satellite in last snapshot
			{
				subid=((prosub[i].prosubid<SubCat.Nsubs)?(prosub[i].prosubid):(sp2pro[prosub[i].prosubid]));
				prosub[i].subid=subid;
				prosub[i].prosubid=SubCat.HaloChains[subid].ProSubID;
				if(prosub[i].subrank=SubCat.SubRank[subid])//still a satellite at this snap,continue
				{	
					flag_continue=1;
					if(prosub[i].DirectInfall)//only update DirectInfall if there has been no evidence of hierachical infall
					{
					prosub[i].DirectInfall=(SubCat.HaloChains[prosub[i].hostsubpro].HostID==SubCat.HaloChains[prosub[i].subid].HostID);
					prosub[i].hostsubpro=SubCat.HaloChains[prosub[i].hostsubpro].ProSubID;//update its last fof's main sub chain;
					}
				}
				else  //ok, becomes a central at this snap
				{
					prosub[i].sublen=SubCat.SubLen[subid];
					prosub[i].foflen=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]):(Cat.Len[SubCat.HaloChains[subid].HostID]));
					prosub[i].SnapInfall=Nsnap;
				}
			}
		}
		free_sub_table(&SubCat);
		free_catalogue(&Cat);
		free_sp2pro(sp2pro,Npro,Nsplitter);
	}

		fprintf(fp,"LsubInfl,LfofInfl,SnapInfl,SubidInfl,Lsub0,DirectInfall\n");
	for(i=0;i<Nsub;i++)
		fprintf(fp,"%d,%d,%d,%d,%d,%d\n",prosub[i].sublen,prosub[i].foflen,prosub[i].SnapInfall,prosub[i].subid,lastlen[i],prosub[i].DirectInfall);
	fclose(fp);
	return 0;
}
