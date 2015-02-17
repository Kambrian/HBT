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
		HBTInt subid;
		HBTInt prosubid;
		HBTInt sublen;
		HBTInt foflen;
		HBTInt subrank;
		HBTInt SnapInfall;
		HBTInt grpidInfall;
		HBTInt IsLastQuasi;
		HBTInt FlagSp; //splinter branch, was a splinter at some time
	} *prosub;
	HBTInt i,subid,Nsnap,Nsub,flag_continue,Ncen,Ncrs;
	HBTInt *sp2pro,*lastlen,*Ncross;
	HBTInt Npro,Nsplitter;
	
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	
	FILE *fp;
	char buf[1024];
	sprintf(outputdir,"%s/anal/steller",SUBCAT_DIR);	
	mkdir(outputdir,0755);		
	logfile=stdout;
	
	if(argc!=2)
	{
		printf(" %s [Nsnap]\n",argv[0]);
		exit(1);
	}
	Nsnap=atoi(argv[1]);
	
	sprintf(buf, "%s/SnapInfall_first_%03d",outputdir,(int)Nsnap);
	myfopen(fp,buf,"w");

	
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	Nsub=SubCat.Nsubs;
	prosub=mymalloc(sizeof(*prosub)*Nsub);
	Ncross=mymalloc(sizeof(HBTInt)*Nsub);
	lastlen=mymalloc(sizeof(HBTInt)*Nsub);
	for(i=0;i<Nsub;i++)
	{
		Ncross[i]=0;
		prosub[i].subid=i;
		prosub[i].prosubid=SubCat.HaloChains[i].ProSubID;
		prosub[i].FlagSp=0; //init
		if(SubCat.HaloChains[i].HostID<0)//quasi
		{
			prosub[i].subrank=0;
			prosub[i].IsLastQuasi=1;
		}
		else
		{
			if(!(prosub[i].subrank=SubCat.SubRank[i]))// a central 
			{
			prosub[i].sublen=SubCat.SubLen[i];
			prosub[i].foflen=((SubCat.HaloChains[i].HostID<0)?(SubCat.SubLen[i]):(Cat.Len[SubCat.HaloChains[i].HostID]));
			prosub[i].SnapInfall=Nsnap;
			prosub[i].grpidInfall=SubCat.HaloChains[i].HostID;
			}
			prosub[i].IsLastQuasi=0;
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
			if(prosub[i].prosubid>=0)//have a progenitor
			{
				if(prosub[i].prosubid>=SubCat.Nsubs) prosub[i].FlagSp=1;
				subid=((prosub[i].prosubid<SubCat.Nsubs)?(prosub[i].prosubid):(sp2pro[prosub[i].prosubid]));
				prosub[i].subid=subid;
				prosub[i].prosubid=SubCat.HaloChains[subid].ProSubID;
				if(prosub[i].prosubid>=0) flag_continue=1;//have something to go on
				if(SubCat.HaloChains[subid].HostID>=0) //not Quasi; skip those quasi-steps
				{
					if(prosub[i].subrank||prosub[i].IsLastQuasi)//a satellite in last snapshot, or a quasi in final snapshots
					{
						if(!(prosub[i].subrank=SubCat.SubRank[subid])) //ok, becomes a central at this snap
						{
							prosub[i].sublen=SubCat.SubLen[subid];
							prosub[i].foflen=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]):(Cat.Len[SubCat.HaloChains[subid].HostID]));
							prosub[i].SnapInfall=Nsnap;
							prosub[i].grpidInfall=SubCat.HaloChains[subid].HostID;
						}
						prosub[i].IsLastQuasi=0;//last-quasi problem solved.
					}
					else if(prosub[i].subrank=SubCat.SubRank[subid]) // a central in last snap and a satellite at this snap
						Ncross[i]++;
				}
			}
		}
		free_sub_table(&SubCat);
		free_catalogue(&Cat);
		free_sp2pro(sp2pro,Npro,Nsplitter);
	}

	fprintf(fp,"LsubInfl,LfofInfl,SnapInfl,GrpIdInfl,Subidborn,Lsub0,FlagSp,Ncross\n");
	Ncen=Ncrs=0;
	for(i=0;i<Nsub;i++)
	{
		fprintf(fp,""HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT"\n",
		prosub[i].sublen,prosub[i].foflen,prosub[i].SnapInfall,prosub[i].grpidInfall,prosub[i].subid,lastlen[i],prosub[i].FlagSp,Ncross[i]);//here the subid is not at SnapInfall,but at born time
		if(Ncross[i])
		{
			Ncen+=Ncross[i];
			Ncrs++;
		}
	}
	fclose(fp);
	printf("%f,%f\n",(float)Ncen/(float)Nsub,(float)Ncrs/(float)Nsub);
	//result:0.365951,0.180203 for 6702
	return 0;
}
