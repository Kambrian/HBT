#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#define CEN_COM
HBTReal *Rmax, *Vmax;
SUBCATALOGUE SubCat;
extern void init_RmaxVmax(HBTInt Nsnap);
int main(int argc,char **argv)
{

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
		HBTReal RmaxInfall;
		HBTReal VmaxInfall;
		HBTInt FlagSp; //splinter branch, was a splinter at some time
	} *prosub;
	HBTInt i,subid,Nsnap,Nsub,flag_continue,Ncen,Ncrs;
	HBTInt *sp2pro,*lastlen,*Ncross;
	HBTInt Npro,Nsplitter;
	
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	
	FILE *fp;
	char buf[1024];
	sprintf(outputdir,"/data/A4700r2d1/kambrain/6113/post/anal");	
	mkdir(outputdir,0755);		
	logfile=stdout;
	
	if(argc!=2)
	{
		printf(" %s [Nsnap]\n",argv[0]);
		exit(1);
	}
	Nsnap=atoi(argv[1]);
	HBTInt Nsnap0=Nsnap;
	
	sprintf(buf, "%s/VmaxInfall_first_%03d",outputdir,(int)Nsnap);
	myfopen(fp,buf,"w");

	
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	init_RmaxVmax(Nsnap);
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
			prosub[i].RmaxInfall=Rmax[i];
			prosub[i].VmaxInfall=Vmax[i];
			prosub[i].grpidInfall=SubCat.HaloChains[i].HostID;
			}
			prosub[i].IsLastQuasi=0;
		}
		lastlen[i]=SubCat.SubLen[i];
	}
	myfree(Rmax);
	myfree(Vmax);
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
		init_RmaxVmax(Nsnap);
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
							prosub[i].RmaxInfall=Rmax[subid];
							prosub[i].VmaxInfall=Vmax[subid];
							prosub[i].grpidInfall=SubCat.HaloChains[subid].HostID;
						}
						prosub[i].IsLastQuasi=0;//last-quasi problem solved.
					}
					else if(prosub[i].subrank=SubCat.SubRank[subid]) // a central in last snap and a satellite at this snap
						Ncross[i]++;
				}
			}
		}
		myfree(Rmax);
		myfree(Vmax);
		free_sub_table(&SubCat);
		free_catalogue(&Cat);
		free_sp2pro(sp2pro,Npro,Nsplitter);
	}

	fprintf(fp,"SubhaloID,host-HaloID,x,y,z,Vx,Vy,Vz,NboundNow,NboundInfall,RmaxInfall,VmaxInfall,SnapInfl,GrpIdInfl,FlagSp,Ncross\n");
	load_sub_table(Nsnap0,&SubCat,SUBCAT_DIR);
	for(i=0;i<Nsub;i++)
	{
		fprintf(fp,HBTIFMT","HBTIFMT",%g,%g,%g,%g,%g,%g,"HBTIFMT","HBTIFMT",%g,%g,"HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT"\n",
		i,SubCat.HaloChains[i].HostID,SubCat.Property[i].CoM[0],SubCat.Property[i].CoM[1],SubCat.Property[i].CoM[2],
	        SubCat.Property[i].VCoM[0],SubCat.Property[i].VCoM[1],SubCat.Property[i].VCoM[2],
	        SubCat.SubLen[i],
	        prosub[i].sublen,prosub[i].RmaxInfall,prosub[i].VmaxInfall,prosub[i].SnapInfall,prosub[i].grpidInfall,prosub[i].FlagSp,Ncross[i]);//here the subid is not at SnapInfall,but at born time
		if(Ncross[i])
		{
			Ncen+=Ncross[i];
			Ncrs++;
		}
	}
	fclose(fp);
	free_sub_table(&SubCat);
	printf("%f,%f\n",(float)Ncen/(float)Nsub,(float)Ncrs/(float)Nsub);
	//result:0.365951,0.180203 for 6702
	return 0;
}

int load_RmaxVmax(HBTReal *rmax,HBTReal *vmax, HBTInt Nsnap)
{
	char buf[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;
	//~ HBTReal *rsig, *r2sig, *r3sig, *rpoisson; /* declare this as input var if you want them!!! */
	
	if(!rmax||!vmax)
	{
		printf("error: allocate rmax , vmax first \n");
		exit(1);
	}
	#ifdef CEN_COM
	sprintf(buf,"%s/../post/profile/RmaxVmax_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
	#else
	sprintf(buf,"%s/../post/profile/RmaxVmax_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	fread(rmax,sizeof(HBTReal),Nsubs,fp);
	fread(vmax,sizeof(HBTReal),Nsubs,fp);
	fseek(fp,sizeof(HBTReal)*Nsubs*5,SEEK_CUR);
// 	fread(rhalf,sizeof(HBTReal),Nsubs,fp);
	//~ fread(rsig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r2sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r3sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(rpoisson,sizeof(HBTReal),Nsubs,fp);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
	return Nsubs;
}
void init_RmaxVmax(HBTInt Nsnap)
{
	Rmax=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	Vmax=mymalloc(sizeof(HBTReal)*SubCat.Nsubs);
	load_RmaxVmax(Rmax,Vmax,Nsnap);
}
