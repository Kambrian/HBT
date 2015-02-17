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
	SUBCATALOGUE SubCat;
	int Nsubs,*pro2dest;
	
	int Nsnap,SnapLoad;
	int i,j,grpid,subid,grplen,*members,*mstbnds;
	char outputdir[1024];
	char buf[1024];
	FILE *fp;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=3)
	{printf("usage: %s [SnapLoad] [groupid]\n",argv[0]);fflush(stdout);exit(1);}
	else
	{
		SnapLoad=atoi(argv[1]);
		grpid=atoi(argv[2]);
	}
	
	sprintf(outputdir,"%s/anal/follow_group",SUBCAT_DIR);
	mkdir(outputdir,0755);
	sprintf(outputdir,"%s/S%dG%d",outputdir,SnapLoad,grpid);
	mkdir(outputdir,0755);
	
	load_sub_catalogue(SnapLoad,&SubCat,SUBCAT_DIR);
	grplen=SubCat.GrpLen_Sub[grpid];
	members=mymalloc(sizeof(int)*grplen);
	mstbnds=mymalloc(sizeof(int)*grplen);
	for(i=0;i<grplen;i++)
	{
		subid=SubCat.GrpOffset_Sub[grpid]+i;
		members[i]=subid;
		mstbnds[i]=SubCat.PSubArr[subid][0];
	}
	
	#define SUBLEN(i) (i<0?0:SubCat.SubLen[i])
	#define SUBCOM(i,j) (i<0?0:SubCat.Property[i].CoM[j])
	FILE *fph;
	sprintf(buf,"%s/host",outputdir);
	myfopen(fph,buf,"w");
	for(Nsnap=SnapLoad;Nsnap<MaxSnap;Nsnap++)
	{
		printf("SNAP%d\n",Nsnap);fflush(stdout);
		load_particle_data(Nsnap,SNAPSHOT_DIR);
		fresh_ID2Index(mstbnds,grplen);
		
		int hostid;
		if(members[0]<0)
		hostid=-1;
		else
		hostid=SubCat.GrpOffset_Sub[SubCat.HaloChains[members[0]].HostID];
		fprintf(fph,"%d,%d,%g,%g,%g,%g\n",Nsnap,SUBLEN(hostid),comoving_virial_radius(SUBLEN(hostid)),
		SUBCOM(hostid,0),SUBCOM(hostid,1),SUBCOM(hostid,2));
		
		sprintf(buf,"%s/%d",outputdir,Nsnap);
		myfopen(fp,buf,"w");
		for(i=0;i<grplen;i++)
		{
			subid=members[i];
			fprintf(fp,"%d",SUBLEN(subid));
			fprintf(fp,",%g",comoving_virial_radius(SUBLEN(subid)));
			for(j=0;j<3;j++)
			fprintf(fp,",%g",SUBCOM(subid,j));
			for(j=0;j<3;j++)
			fprintf(fp,",%g",Pdat.Pos[mstbnds[i]][j]);
			fprintf(fp,"\n");
		}
		fclose(fp);
		erase_sub_catalogue(&SubCat);
		for(i=0;i<grplen;i++)//restore mstbound particle id
		#ifdef PID_ORDERED
			mstbnds[i]++;
		#else
			mstbnds[i]=Pdat.PID[mstbnds[i]];
		#endif
		
		if(Nsnap<MaxSnap-1)
		{
			load_pro2dest(Nsnap,&pro2dest,&Nsubs,SUBCAT_DIR);
			load_sub_catalogue(Nsnap+1,&SubCat,SUBCAT_DIR);			
			for(i=0;i<grplen;i++)
			{
				subid=pro2dest[members[i]];
				members[i]=subid;
				if(subid>=0&&SubCat.SubLen[subid]>0)
				mstbnds[i]=SubCat.PSubArr[subid][0];
			}
			free_pro2dest(pro2dest);
		}
		
	}
	fclose(fph);
	
	return 0;
}
