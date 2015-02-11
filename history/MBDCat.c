/*To output Most-bound particle catalogue (including orphan-galaxies) */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "history_vars.h"

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;	
MBDCATALOGUE MbdCat;
HBTInt Nsnap,i,pid,subid,*Sub2Hist,Nsubs;
double sqa;
char buf[1024];
HBTInt SnapBegin;

SnapBegin=0; //default
if(argc>1)	SnapBegin=atoi(argv[1]);

sprintf(buf,"%s/mostbound",SUBCAT_DIR);
mkdir(buf,0755);
sprintf(buf,"%s/mostbound/logfile",SUBCAT_DIR);
myfopen(logfile,buf,"w");

if(0==SnapBegin)
{
MbdCat.NSubs=0;
MbdCat.NOrphans=0;
MbdCat.NQuasi=0;
MbdCat.Nodes=NULL;
MbdCat.GrpLen_Sub=NULL;
MbdCat.GrpOffset_Sub=NULL;
MbdCat.GrpLen_Orphan=NULL;
MbdCat.GrpOffset_Orphan=NULL;
myfopen(logfile,buf,"w");
}
else if(SnapBegin>0)
{
	load_mbd_catalogue(SnapBegin-1,&MbdCat);
	myfopen(logfile,buf,"a");
}
else
{
	fprintf(stderr,"usage: %s [SnapBegin]\n", argv[0]);
	exit(1);
}

for(Nsnap=SnapBegin;Nsnap<MaxSnap;Nsnap++)
{
	fprintf(logfile,"loading Nsnap=%d...", (int)Nsnap);fflush(logfile);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fill_PIDHash();
	fprintf(logfile,"loaded\n");fflush(logfile);
	fresh_MBDID2Index(&MbdCat);
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(header.time);
	#endif
	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	fresh_ID2Index(&SubCat,FRSH_SUBCAT);
	
	CollectOrphans(Nsnap,&SubCat,&MbdCat);
	
	//convert subcat to mbdcat
	load_sub2hist(Nsnap,&Sub2Hist,&Nsubs,SUBCAT_DIR);
	for(subid=0;subid<SubCat.Nsubs;subid++)
	{
	if(SubCat.SubLen[subid])	
	MbdCat.Nodes[subid].MBD_PID=SubCat.PSubArr[subid][0];
	else
	MbdCat.Nodes[subid].MBD_PID=-1;/////////////////////////////
	MbdCat.Nodes[subid].HistID=Sub2Hist[subid];
	MbdCat.Nodes[subid].HostID=SubCat.HaloChains[subid].HostID;
	}
	MbdCat.NSubs=SubCat.Nsubs;
	MbdCat.NQuasi=SubCat.NQuasi;
	myfree(Sub2Hist);
	
	//fill Len and Offset
	MbdCat.Ngroups=SubCat.Ngroups;
	myfree(MbdCat.GrpLen_Sub);
	myfree(MbdCat.GrpOffset_Sub);
	myfree(MbdCat.GrpLen_Orphan);
	myfree(MbdCat.GrpOffset_Orphan);
	MbdCat.GrpLen_Sub=mymalloc(sizeof(HBTInt)*SubCat.Ngroups);
	MbdCat.GrpOffset_Sub=mymalloc(sizeof(HBTInt)*SubCat.Ngroups);
	MbdCat.GrpLen_Orphan=mymalloc(sizeof(HBTInt)*MbdCat.Ngroups);
	MbdCat.GrpOffset_Orphan=mymalloc(sizeof(HBTInt)*MbdCat.Ngroups);
	memcpy(MbdCat.GrpLen_Sub,SubCat.GrpLen_Sub,sizeof(HBTInt)*SubCat.Ngroups);
	memcpy(MbdCat.GrpOffset_Sub,SubCat.GrpOffset_Sub,sizeof(HBTInt)*SubCat.Ngroups);
	HBTInt grpid,newgrpid,GrpLen,GrpOffset;
	grpid=0;GrpLen=0;GrpOffset=MbdCat.NSubs;
	for(subid=MbdCat.NSubs;subid<MbdCat.NSubs+MbdCat.NOrphans;subid++)
	{
		newgrpid=MbdCat.Nodes[subid].HostID;
		if(newgrpid<0) break; //stop at OrphanQuasi
		if(grpid==newgrpid)
			GrpLen++;
		else if(grpid>newgrpid)
		{
			fprintf(logfile,"Err: grpid not sorted! %d, %d\n",grpid,newgrpid);
			fflush(logfile);
			exit(1);
		}
		else  //new host
		{
			MbdCat.GrpLen_Orphan[grpid]=GrpLen;
			MbdCat.GrpOffset_Orphan[grpid]=GrpOffset;
			GrpOffset+=GrpLen;
			for(grpid++;grpid<newgrpid;grpid++) //skip empty groups
			{
				MbdCat.GrpLen_Orphan[grpid]=0;
				MbdCat.GrpOffset_Orphan[grpid]=GrpOffset;
			}
			//init new grplen
			GrpLen=1;
		}
	}
	//finalize
	if(GrpLen)  //only do this when you do have some groups
	{
	MbdCat.GrpLen_Orphan[grpid]=GrpLen;
	MbdCat.GrpOffset_Orphan[grpid]=GrpOffset;
	GrpOffset+=GrpLen;
	}
	if(GrpOffset!=MbdCat.NSubs+MbdCat.NOrphans-MbdCat.NOrphanQuasi)
	{
		fprintf(logfile,"error:# of subhalos in groups not expected, "HBTIFMT", "HBTIFMT, GrpOffset,MbdCat.NSubs+MbdCat.NOrphans-MbdCat.NOrphanQuasi);
		fflush(logfile);
		exit(2);
	}
	for(grpid++;grpid<MbdCat.Ngroups;grpid++)//skip remaining groups
	{
		MbdCat.GrpLen_Orphan[grpid]=0;
		MbdCat.GrpOffset_Orphan[grpid]=GrpOffset;
	}
	erase_sub_catalogue(&SubCat);
	free_PIDHash();
	
	//fill pos and vel, and restore Ind2ID
	for(subid=0;subid<MbdCat.NSubs+MbdCat.NOrphans;subid++)
	{
		pid=MbdCat.Nodes[subid].MBD_PID;
		if(pid<0)
		{
			for(i=0;i<3;i++)
				MbdCat.Nodes[subid].Pos[i]=0.;
			for(i=0;i<3;i++)
				MbdCat.Nodes[subid].Vel[i]=0.;
			MbdCat.Nodes[subid].MBD_PID=-1;		
		}
		else
		{
			for(i=0;i<3;i++)
				MbdCat.Nodes[subid].Pos[i]=Pdat.Pos[pid][i];
			for(i=0;i<3;i++)
				MbdCat.Nodes[subid].Vel[i]=Pdat.Vel[MbdCat.Nodes[subid].MBD_PID][i]*sqa;
			MbdCat.Nodes[subid].MBD_PID=Pdat.PID[MbdCat.Nodes[subid].MBD_PID];	
		}	
	}
	
	save_mbd_catalogue(Nsnap, &MbdCat);		
	free_particle_data();

}

free_mbd_catalogue(&MbdCat);
fclose(logfile);

return 0;

}

static int comp_NodeHost(const void *a, const void *b)//used to sort HostID in ascending order; 
{//put HostID<0 to the end (as if HostID=+inf)
//or better put it to the beginning, together with quasi-subs???????????????????????????????
  if(((MBDNode *)a)->HostID == ((MBDNode *) b)->HostID)
    return 0;
	
  if(((MBDNode *)a)->HostID <0 ) 
    return +1;

  if(((MBDNode *) b)->HostID <0) 
    return -1;

  if(((MBDNode *)a)->HostID > ((MBDNode *) b)->HostID)
    return +1;
	
  //if(((MBDNode *)a)->HostID < ((MBDNode *) b)->HostID)	
   return -1;
}

HBTInt CollectOrphans(HBTInt Nsnap, SUBCATALOGUE *SubCat, MBDCATALOGUE *MbdCat)
{//to find orphans subs for MbdCat inherited from Nsnap-1, using subcat and cat at Nsnap.
//update MbdCat.Nodes: MBD_PID,HistID,HostID, and NSubs,NOrphans
//Note: MBD_PID are PIndex rather than PID when returned
	CATALOGUE Cat;
	HBTInt *pro2dest, Npro;
	HBTInt subid,NOrphans,NOrphans_new;
	MBDNode *Orphans,*NewNodes; 
	
	//find new orphans
	load_pro2dest(Nsnap-1,&pro2dest,&Npro,SUBCAT_DIR);
	NewNodes=mymalloc(sizeof(MBDNode)*(SubCat->Nsubs+MbdCat->NSubs+MbdCat->NOrphans));//leave space for new living subs
	Orphans=NewNodes+SubCat->Nsubs;//pointer for orphans
	NOrphans=0;
	for(subid=0;subid<MbdCat->NSubs;subid++)
	{
		if(pro2dest[subid]<0)
		{
			Orphans[NOrphans].MBD_PID=MbdCat->Nodes[subid].MBD_PID;
			Orphans[NOrphans].HistID=MbdCat->Nodes[subid].HistID;
			NOrphans++;
		}
	}
	free_pro2dest(pro2dest);
	//add existing orphans
	NewNodes=realloc(NewNodes,sizeof(MBDNode)*(SubCat->Nsubs+NOrphans+MbdCat->NOrphans));
	Orphans=NewNodes+SubCat->Nsubs;//pointer for orphans
	for(subid=MbdCat->NSubs;subid<MbdCat->NSubs+MbdCat->NOrphans;subid++)
	{
		Orphans[NOrphans].MBD_PID=MbdCat->Nodes[subid].MBD_PID;
		Orphans[NOrphans].HistID=MbdCat->Nodes[subid].HistID;
		NOrphans++;
	}
	//find hosts
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	fresh_ID2Index(&Cat,FRSH_GRPCAT); 	
	prepare_ind2halo(&Cat);
	for(subid=0;subid<NOrphans;subid++)
		Orphans[subid].HostID=Cat.ID2Halo[Orphans[subid].MBD_PID];
	free_catalogue(&Cat);
	qsort(Orphans,NOrphans,sizeof(MBDNode),comp_NodeHost);
	for(subid=NOrphans-1;subid>=0;subid--)
	{
		if(Orphans[subid].HostID>=0) break; //bug-fix: >0 to >=0; 10-02-2015.
	}
	MbdCat->NOrphanQuasi=NOrphans-1-subid;
	myfree(MbdCat->Nodes);
	MbdCat->Nodes=NewNodes;
	MbdCat->NSubs=SubCat->Nsubs;
	MbdCat->NOrphans=NOrphans;	
	return MbdCat->NOrphans;
}
