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

#define MassType 0 //mass definition for halo: 0~2 (vir); 3 (fof); 4 (bound mass)
#define SubMassType MassType
#define MassTolerance 0.3
#define DistanceTolerance 0.15
#define NTarget 4

void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
{
	char buf[1024];
	FILE *fp;
	int Nvir[3],i,j;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,14*4L,SEEK_CUR);
		//fread(Mvir,sizeof(int),3,fp);
		fread(Nvir,sizeof(int),3,fp);
		for(j=0;j<3;j++)
		Mvir[i][j]=Nvir[j]*partmass;
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}

struct
{
	HBTInt subid;
	HBTInt prosubid;
	HBTInt sublen;
// 	HBTInt foflen;
	HBTReal Mvir;
	HBTInt subrank;
	HBTInt SnapInfall;
	HBTInt DirectInfall;// 1 means direct infall to the final fof; 0 means hierachical infall
	HBTInt FlagSp; //1 means splinter branch
	HBTInt hostsubpro;
} *prosubLast;

struct
{
	HBTInt subid;
	HBTInt prosubid;
	HBTInt sublen;
// 	HBTInt foflen;
	HBTReal Mvir;
	HBTInt subrank;
	HBTInt SnapInfall;
	HBTInt grpidInfall;
	HBTInt IsLastQuasi;
	HBTInt FlagSp; //splinter branch, was a splinter at some time
} *prosubFirst;
void track_last_infall(int Nsnap, HBTInt *SubList, HBTInt NumSub)
{
	SUBCATALOGUE SubCat;
//	CATALOGUE Cat;
	
	HBTInt i,subid,Nsub,flag_continue;
	HBTInt *sp2pro,*lastlen,Npro,Nsplitter;
	
	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
//	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	HBTReal (*MvirHost)[3], (*RvirHost)[3];
	#if SubMassType<3
	MvirHost=mymalloc(sizeof(HBTReal)*3*SubCat.Ngroups);
	RvirHost=mymalloc(sizeof(HBTReal)*3*SubCat.Ngroups);
	load_halo_virial_size(MvirHost,RvirHost,header.mass[1],SubCat.Ngroups,Nsnap);
	#endif
	Nsub=NumSub;
	prosubLast=mymalloc(sizeof(*prosubLast)*Nsub);
	lastlen=mymalloc(sizeof(HBTInt)*Nsub);
	for(i=0;i<Nsub;i++)
	{
	  subid=SubList[i];
		prosubLast[i].subid=subid;
		prosubLast[i].prosubid=SubCat.HaloChains[subid].ProSubID;
		prosubLast[i].DirectInfall=1;
		prosubLast[i].FlagSp=0;
		if(!(prosubLast[i].subrank=SubCat.SubRank[subid]))// a central 
		{
		prosubLast[i].sublen=SubCat.SubLen[subid];
//		prosubLast[i].foflen=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]):(Cat.Len[SubCat.HaloChains[subid].HostID]));
		#if SubMassType<3
		prosubLast[i].Mvir=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]*header.mass[1]):(MvirHost[SubCat.HaloChains[subid].HostID][SubMassType]));
		#else 
		prosubLast[i].Mvir=0.;
		#endif
		prosubLast[i].SnapInfall=Nsnap;
		}
		else
		{
			prosubLast[i].hostsubpro=SubCat.HaloChains[SubCat.GrpOffset_Sub[SubCat.HaloChains[subid].HostID]].ProSubID;//its last fof's main sub chain;
		}
		lastlen[i]=SubCat.SubLen[subid];
	}
	free_sub_table(&SubCat);
//	free_catalogue(&Cat);
	#if SubMassType<3
	myfree(MvirHost);
	myfree(RvirHost);
	#endif
	flag_continue=Nsub;
	while(flag_continue)
	{
		load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);
		Nsnap--;
		flag_continue=0;
		load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
		printf("Nsnap=%d\n",Nsnap);fflush(stdout);
//		load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
		#if SubMassType<3
		MvirHost=mymalloc(sizeof(HBTReal)*3*SubCat.Ngroups);
		RvirHost=mymalloc(sizeof(HBTReal)*3*SubCat.Ngroups);
		load_halo_virial_size(MvirHost,RvirHost,header.mass[1],SubCat.Ngroups,Nsnap);
		#endif
		for(i=0;i<Nsub;i++)
		{
			if(prosubLast[i].subrank)//a satellite in last snapshot
			{
				if(prosubLast[i].prosubid>=SubCat.Nsubs) prosubLast[i].FlagSp=1; //tag as a splinter
				subid=((prosubLast[i].prosubid<SubCat.Nsubs)?(prosubLast[i].prosubid):(sp2pro[prosubLast[i].prosubid]));
				prosubLast[i].subid=subid;
				prosubLast[i].prosubid=SubCat.HaloChains[subid].ProSubID;
				if(prosubLast[i].subrank=SubCat.SubRank[subid])//still a satellite at this snap,continue
				{	
					flag_continue=1;
					if(prosubLast[i].DirectInfall)//only update DirectInfall if there has been no evidence of hierachical infall
					{
					prosubLast[i].DirectInfall=(SubCat.HaloChains[prosubLast[i].hostsubpro].HostID==SubCat.HaloChains[prosubLast[i].subid].HostID);
					prosubLast[i].hostsubpro=SubCat.HaloChains[prosubLast[i].hostsubpro].ProSubID;//update its last fof's main sub chain;
					}
				}
				else  //ok, becomes a central at this snap
				{
					prosubLast[i].sublen=SubCat.SubLen[subid];
// 					prosubLast[i].foflen=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]):(Cat.Len[SubCat.HaloChains[subid].HostID]));
					#if SubMassType<3
					prosubLast[i].Mvir=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]*header.mass[1]):(MvirHost[SubCat.HaloChains[subid].HostID][SubMassType]));
					#else
					prosubLast[i].Mvir=0.;
					#endif
					prosubLast[i].SnapInfall=Nsnap;
				}
			}
		}
		free_sub_table(&SubCat);
// 		free_catalogue(&Cat);
#if SubMassType<3
		myfree(MvirHost);
		myfree(RvirHost);
#endif
		free_sp2pro(sp2pro,Npro,Nsplitter);
	}
/*
		fprintf(fp,"LsubInfl,LfofInfl,SnapInfl,SubidInfl,Lsub0,DirectInfall\n");
	for(i=0;i<Nsub;i++)
		fprintf(fp,"%d,%d,%d,%d,%d,%d\n",prosubLast[i].sublen,prosubLast[i].foflen,prosubLast[i].SnapInfall,prosubLast[i].subid,lastlen[i],prosubLast[i].DirectInfall);
	fclose(fp);*/
// 	return prosubLast;
}


void track_first_infall(int Nsnap, HBTInt *SubList, HBTInt NumSub)
{
	SUBCATALOGUE SubCat;
// 	CATALOGUE Cat;

	HBTInt i,subid,Nsub,flag_continue,Ncen,Ncrs;
	HBTInt *sp2pro,*lastlen,*Ncross;
	HBTInt Npro,Nsplitter;

	load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
// 	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
#if SubMassType<3
	HBTReal (*MvirHost)[3], (*RvirHost)[3];
	MvirHost=mymalloc(sizeof(HBTReal)*3*SubCat.Ngroups);
	RvirHost=mymalloc(sizeof(HBTReal)*3*SubCat.Ngroups);
	load_halo_virial_size(MvirHost,RvirHost,header.mass[1],SubCat.Ngroups,Nsnap);
#endif
	Nsub=NumSub;
	prosubFirst=mymalloc(sizeof(*prosubFirst)*Nsub);
	Ncross=mymalloc(sizeof(HBTInt)*Nsub);
	lastlen=mymalloc(sizeof(HBTInt)*Nsub);
	for(i=0;i<Nsub;i++)
	{
	  subid=SubList[i];
		Ncross[i]=0;
		prosubFirst[i].subid=subid;
		prosubFirst[i].prosubid=SubCat.HaloChains[subid].ProSubID;
		prosubFirst[i].FlagSp=0; //init
		if(SubCat.HaloChains[subid].HostID<0)//quasi
		{
			prosubFirst[i].subrank=0;
			prosubFirst[i].IsLastQuasi=1;
		}
		else
		{
			if(!(prosubFirst[i].subrank=SubCat.SubRank[subid]))// a central 
			{
			prosubFirst[i].sublen=SubCat.SubLen[subid];
// 			prosubFirst[i].foflen=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]):(Cat.Len[SubCat.HaloChains[subid].HostID]));
#if SubMassType<3
			prosubFirst[i].Mvir=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]*header.mass[1]):(MvirHost[SubCat.HaloChains[subid].HostID][SubMassType]));
#else
			prosubFirst[i].Mvir=0.;
#endif
			prosubFirst[i].SnapInfall=Nsnap;
			prosubFirst[i].grpidInfall=SubCat.HaloChains[subid].HostID;
			}
			prosubFirst[i].IsLastQuasi=0;
		}
		lastlen[i]=SubCat.SubLen[subid];
	}
	free_sub_table(&SubCat);
// 	free_catalogue(&Cat);
#if SubMassType<3
	myfree(MvirHost);
	myfree(RvirHost);
#endif	
	flag_continue=Nsub;
	while(flag_continue)
	{
		load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);
		Nsnap--;
		flag_continue=0;
		load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
		printf("Nsnap=%d\n",Nsnap);fflush(stdout);
// 		load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
#if SubMassType<3
		MvirHost=mymalloc(sizeof(HBTReal)*3*SubCat.Ngroups);
		RvirHost=mymalloc(sizeof(HBTReal)*3*SubCat.Ngroups);
		load_halo_virial_size(MvirHost,RvirHost,header.mass[1],SubCat.Ngroups,Nsnap);
#endif
		for(i=0;i<Nsub;i++)
		{				
			if(prosubFirst[i].prosubid>=0)//have a progenitor
			{
				if(prosubFirst[i].prosubid>=SubCat.Nsubs) prosubFirst[i].FlagSp=1;
				subid=((prosubFirst[i].prosubid<SubCat.Nsubs)?(prosubFirst[i].prosubid):(sp2pro[prosubFirst[i].prosubid]));
				prosubFirst[i].subid=subid;
				prosubFirst[i].prosubid=SubCat.HaloChains[subid].ProSubID;
				if(prosubFirst[i].prosubid>=0) flag_continue=1;//have something to go on
				if(SubCat.HaloChains[subid].HostID>=0) //not Quasi; skip those quasi-steps
				{
					if(prosubFirst[i].subrank||prosubFirst[i].IsLastQuasi)//a satellite in last snapshot, or a quasi in final snapshots
					{
						if(!(prosubFirst[i].subrank=SubCat.SubRank[subid])) //ok, becomes a central at this snap
						{
							prosubFirst[i].sublen=SubCat.SubLen[subid];
// 							prosubFirst[i].foflen=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]):(Cat.Len[SubCat.HaloChains[subid].HostID]));
#if SubMassType<3
							prosubFirst[i].Mvir=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]*header.mass[1]):(MvirHost[SubCat.HaloChains[subid].HostID][SubMassType]));
#else
							prosubFirst[i].Mvir=0.;
#endif
							prosubFirst[i].SnapInfall=Nsnap;
							prosubFirst[i].grpidInfall=SubCat.HaloChains[subid].HostID;
						}
						prosubFirst[i].IsLastQuasi=0;//last-quasi problem solved.
					}
					else if(prosubFirst[i].subrank=SubCat.SubRank[subid]) // a central in last snap and a satellite at this snap
						Ncross[i]++;
				}
			}
		}
		free_sub_table(&SubCat);
// 		free_catalogue(&Cat);
#if SubMassType<3
		myfree(MvirHost);
		myfree(RvirHost);
#endif
		free_sp2pro(sp2pro,Npro,Nsplitter);
	}
/*
	fprintf(fp,"LsubInfl,LfofInfl,SnapInfl,GrpIdInfl,Subidborn,Lsub0,FlagSp,Ncross\n");
	Ncen=Ncrs=0;
	for(i=0;i<Nsub;i++)
	{
		fprintf(fp,""HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT"\n",
		prosubFirst[i].sublen,prosubFirst[i].foflen,prosubFirst[i].SnapInfall,prosubFirst[i].grpidInfall,prosubFirst[i].subid,lastlen[i],prosubFirst[i].FlagSp,Ncross[i]);//here the subid is not at SnapInfall,but at born time
		if(Ncross[i])
		{
			Ncen+=Ncross[i];
			Ncrs++;
		}
	}
	fclose(fp);
	printf("%f,%f\n",(float)Ncen/(float)Nsub,(float)Ncrs/(float)Nsub);
	//result:0.365951,0.180203 for 6702
	return 0;*/
}

int main(int argc,char **argv)
{
SUBCATALOGUE SubCatHost;
char buf[1024];FILE *fp;
HBTInt Nsnap,NsnapHost,subid,grpid,i,j; 
HBTReal zhost=0.3;
HBTReal TargetMass[NTarget]={1742, 883, 466, 304}, TargetDistance[2]={77.5,30}, r;
logfile=stdout;

HBTReal z[MaxSnap],zdiff[MaxSnap];
sprintf(buf,"%s/Redshift.dat",SUBCAT_DIR);
myfopen(fp,buf,"r");
for(i=0;i<IniSnap;i++)zdiff[i]=1000.;
for(i=IniSnap;i<MaxSnap;i++)
{
	fscanf(fp,"%*d\t%g\n",&(z[i]));
	z[i]=1/z[i]-1.;
	zdiff[i]=fabs(z[i]-zhost);
}
fclose(fp);
load_particle_header(IniSnap,SNAPSHOT_DIR);
// for(i=0;i<MaxSnap;i++)
// {
// 	load_particle_header(i,SNAPSHOT_DIR);
// 	z[i]=header.ztp;
// 	zdiff[i]=fabs(z[i]-zhost);
// }
NsnapHost=min_of_vec(zdiff,MaxSnap);

printf("Target Snapshot=%d, z=%.2f\n",(int)NsnapHost,z[NsnapHost]);

load_sub_table(NsnapHost,&SubCatHost,SUBCAT_DIR);
HBTReal *Mhost=mymalloc(sizeof(HBTReal)*SubCatHost.Ngroups);
#if MassType<3
HBTReal (*MvirHost)[3], (*RvirHost)[3];
MvirHost=mymalloc(sizeof(HBTReal)*3*SubCatHost.Ngroups);
RvirHost=mymalloc(sizeof(HBTReal)*3*SubCatHost.Ngroups);
load_halo_virial_size(MvirHost,RvirHost,header.mass[1],SubCatHost.Ngroups,NsnapHost);
for(i=0;i<SubCatHost.Ngroups;i++) Mhost[i]=MvirHost[i][MassType];
myfree(MvirHost);
myfree(RvirHost);
#elif MassType==3
CATALOGUE Cat;
load_group_catalogue(NsnapHost,&Cat,GRPCAT_DIR);
for(i=0;i<Cat.Ngroups;i++) Mhost[i]=Cat.Len[i]*header.mass[1];
free_catalogue(&Cat);
#elif MassType==4
for(i=0;i<SubCatHost.Ngroups;i++) Mhost[i]=SubCatHost.SubLen[SubCatHost.GrpOffset_Sub[i]]*header.mass[1];
#else
printf("error: wrong MassType\n");exit(1);
#endif
HBTInt *HostList=mymalloc(sizeof(HBTInt)*SubCatHost.Ngroups);
HBTInt NumHost=0;
#define SelectMass(i,j) (fabs(Mhost[i]-TargetMass[j])<TargetMass[j]*MassTolerance)
#define SelectDist(j) (fabs(r-TargetDistance[j])<TargetDistance[j]*DistanceTolerance)
//pick hosts
for(i=0;i<SubCatHost.Ngroups;i++)
{
 if(SelectMass(i,0)||SelectMass(i,1)||SelectMass(i,2)||SelectMass(i,3))
//     if((Mhost[i]<TargetMass[0])&&(Mhost[i]>TargetMass[1]))
      HostList[NumHost++]=i;
}
HBTInt *SubList=mymalloc(sizeof(HBTInt)*SubCatHost.GrpLen_Sub[HostList[0]]*NumHost);
HBTReal *SubDist=mymalloc(sizeof(HBTReal)*SubCatHost.GrpLen_Sub[HostList[0]]*NumHost);
HBTInt NumSub=0, hostid, hostsubid;
for(i=0;i<NumHost;i++)
{
  hostid=HostList[i];
  hostsubid=SubCatHost.GrpOffset_Sub[hostid];
  for(j=hostsubid+1;j<hostsubid+SubCatHost.GrpLen_Sub[hostid];j++)
  {
    r=distance(SubCatHost.Property[j].CoM,SubCatHost.Property[hostsubid].CoM)/(1+z[NsnapHost]);
//     if(SelectDist(0)||SelectDist(1))  
    if(r<90)  
    {
      SubDist[NumSub]=r;
      SubList[NumSub]=j;
      NumSub++;
      if(NumSub>SubCatHost.GrpLen_Sub[HostList[0]]*NumHost)
      {
	printf("error: too many subhaloes to hold\n");exit(1);
      }
    }
  }
}
    
printf("Found %d hosts, %d subs\n", (int)NumHost, (int)NumSub);fflush(stdout);
track_last_infall(NsnapHost, SubList, NumSub);  
track_first_infall(NsnapHost, SubList, NumSub);

sprintf(buf,"%s/anal/MassStrip_Chunyan3_%d.Mhost%d",SUBCAT_DIR,NsnapHost,MassType);	
myfopen(fp,buf,"w");

fprintf(fp, "HostID, SubID, HostMass[1e10Msun/h], SubMass, Distance[kpc/h], MInf_Last, MInfVir_Last, ZInf_Last, MInf_First, MInfVir_First, ZInf_First, IsDirectInfall\n");
for(i=0;i<NumSub;i++)
{
  subid=SubList[i];
  hostid=SubCatHost.HaloChains[subid].HostID;
  hostsubid=SubCatHost.GrpOffset_Sub[hostid];
  if(!prosubFirst[i].FlagSp)
  fprintf(fp,"%d, %d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %d\n", hostid, subid, 
	Mhost[hostid], SubCatHost.SubLen[subid]*header.mass[1], SubDist[i],
	prosubLast[i].sublen*header.mass[1], prosubLast[i].Mvir, z[prosubLast[i].SnapInfall],
	prosubFirst[i].sublen*header.mass[1],prosubFirst[i].Mvir, z[prosubFirst[i].SnapInfall],
	prosubLast[i].DirectInfall);
}
fclose(fp);

free_sub_table(&SubCatHost);
myfree(Mhost);

  return 0;
}

