#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

//#define SUBFIND_DIR "/home/kambrain/data/6702DM/subcatS"   //define this properly to apply on subfind, 
															//undef for HBT
															
#define NDIV 256
#define NBIN 120
#define Factor_relax 3
LINKLIST ll;
#ifdef PERIODIC_BDR
#define MAKELL(pos,np,ndiv) make_linklist_box(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_shift_gridid(x,NDIV)
#else
#define MAKELL(pos,np,ndiv) make_linklist(&ll,pos,np,ndiv)
#define FIXGRID(x) linklist_round_gridid(x,NDIV)
#endif

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue
#undef SUBCAT_DIR
#define SUBCAT_DIR SUBFIND_DIR
#endif

int main()
{
	char inputdir[512]=SUBCAT_DIR;
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	
CATALOGUE Cat;
SUBCATALOGUE SubCat;
int Nsnap=99,i,j,k,pid,subid,hostsubid;

int Ngas,Ndm,Nother;
float *cen,rvir,rcen,r_this[NBIN],r_other[NBIN],r_back[NBIN];
int n_this[NBIN],n_other[NBIN],n_back[NBIN];
int *PartHost;
int subbox_grid[3][2],dr_bin;
float dr;

char buf[1024];
FILE *fpsubprof;
logfile=stdout;
sprintf(outputdir,"%s/anal",SUBCAT_DIR);	

	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	load_particle_data(Nsnap,snapdir);
	fill_PIDHash();
	fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCat,-2);
	free_PIDHash();
	Ngas=header.npart[0];	Ndm=header.npart[1];	Nother=header.npart[2]+header.npart[3]+header.npart[4]+header.npart[5];
	PartHost=mymalloc(sizeof(int)*Ndm);
	//~ PartHost--;//so that it's accessed through [1,Ndm]
	for(i=0;i<Ndm;i++)
		PartHost[i]=-1;
for(subid=0;subid<SubCat.Nsubs;subid++)
{
	for(i=0;i<SubCat.SubLen[subid];i++)
		PartHost[SubCat.PSubArr[subid][i]]=subid;
}
	sprintf(buf,"%s/subprof_dm_%03d.%01dvir",outputdir,Nsnap,Factor_relax);
	myfopen(fpsubprof,buf,"w");

	MAKELL(Pdat.Pos,NP_DM,NDIV);
	//printf("%d,%d\n",hoc[10][9][133],hoc[77][2][1]);
	printf("finished linkedlist\n");

for(subid=0;subid<SubCat.Nsubs;subid++)
{
	for(i=0;i<NBIN;i++)	
	{
		n_this[i]=0;
		n_other[i]=0;
		n_back[i]=0;
		r_this[i]=0;
		r_other[i]=0;
		r_back[i]=0;
	}
	if(SubCat.SubLen[subid])
	{
	cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
	rvir=pow(G*SubCat.SubLen[subid]*header.mass[1]/100/header.Hz/header.Hz,1.0/3)/header.time;
	rvir*=Factor_relax;	
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rvir-ll.range[i][0])/ll.step[i]);
	subbox_grid[i][1]=floor((cen[i]+rvir-ll.range[i][0])/ll.step[i]);
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=linklist_get_hoc(&ll,FIXGRID(i),FIXGRID(j),FIXGRID(k));
				while(pid>=0)
				{
					dr=distance(Pdat.Pos[pid],cen);
					if(dr<rvir)
					{
						dr_bin=floor(dr/rvir*NBIN);
						if(PartHost[pid]==subid)
						{
							r_this[dr_bin]+=dr;
							n_this[dr_bin]++;
						}
						else if(PartHost[pid]==-1)
						{
							r_back[dr_bin]+=dr;
							n_back[dr_bin]++;
						}
						else
						{
							r_other[dr_bin]+=dr;
							n_other[dr_bin]++;
						}
					}
					pid=ll.list[pid];
				}
			}
	for(i=0;i<NBIN;i++)
	{
		if(n_this[i])
		{
		r_this[i]=r_this[i]/n_this[i]/rvir;
		}
		else
		{
			r_this[i]=(i+0.5)/NBIN;
		}
		if(n_back[i])
		{
		r_back[i]=r_back[i]/n_back[i]/rvir;
		}
		else
		{
			r_back[i]=(i+0.5)/NBIN;
		}
		if(n_other[i])
		{
		r_other[i]=r_other[i]/n_other[i]/rvir;
		}
		else
		{
			r_other[i]=(i+0.5)/NBIN;
		}
	}
}
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",r_this[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%d\t",n_this[i]);//number of particles that belong to this sub
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",r_other[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%d\t",n_other[i]);//other subs
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",r_back[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%d\t",n_back[i]);//backgroud particles
		hostsubid=(SubCat.HaloChains[subid].HostID<0)?(subid):(SubCat.GrpOffset_Sub[SubCat.HaloChains[subid].HostID]);
		rcen=(SubCat.SubLen[subid]>0)?distance(cen,Pdat.Pos[SubCat.PSubArr[hostsubid][0]]):0;
		fprintf(fpsubprof,"%g\t%g\t%d\t%d\t%d\n",rvir/Factor_relax,rcen,SubCat.SubLen[subid],subid,SubCat.SubRank[subid]);
		//~ fflush(fpsubprof);
		//~ printf("%d\n",subid);
}

		
fclose(fpsubprof);
return 0;
}
