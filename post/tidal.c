#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include "globals.c"
#include "mymath.c"
#include "load_group.c"

#define NDIV 200
#define NBIN 50
#define Factor_relax 4
int hoc[NDIV][NDIV][NDIV],ll[NP_DM];
#include "linkedlist.c"


int main()
{
//~ char inputdir[512]="/SANdisk5/kambrain/Sim6501/SubCat";
//~ char fofdir[512]="/raid1/hywang/ReSim/SIM6501/group_catalogue"; //"/home/kambrain/fof_hy";
//~ char snapdir[512]="/SANdisk5/SIM6501";
//~ char outputdir[1024]="/SANdisk5/kambrain/Sim6501/SubCat/anal";
char inputdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5";
char fofdir[512]="/raid1/hywang/ReSim/SIM6702/group_catalogue"; //"/home/kambrain/fof_hy";
char snapdir[512]="/SANdisk4/data/NewR/SIM6702";
char outputdir[1024]="/SANdisk5/kambrain/Sim6702/SubCat5/anal";

CATALOGUE Cat;
SUBCATALOGUE SubCat;
int Nsnap=99,i,j,k,pid,subid,hostsubid;

int Ngas,Ndm,Nother;
float *cen,rvir,rcen,r_this[NBIN],r_other[NBIN],r_back[NBIN];
int n_this[NBIN],n_other[NBIN],n_back[NBIN];
int *PartHost;
int subbox_grid[3][2],dr_bin;
float range[3][2],step[3],dr;


char buf[1024];
FILE *fpsubprof;
logfile=stdout;

	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,inputdir);
	load_particle_data(Nsnap,&Cat,&SubCat,snapdir);
	Ngas=headerA.npart[0];	Ndm=headerA.npart[1];	Nother=headerA.npart[2]+headerA.npart[3]+headerA.npart[4]+headerA.npart[5];
	PartHost=mymalloc(sizeof(int)*Ndm);
	PartHost--;//so that it's accessed through [1,Ndm]
	for(i=0;i<Ndm;i++)
		PartHost[i]=-1;
for(subid=0;subid<SubCat.Nsubs;subid++)
{
	for(i=0;i<SubCat.SubLen[subid];i++)
		PartHost[SubCat.PSubArr[subid][i]]=subid;
}
	sprintf(buf,"%s/subprof_dm_%03d",outputdir,Nsnap);
	if((fpsubprof=fopen(buf,"w"))==NULL)
	{
		printf("error: file open failed for %s!\n",buf);
		exit(1);
	}

	makell(Pdat.Pos,Ndm,NDIV,range,step);
	printf("%d,%d,%g\n",hoc[10][9][133],hoc[77][2][1],Hubble);
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
	rvir=pow(G*SubCat.SubLen[subid]*headerA.mass[1]/100/headerA_Hz/headerA_Hz,1.0/3);
	rvir*=Factor_relax;	
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rvir-range[i][0])/step[i]);
	subbox_grid[i][1]=floor((cen[i]+rvir-range[i][0])/step[i]);
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=hoc[i][j][k];
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
					pid=ll[pid];
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
		hostsubid=(SubCat.HostID[subid]<0)?(subid):(SubCat.GrpOffset_Sub[SubCat.HostID[subid]]);
		rcen=(SubCat.SubLen[subid]>0)?distance(cen,Pdat.Pos[SubCat.PSubArr[hostsubid][0]]):0;
		fprintf(fpsubprof,"%g\t%g\t%d\t%d\t%d\n",rvir/Factor_relax,rcen,SubCat.SubLen[subid],subid,SubCat.SubRank[subid]);
		//~ fflush(fpsubprof);
		//~ printf("%d\n",subid);
}

		
fclose(fpsubprof);
return 0;
}
