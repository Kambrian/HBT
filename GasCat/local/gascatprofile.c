#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include "../globals.c"
#include "../mymath.c"
#include "../load_group.c"
#include "../tree.c"


#define NBIN 10
#include "gas_stuff.c"

int main()
{
//~ char inputdir[512]="/SANdisk5/kambrain/Sim6501/SubCat";
//~ char fofdir[512]="/raid1/hywang/ReSim/SIM6501/group_catalogue"; //"/home/kambrain/fof_hy";
//~ char snapdir[512]="/SANdisk5/SIM6501";
//~ char outputdir[1024]="/SANdisk5/kambrain/Sim6501/SubCat/anal";
char subdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5";
char gasdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5/GasCat/nonthermal_exclusive";
char fofdir[512]="/raid1/hywang/ReSim/SIM6702/group_catalogue"; //"/home/kambrain/fof_hy";
char snapdir[512]="/SANdisk4/data/NewR/SIM6702";
char outputdir[1024]="/SANdisk5/kambrain/Sim6702/SubCat5/GasCat/nonthermal_exclusive";

CATALOGUE Cat;
SUBCATALOGUE SubCat;
GASCATALOGUE GasCat;
int Nsnap=99,i,j,k,pid,subid,grpid;

float *cen,L_div[NBIN],Rgho_div[NBIN],T_div[NBIN],r_div[NBIN];
int n_div[NBIN];

int subbox_grid[3][2],dr_bin;
float rscale,dr;


float lowerlim,upperlim,Lx;
double rfun(double,void *);
double integ,err;
gsl_function F;
size_t neval;

char buf[1024];
FILE *fpsubprof;
logfile=stdout;

F.function=&rfun;
	
	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_sub_catalogue(Nsnap,&SubCat,subdir);
	load_particle_data(Nsnap,&Cat,&SubCat,snapdir);
	load_gas_data(Nsnap,snapdir);
	Rtidal=mymalloc(sizeof(float)*SubCat.Nsubs);
	load_tidal_radius(Nsnap,Rtidal,SubCat.Nsubs,subdir);
	load_gas_cat(Nsnap,&GasCat,gasdir);
	
	//~ makell(Gdat.Pos,NP_gas);

	sprintf(buf,"%s/gasprof_%03d",outputdir,Nsnap);
	if((fpsubprof=fopen(buf,"w"))==NULL)
	{
		printf("error: file open failed for %s!\n",buf);
		exit(1);
	}

for(subid=0;subid<SubCat.Nsubs;subid++)
{
	for(i=0;i<10;i++)	
	{
		L_div[i]=0;
		T_div[i]=0;
		Rgho_div[i]=0;
		n_div[i]=0;
		r_div[i]=0;
	}
	if(SubCat.SubLen[subid])
	{
	cen=Pdat.Pos[SubCat.PSubArr[subid][0]];
	rscale=Rtidal[i]*FACT_Relax;
for(i=0;i<GasCat.SubLen[subid];i++)
{
					pid=GasCat.PSubArr[subid][i];
					dr=distance(Gdat.Pos[pid],cen);
					if(dr<rscale)
					{
					dr_bin=floor(dr/rscale*NBIN);
						if(Gdat.U[pid]<5e3)//[lower,upper]~[20,80],the integral is 10 order of magnitude lower than that in [1,4],so omit these non-X-ray gas (T<2e5K)(~20eV)
							Lx=0;
						else
						{
							lowerlim=0.5e6/3.48/Gdat.U[pid]; //0.5~2kev
							upperlim=2e6/3.48/Gdat.U[pid];
							gsl_integration_qng(&F,lowerlim,upperlim,0,0.001,&integ,&err,&neval);
							Lx=headerA.mass[0]*Gdat.Rho[pid]/(headerA.time*headerA.time*headerA.time)*sqrt(Gdat.U[pid])*integ*3040.7;
						}
						L_div[dr_bin]+=Lx;
						Rgho_div[dr_bin]+=Gdat.Rho[pid];
						T_div[dr_bin]+=Gdat.U[pid];
						r_div[dr_bin]+=dr;
						n_div[dr_bin]++;
					}
}
}
	for(i=0;i<NBIN;i++)
	{
		if(n_div[i])
		{
		L_div[i]/=n_div[i];
		Rgho_div[i]/=n_div[i];
		T_div[i]=40.38*T_div[i]/n_div[i];
		r_div[i]=r_div[i]/n_div[i]/rscale;
		}
		else
		{
			r_div[i]=(i+0.5)/NBIN;
		}
	}
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",L_div[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",Rgho_div[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",T_div[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%g\t",r_div[i]);
		for(i=0;i<NBIN;i++)	fprintf(fpsubprof,"%d\t",n_div[i]);
		fprintf(fpsubprof,"%d\t%d\t%d\n",GasCat.SubLen[subid],SubCat.SubLen[subid],subid);
}
		
fclose(fpsubprof);
return 0;
}

double rfun(double x,void * param)
{
	return 1.0/(pow(x,0.4)*exp(x));
}


//~ float rfun(float x)
//~ {
	//~ return 1.0/(pow(x,0.4)*exp(x));
//~ }
