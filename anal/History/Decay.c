//decay time, the time when the current mass has decayed to only a small portion of the infall mass
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define RUN_AqA4  //need special match of snapnum, since AqA4 and AqA4_RelaxInf use different number of snaps

#ifdef RUN_AqA4
#define SUBCAT_DIR_FREE "/home/jvbq85/data/HBT/data/AqA4/subcatRelaxInf"
#define MaxSnapFree 128
HBTInt snaplistfree[MaxSnapFree]={7,15,23,31,39,47,55,63,71,79,87,95,103,111,119,127,
	135,143,151,159,167,175,183,191,199,207,215,223,231,239,247,255,263,
	271,279,287,295,303,311,319,327,335,343,351,359,367,375,383,391,399,
	407,415,423,431,439,447,455,463,471,479,487,495,503,511,519,527,535,
	543,551,559,567,575,583,591,599,607,615,623,631,639,647,655,663,671,
	679,687,695,703,711,719,727,735,743,751,759,767,775,783,791,799,807,
	815,823,831,839,847,855,863,871,879,887,895,903,911,919,927,935,943,
	951,959,967,975,983,991,999,1007,1015,1023};
#else
#define SUBCAT_DIR_FREE "/home/jvbq85/data/HBT/data/AqA3/subcatRelaxInf"
#define MaxSnapFree MaxSnap	
#endif
	
#define NDECAY 4

void load_sub_table_free(const HBTInt Nsnap,SUBCATALOGUE *SubCatFree,char * dir);
void free_sub_table_free(const HBTInt Nsnap,SUBCATALOGUE *SubCatFree);
int main(int argc,char **argv)
{
	SUBCATALOGUE SubCatBnd, SubCatFree;
	HBTReal DecayBarrier[NDECAY]={0.5, 0.1, 0.01, 0.001};
	struct
	{
		HBTInt subid;
		HBTInt grpidInfall;
		HBTInt sublen;
		HBTInt foflen;
		HBTInt SnapInfall;
		HBTInt FlagSp;
		HBTInt SnapDecay[NDECAY]; //decay to DecayBarrier times the infall fof mass
	} *prosub;
	HBTInt i,j,subid,fofid,Nsnap,Nsub,trash;
	HBTInt *pro2dest,Npro;	
	
	FILE *fp;
	char buf[1024],cmd[1024];	
	logfile=stdout;
		
	sprintf(buf, "%s/anal/steller/SnapInfall_first_%03d",SUBCAT_DIR_FREE,MaxSnapFree-1);
	Nsub=count_lines(buf)-1;
	printf(""HBTIFMT" subs\n",Nsub);
	myfopen(fp,buf,"r");
	fscanf(fp,"LsubInfl,LfofInfl,SnapInfl,GrpIdInfl,Subidborn,Lsub0,FlagSp,Ncross\n");
	prosub=mymalloc(sizeof(*prosub)*Nsub);
	i=0;
	while(1)
	{
		HBTInt nread;
		nread=fscanf(fp,""HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT"\n",&prosub[i].sublen,&prosub[i].foflen,&prosub[i].SnapInfall,&prosub[i].grpidInfall,&trash, &trash, &prosub[i].FlagSp, &trash);
		if(EOF==nread) break;
		#ifdef RUN_AqA4
		prosub[i].SnapInfall=snaplistfree[prosub[i].SnapInfall];
		#endif
		for(j=0;j<NDECAY;j++) prosub[i].SnapDecay[j]=-1;
		prosub[i].subid=-1;
		i++;
	}
	//~ printf(""HBTIFMT","HBTIFMT"\n", i, prosub[i-1].sublen);
	fclose(fp);
	if(i!=Nsub)
	{
		printf("error reading %s, nlines mismatch,"HBTIFMT", "HBTIFMT"\n", buf, Nsub, i);
		exit(1);
	}
	
	printf("0000");
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
		printf("\b\b\b\b%04d",(int)Nsnap);fflush(stdout);
		load_sub_table(Nsnap,&SubCatBnd,SUBCAT_DIR);
		//~ load_sub_table_free(Nsnap,&SubCatFree,SUBCAT_DIR_FREE);
		load_pro2dest(Nsnap-1,&pro2dest,&Npro,SUBCAT_DIR);
		for(i=0;i<Nsub;i++)
		{
			if(Nsnap==prosub[i].SnapInfall) //at infall
			{
				fofid=prosub[i].grpidInfall;
				prosub[i].subid=SubCatBnd.GrpOffset_Sub[fofid];
				for(j=0;j<NDECAY;j++)
				{
					//~ if(prosub[i].SnapDecay[j]<0)//not found yet
						if(SubCatBnd.SubLen[prosub[i].subid]<prosub[i].foflen*DecayBarrier[j])//has decayed
							prosub[i].SnapDecay[j]=Nsnap;
				}
			}
			if(Nsnap>prosub[i].SnapInfall) //after infall
			{
				prosub[i].subid=pro2dest[prosub[i].subid];//update to current subid
				for(j=0;j<NDECAY;j++)
				{
					if(prosub[i].SnapDecay[j]<0)//not found yet
					{
						if(prosub[i].subid<0||SubCatBnd.SubLen[prosub[i].subid]<prosub[i].foflen*DecayBarrier[j])//has decayed
							prosub[i].SnapDecay[j]=Nsnap;
					}
				}
			}
		}
		free_sub_table(&SubCatBnd);
		//~ free_sub_table_free(Nsnap,&SubCatFree);
		free_pro2dest(pro2dest);
	}
	printf("\n");
	

	sprintf(buf, "%s/anal/steller/Decay_%03d",SUBCAT_DIR_FREE,MaxSnapFree-1);
	myfopen(fp,buf,"w");
	fprintf(fp,"#%f,%f,%f,%f\n",DecayBarrier[0],DecayBarrier[1],DecayBarrier[2],DecayBarrier[3]);
	for(i=0;i<Nsub;i++)
		fprintf(fp,HBTIFMT","HBTIFMT","HBTIFMT","HBTIFMT"\n",prosub[i].SnapDecay[0],prosub[i].SnapDecay[1],prosub[i].SnapDecay[2],prosub[i].SnapDecay[3]);
	fclose(fp);
	
	double z[MaxSnap];
	for(i=0;i<MaxSnap;i++)
	{
		load_particle_header(i,SNAPSHOT_DIR);
		z[i]=1./header.time-1;
	}
	load_sub_table(MaxSnap-1,&SubCatBnd,SUBCAT_DIR);
	
	sprintf(buf, "%s/anal/steller/Decay_%03d.redshift",SUBCAT_DIR_FREE,MaxSnapFree-1);
	myfopen(fp,buf,"w");
	fprintf(fp,"#flagsp,m0,m,1,%g,%g,%g,%g\n",DecayBarrier[0],DecayBarrier[1],DecayBarrier[2],DecayBarrier[3]);
	#define REDSHIFT(x) ((x)<0?-1:z[(x)])
	#define MASS(x) ((HBTInt)((x)<0?0:SubCatBnd.SubLen[(x)]))
	for(i=0;i<Nsub;i++)
		fprintf(fp,HBTIFMT","HBTIFMT","HBTIFMT",%.2f,%.2f,%.2f,%.2f,%.2f\n",
		         prosub[i].FlagSp,prosub[i].foflen,MASS(prosub[i].subid),
		        REDSHIFT(prosub[i].SnapInfall),REDSHIFT(prosub[i].SnapDecay[0]),
				REDSHIFT(prosub[i].SnapDecay[1]),REDSHIFT(prosub[i].SnapDecay[2]),
				REDSHIFT(prosub[i].SnapDecay[3]));
	fclose(fp);
	return 0;
}

void load_sub_table_free(const HBTInt Nsnap,SUBCATALOGUE *SubCatFree,char * dir)
{
	#ifdef RUN_AqA4
	HBTInt i;
	
	for(i=0;i<MaxSnapFree;i++)
	{
		if(Nsnap==snaplistfree[i])
		{
			load_sub_table(i,SubCatFree,dir);
			break;
		}
	}
	#else
	load_sub_table(Nsnap,SubCatFree,dir);
	#endif	
}
void free_sub_table_free(const HBTInt Nsnap,SUBCATALOGUE *SubCatFree)
{
	#ifdef RUN_AqA4
	HBTInt i;
	
	for(i=0;i<MaxSnapFree;i++)
	{
		if(Nsnap==snaplistfree[i])
		{
			free_sub_table(SubCatFree);
			break;
		}
	}	
	#else
	free_sub_table(SubCatFree);
	#endif
}
