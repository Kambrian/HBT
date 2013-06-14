// test how much local matter can be accreted onto fof's adaptive remnant 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "binding.c"

void makell_DM();
void collect_dm_particles(int subid,SUBCATALOGUE*SubCat,int *P2DMSrcLen,int **P2DMSrc);

int main(int nargc,char **argv)
{
int SnapLoad,SnapBind;
CATALOGUE Cat;
SUBCATALOGUE SubCat;
int *DMLoc,DMLocLen,BndLocLen,GrpLen,*GrpPIDs,Lmain_removed,*PInd_main_removed,GrpLen2,*GrpPIDs2,flag_Grp2=0;
int fofid,subid,i;
int *Pro2dest,Npro;
struct SubProperty Prop;
FILE *fp;
char buf[1024];

float s[3],vs[3];
float CoreFrac, MassRelax_Factor;

logfile=stdout;	
if(nargc!=4)
{
	printf("usage: %s [SnapLoad] [fofid] [MassRelax]\n",argv[0]);
	exit(1);
}
SnapLoad=atoi(argv[1]);
fofid=atoi(argv[2]);
MassRelax_Factor=atof(argv[3]);

CoreFrac=1.0/MassRelax_Factor/MassRelax_Factor;

sprintf(buf, "%s/anal/follow/follow_%03d_%d.adpt%2.1f.loc",SUBCAT_DIR,SnapLoad,fofid,MassRelax_Factor);

myfopen(fp,buf,"w");

		load_group_catalogue(SnapLoad,&Cat,GRPCAT_DIR);
		/*==back up the fof==*/
		GrpLen=Cat.Len[fofid];
		GrpPIDs=mymalloc(sizeof(int)*GrpLen);
		memcpy(GrpPIDs,Cat.PIDorIndex+Cat.Offset[fofid],sizeof(int)*GrpLen);
		#ifdef GRPINPUT_INDEX  //convert index to id
		load_particle_data(SnapLoad,SNAPSHOT_DIR);
		for(i=0;i<GrpLen;i++)
		#ifdef PID_ORDERED
			GrpPIDs[i]++;
		#else
			GrpPIDs[i]=Pdat.PID[GrpPIDs[i]];
		#endif
		#endif
		/*==initialize GrpPIDs2 ==*/
		GrpLen2=GrpLen;
		GrpPIDs2=mymalloc(sizeof(int)*GrpLen2);
		memcpy(GrpPIDs2,GrpPIDs,sizeof(int)*GrpLen2);
		/*==compact Cat to contain just that fof==*/
		free_catalogue(&Cat);
		Cat.Nids=GrpLen;
		Cat.PIDorIndex=mymalloc(sizeof(int)*GrpLen);
		memcpy(Cat.PIDorIndex,GrpPIDs,sizeof(int)*GrpLen);
			
		for(SnapBind=SnapLoad;SnapBind<MaxSnap;SnapBind++)
		{
			printf("Snap %d\n",SnapBind);fflush(stdout);
			load_particle_data(SnapBind,SNAPSHOT_DIR);
			fresh_ID2Index(Cat.PIDorIndex,Cat.Nids); 	
			unbind_core_ext(&(Cat.Nids),&(Cat.PIDorIndex),s,&Lmain_removed,&PInd_main_removed,CoreFrac,vs);	
			//~ printf("CoreFrac: %f\n",CoreFrac);
			free(PInd_main_removed);
			if(flag_Grp2)  
			{
				if(Cat.Nids>GrpLen2) //update GrpPIDs2 when subhalo has grown big enough
				{
				GrpLen2=Cat.Nids;
				free(GrpPIDs2);
				GrpPIDs2=mymalloc(sizeof(int)*GrpLen2);
				for(i=0;i<GrpLen2;i++)
				#ifdef PID_ORDERED
					GrpPIDs2[i]=Cat.PIDorIndex[i]+1;
				#else
					GrpPIDs2[i]=Pdat.PID[Cat.PIDorIndex[i]];
				#endif
				CoreFrac=((float)GrpLen2)/((float)GrpLen)/MassRelax_Factor;//?? necessary to change this?
		 		}	
				else if(Cat.Nids<GrpLen*CoreFrac)//update GrpPIDs with GrpPIDs2, and update GrpPIDs2 with the current sub (or its nearest progenitor?)
				{
				GrpLen=GrpLen2;
				free(GrpPIDs);
				GrpPIDs=GrpPIDs2;
				GrpLen2=Cat.Nids;
				GrpPIDs2=mymalloc(sizeof(int)*GrpLen2);
				for(i=0;i<GrpLen2;i++)
				#ifdef PID_ORDERED
					GrpPIDs2[i]=Cat.PIDorIndex[i]+1;
				#else
					GrpPIDs2[i]=Pdat.PID[Cat.PIDorIndex[i]];
				#endif
				CoreFrac=((float)GrpLen2)/((float)GrpLen)/MassRelax_Factor;/*every time when GrpLen2 changes, 
																		* update CoreFrac to ensure we will be 
																		* operating in the same range ,
																		* that is, the lower limit sub will be unbind() with
																		* near real CoreFrac, rather than a possibly much higher
																		* CoreFrac if not updated*/
				}
			}
			else if(Cat.Nids<=GrpLen*CoreFrac*MassRelax_Factor)// register GrpPIDs2 when subhalo has reduced enough, and adapt CoreFrac
			{
				GrpLen2=Cat.Nids;
				free(GrpPIDs2);
				GrpPIDs2=mymalloc(sizeof(int)*GrpLen2);
				for(i=0;i<GrpLen2;i++)
				#ifdef PID_ORDERED
					GrpPIDs2[i]=Cat.PIDorIndex[i]+1;
				#else
					GrpPIDs2[i]=Pdat.PID[Cat.PIDorIndex[i]];
				#endif
				flag_Grp2=1;
				CoreFrac=((float)GrpLen2)/((float)GrpLen)/MassRelax_Factor;//?? necessary to change this?
			}
					
			/*==copy bound mass to subcat==*/
			SubCat.Ngroups=1;SubCat.GrpOffset_Sub=NULL;SubCat.GrpLen_Sub=NULL;
			SubCat.Nsubs=1;	create_sub_cat(&SubCat);
			SubCat.SubLen[0]=Cat.Nids;
			SubCat.PSubArr[0]=mymalloc(sizeof(int)*Cat.Nids);
			for(i=0;i<3;i++)
			{
			SubCat.Property[0].CoM[i]=s[i];
			SubCat.Property[0].VCoM[i]=vs[i];
			}
			memcpy(SubCat.PSubArr[0],Cat.PIDorIndex,sizeof(int)*Cat.Nids);
			subid=0;
			/*==local cat==*/
			makell_DM();
			collect_dm_particles(subid,&SubCat,&DMLocLen,&DMLoc);
			BndLocLen=DMLocLen;
			unbind_add(&BndLocLen,&DMLoc,&Prop,SubCat.SubLen[subid],SubCat.PSubArr[subid],s,vs);	
			
			fprintf(fp,"%d,%d,%d,%d,%f\n",SnapBind,BndLocLen+SubCat.SubLen[subid],DMLocLen+SubCat.SubLen[subid],SubCat.SubLen[subid],CoreFrac);fflush(fp);

			myfree(DMLoc);
			erase_sub_catalogue(&SubCat);
			Cat.Nids=GrpLen;
			free(Cat.PIDorIndex);
			Cat.PIDorIndex=mymalloc(sizeof(int)*GrpLen);
			memcpy(Cat.PIDorIndex,GrpPIDs,sizeof(int)*GrpLen);
		}
		fclose(fp);	
		return 0;
}

#define DMCollect_ScaleRelax 2
#define NDIV 200
static int hoc_DM[NDIV][NDIV][NDIV],ll_DM[NP_DM];
static float range_DM[3][2],step_DM[3];
void makell_DM()
{
	float (*pos)[3];
	int np;
	int i,j,grid[3];
	pos=Pdat.Pos;
	np=NP_DM;
	printf("creating linked list..\n");
	
	/*determining enclosing cube*/
	for(i=0;i<3;i++)
		for(j=0;j<2;j++)
			range_DM[i][j]=pos[0][i];
	for(i=1;i<np;i++)
		for(j=0;j<3;j++)
		{
			if(pos[i][j]<range_DM[j][0])
				range_DM[j][0]=pos[i][j];
			else if(pos[i][j]>range_DM[j][1])
				range_DM[j][1]=pos[i][j];
		}
	for(j=0;j<3;j++)
		step_DM[j]=(range_DM[j][1]-range_DM[j][0])/NDIV;
	
	/*initialize hoc*/
	int *phoc=&(hoc_DM[0][0][0]);
	for(i=0;i<NDIV*NDIV*NDIV;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((pos[i][j]-range_DM[j][0])/step_DM[j]);
			if(grid[j]<0) 
				grid[j]=0;
			else if(grid[j]>=NDIV)
				grid[j]=NDIV-1;
		}
		ll_DM[i]=hoc_DM[grid[0]][grid[1]][grid[2]];
		hoc_DM[grid[0]][grid[1]][grid[2]]=i; /*use hoc(floor(xsat)+1) as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
}


void collect_dm_particles(int subid,SUBCATALOGUE*SubCat,int *P2DMSrcLen,int **P2DMSrc)
{
	int DMlen,i,j,k,pid,subbox_grid[3][2],maxlen,*DMSrc;
	float *cen,rvir,rscale,dr;
	short *SubMask;
	
	SubMask=malloc(sizeof(short)*NP_DM);
	for(i=0;i<NP_DM;i++)
		SubMask[i]=1;
	for(i=0;i<SubCat->SubLen[subid];i++)
	{
		pid=SubCat->PSubArr[subid][i];
		SubMask[pid]=0;
	}	
	
	DMlen=0;
	if(SubCat->SubLen[subid])
	{
	maxlen=SubCat->SubLen[subid]*2;
	maxlen=((maxlen<NP_DM)?maxlen:NP_DM);
	DMSrc=mymalloc(sizeof(int)*maxlen);
	cen=SubCat->Property[subid].CoM;
	rvir=comoving_virial_radius(SubCat->SubLen[subid]);
	rscale=rvir*DMCollect_ScaleRelax;
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rscale-range_DM[i][0])/step_DM[i]);
	if(subbox_grid[i][0]<0)subbox_grid[i][0]=0;
	subbox_grid[i][1]=floor((cen[i]+rscale-range_DM[i][0])/step_DM[i]);
	if(subbox_grid[i][1]>=NDIV)subbox_grid[i][1]=NDIV-1;
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=hoc_DM[i][j][k];
				while(pid>=0)
				{
					dr=distance(Pdat.Pos[pid],cen);
					if(dr<rscale&&SubMask[pid])//within rscale and not belong to sub
					{
						DMSrc[DMlen]=pid;
						DMlen++;
						if(DMlen>=maxlen)
						{
							maxlen*=2;
							DMSrc=realloc(DMSrc,sizeof(int)*maxlen);
						}
					}
					pid=ll_DM[pid];
				}
			}
			DMSrc=realloc(DMSrc,sizeof(int)*DMlen);
	}
	else
	{
		DMSrc=NULL;
	}
	*P2DMSrcLen=DMlen;
	*P2DMSrc=DMSrc;
	free(SubMask);
}
