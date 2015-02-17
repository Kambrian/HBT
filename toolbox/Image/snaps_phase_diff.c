//to make snapshots of subhalo particles in phase space: (v,r)
//output them into three parts: SubFind particles, HBT-only particles, and background.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define DMCollect_ScaleRelax 0.5;

#define SUBFIND_DIR "/home/kambrain/data/6702DM/subcatS"

#ifdef SUBFIND_DIR
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
//~ #define load_sub_catalogue load_subfind_catalogue
//~ #undef SUBCAT_DIR
//~ #define SUBCAT_DIR SUBFIND_DIR
#endif

void makell_DM();
void collect_dm_particles(int subid,SUBCATALOGUE*SubCat,int *P2DMSrcLen,int **P2DMSrc,short *SubMask);

int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat,SubCatS;
	int Nsnap,i,j,pid;
	int subid;
	FILE *fp;
	char buf[1024];

	int *PIndexBT,NBT,NS,NBk;
	int *PIndexS,*PIndexBk;
	float sqa,tmp;
	char outputdir[1024];
	sprintf(outputdir,"%s/anal/image",SUBCAT_DIR);
	mkdir(outputdir,0755);
	logfile=stdout;
	
	Nsnap=atoi(argv[1]);
	subid=atoi(argv[2]);
	
//prepare particles
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_subfind_catalogue(Nsnap,&SubCatS,SUBFIND_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fresh_ID2Index(&SubCat,FRSH_SUBCAT);
	fresh_ID2Index(&SubCatS,FRSH_SUBCAT);
	
	PIndexS=SubCatS.PSubArr[subid];
	NBT=SubCat.SubLen[subid];
	NS=SubCatS.SubLen[subid];
	
	short *Mask;
	Mask=mymalloc(sizeof(short)*NP_DM);
	for(i=0;i<NP_DM;i++)
	Mask[i]=0;
	
	for(i=0;i<NS;i++)
	Mask[PIndexS[i]]=1;
	
	PIndexBT=mymalloc(sizeof(int)*NBT);
	for(i=0,j=0;i<NBT;i++)
	{
		pid=SubCat.PSubArr[subid][i];
		if(!Mask[pid])	
		{
			PIndexBT[j]=pid;
			j++;
			Mask[pid]=1;
		}
	}
	NBT=j;
	PIndexBT=realloc(PIndexBT,sizeof(int)*NBT);
	
	makell_DM();
	collect_dm_particles(subid,&SubCat,&NBk,&PIndexBk,Mask);
	printf("NSF=%d\nNBT=%d\nNBK=%d\n",NS,NBT,NBk);
	
		#ifdef VEL_INPUT_PHYSICAL
		sqa=1.0;
		#else
		 sqa = sqrt(header.time);
		 #endif
		 	
		sprintf(buf,"%s/posmap_%03d_%d.SF",outputdir,Nsnap,subid);
		myfopen(fp,buf,"w");
		fwrite(SubCatS.Property[subid].CoM,sizeof(float),3,fp);
		for(i=0;i<NS;i++)
			fwrite(Pdat.Pos[PIndexS[i]],sizeof(float),3,fp);
		fclose(fp);
		
		sprintf(buf,"%s/velmap_%03d_%d.SF",outputdir,Nsnap,subid);
		myfopen(fp,buf,"w");
		fwrite(SubCatS.Property[subid].VCoM,sizeof(float),3,fp);
		for(i=0;i<NS;i++)
		{
			for(j=0;j<3;j++)
			{
			tmp=Pdat.Vel[PIndexS[i]][j]*sqa;
			fwrite(&tmp,sizeof(float),1,fp);
			}
		}
		fclose(fp);

		sprintf(buf,"%s/posmap_%03d_%d.BT",outputdir,Nsnap,subid);
		myfopen(fp,buf,"w");
		fwrite(SubCat.Property[subid].CoM,sizeof(float),3,fp);
		for(i=0;i<NBT;i++)
			fwrite(Pdat.Pos[PIndexBT[i]],sizeof(float),3,fp);
		fclose(fp);
		
		sprintf(buf,"%s/velmap_%03d_%d.BT",outputdir,Nsnap,subid);
		myfopen(fp,buf,"w");
		fwrite(SubCat.Property[subid].VCoM,sizeof(float),3,fp);
		for(i=0;i<NBT;i++)
		{
			for(j=0;j<3;j++)
			{
			tmp=Pdat.Vel[PIndexBT[i]][j]*sqa;
			fwrite(&tmp,sizeof(float),1,fp);
			}
		}
		fclose(fp);
		
		sprintf(buf,"%s/posmap_%03d_%d.BK",outputdir,Nsnap,subid);
		myfopen(fp,buf,"w");
		fwrite(SubCat.Property[subid].CoM,sizeof(float),3,fp);
		for(i=0;i<NBk;i++)
			fwrite(Pdat.Pos[PIndexBk[i]],sizeof(float),3,fp);
		fclose(fp);
		
		sprintf(buf,"%s/velmap_%03d_%d.BK",outputdir,Nsnap,subid);
		myfopen(fp,buf,"w");
		fwrite(SubCat.Property[subid].VCoM,sizeof(float),3,fp);
		for(i=0;i<NBk;i++)
		{
			for(j=0;j<3;j++)
			{
			tmp=Pdat.Vel[PIndexBk[i]][j]*sqa;
			fwrite(&tmp,sizeof(float),1,fp);
			}
		}
		fclose(fp);
		
		free_sub_catalogue(&SubCat);
				
return 0;
}


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

void collect_dm_particles(int subid,SUBCATALOGUE*SubCat,int *P2DMSrcLen,int **P2DMSrc,short *SubMask)
{
	int DMlen,i,j,k,pid,subbox_grid[3][2],maxlen,*DMSrc;
	float *cen,rvir,rscale,dr;
	
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
					if(dr<rscale&&!SubMask[pid])//within rscale and not belong to sub
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
}
