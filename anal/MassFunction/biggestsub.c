//extract biggest subhalos
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define RMIN 0
#define RMAX 1
#define VirType 0

struct MassList
{
	float *list;
	int Len;
};
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
		fread(Nvir,sizeof(int),3,fp);
		for(j=0;j<3;j++)
		Mvir[i][j]=Nvir[j]*partmass;
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}
void makell_sub();
void collect_submass(int grpid,struct MassList *Mlist,float rmin,float rmax);
void freell_sub();

SUBCATALOGUE SubCat;
double partmass;
float (*Mvir)[3],(*Rvir)[3];

int main(int argc, char** argv)
{
	
	int Nsubs,*pro2dest,Npro,Nsplitter,*sp2pro;
	int i,Nsnap=99;
	int grpid,subid,pid;
	struct MassList GrpList;
	
	FILE *fp;
	char buf[1024];
	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	partmass=header.mass[1];
	Mvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
	load_halo_virial_size(Mvir,Rvir,(float)partmass,SubCat.Ngroups,Nsnap);
	
	makell_sub();
	
	sprintf(buf,"%s/anal/massfun/biggestsub_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"w");
	fwrite(&SubCat.Nsubs,sizeof(int),1,fp);
	for(i=0;i<SubCat.Ngroups;i++)
	{
		float m;
		subid=SubCat.GrpOffset_Sub[i];
		m=SubCat.SubLen[subid]*partmass;
		fwrite(&m,sizeof(float),1,fp);
		collect_submass(i,&GrpList,RMIN*Rvir[i][VirType],RMAX*Rvir[i][VirType]);
		if(GrpList.Len)
		m=GrpList.list[Fmax_of_vec(GrpList.list,GrpList.Len)];
		else
		m=0;
		myfree(GrpList.list);
		fwrite(&m,sizeof(float),1,fp);
		fwrite(&Mvir[i][0],sizeof(float),1,fp);
	}
	fwrite(&SubCat.Nsubs,sizeof(int),1,fp);
	fclose(fp);
	freell_sub();
	erase_sub_catalogue(&SubCat);
	myfree(Mvir);
	myfree(Rvir);
	return 0;
}
#define NDIV 200
static int hoc[NDIV][NDIV][NDIV],*ll;
static float range[3][2], step[3];
void makell_sub()
{
	int i,j,grid[3],np;
	
	#define POS(i,j) SubCat.Property[i].CoM[j]
	np=SubCat.Nsubs;
	printf("creating linked list..\n");
	ll=mymalloc(sizeof(int)*np);
	/*determining enclosing cube*/
	for(i=0;i<3;i++)
		for(j=0;j<2;j++)
			range[i][j]=POS(0,i);
	for(i=1;i<np;i++)
		for(j=0;j<3;j++)
		{
			if(POS(i,j)<range[j][0])
				range[j][0]=POS(i,j);
			else if(POS(i,j)>range[j][1])
				range[j][1]=POS(i,j);
		}
	for(j=0;j<3;j++)
		step[j]=(range[j][1]-range[j][0])/NDIV;
	
	/*initialize hoc*/
	int *phoc=&(hoc[0][0][0]);
	for(i=0;i<NDIV*NDIV*NDIV;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((POS(i,j)-range[j][0])/step[j]);
			if(grid[j]<0) 
				grid[j]=0;
			else if(grid[j]>=NDIV)
				grid[j]=NDIV-1;
		}
		ll[i]=hoc[grid[0]][grid[1]][grid[2]];
		hoc[grid[0]][grid[1]][grid[2]]=i; /*use hoc(floor(xsat)+1) as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
	#undef POS(i,j)
}
void freell_sub()
{
	myfree(ll);
}
void collect_submass(int grpid,struct MassList *Mlist,float rmin,float rmax)
{
	int i,j,k,cenid,pid,subbox_grid[3][2],maxlen,len;
	float *cen,rscale,dr,*Src;
	
	len=0;
	cenid=SubCat.GrpOffset_Sub[grpid];
	if(SubCat.SubLen[cenid])
	{
	maxlen=SubCat.GrpLen_Sub[grpid]*2;
	maxlen=((maxlen<SubCat.Nsubs)?maxlen:SubCat.Nsubs);
	Src=mymalloc(sizeof(float)*maxlen);
	cen=SubCat.Property[cenid].CoM;
	rscale=rmax;
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rscale-range[i][0])/step[i]);
	if(subbox_grid[i][0]<0)subbox_grid[i][0]=0;
	subbox_grid[i][1]=floor((cen[i]+rscale-range[i][0])/step[i]);
	if(subbox_grid[i][1]>=NDIV)subbox_grid[i][1]=NDIV-1;
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=hoc[i][j][k];
				while(pid>=0)
				{
					dr=distance(SubCat.Property[pid].CoM,cen);
					if(pid!=cenid&&SubCat.SubLen[pid]&&dr<rmax&&dr>rmin)
					{
						#ifdef NORM
						Src[len]=SubCat.SubLen[pid]*partmass/Mvir[grpid][VirType];
						#else
						Src[len]=SubCat.SubLen[pid]*partmass;
						#endif
						len++;
						if(len>=maxlen)
						{
							maxlen*=2;
							Src=realloc(Src,sizeof(float)*maxlen);
						}
					}
					pid=ll[pid];
				}
			}
			Src=realloc(Src,sizeof(float)*len);
	}
	else
	{
		Src=NULL;
	}
	Mlist->Len=len;
	Mlist->list=Src;
}
