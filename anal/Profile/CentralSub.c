//To find the neareast subhalo to halo center
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define CEN_MINPOT 0
#define CEN_MAXDEN 1
#define CEN_DMCOM 2
#define CEN_GALCOM 3

#define CENTYPE CEN_GALCOM

SUBCATALOGUE SubCat;
CATALOGUE Cat;
int *IDrhomax,*IDpotmin;
float (*Mvir)[3], (*Rvir)[3];

void load_haloCenID(int *IDrhomax, int *IDpotmin, int Nsnap,int Ngroups);
void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);

int main(int argc,char **argv)
{

	float *DCen;  //distance of assumed central to the center
	float *DMin;  //distance of neareast sub to center
	int *RankMin; //rank of nearest sub
	int Nsnap,i,j,pid,grpid;
	FILE *fp;
	char buf[1024];
	char outputdir[1024];
	
	sprintf(outputdir,"%s/anal/",SUBCAT_DIR);
	mkdir(outputdir,0755);
	logfile=stdout;
	
	if(argc!=2)
	{
		printf("To find the nearest subhalo at halo center.\n usage: %s [Nsnap] \n", argv[0]);
		exit(1);
	}
	Nsnap=atoi(argv[1]);
	
		sprintf(buf,"%s/anal/CentralSub_%03d.CEN%d",SUBCAT_DIR,Nsnap,CENTYPE);
		myfopen(fp,buf,"w");
		
		load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
		load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
		load_particle_data(Nsnap,SNAPSHOT_DIR);
		fresh_ID2Index(&Cat,FRSH_GRPCAT);
		fresh_ID2Index(&SubCat,FRSH_SUBCAT);
		IDrhomax=mymalloc(sizeof(int)*SubCat.Ngroups);
		IDpotmin=mymalloc(sizeof(int)*SubCat.Ngroups);
		#if CENTYPE == CEN_MINPOT || CENTYPE == CEN_MAXDEN
		load_haloCenID(IDrhomax, IDpotmin,Nsnap,SubCat.Ngroups);
		fresh_ID2Index(IDpotmin,SubCat.Ngroups); 
		fresh_ID2Index(IDrhomax,SubCat.Ngroups); 
		#endif
		Mvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
		Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
		load_halo_virial_size(Mvir,Rvir,header.mass[1],SubCat.Ngroups,Nsnap);
		DCen=mymalloc(sizeof(float)*SubCat.Ngroups);
		DMin=mymalloc(sizeof(float)*SubCat.Ngroups);
		RankMin=mymalloc(sizeof(int)*SubCat.Ngroups);
		for(grpid=0;grpid<SubCat.Ngroups;grpid++)
		{
			float d,dmin;
			int imin,subid;
			float cen[3];
			double mx,my,mz,dx,dy,dz,sx,sy,sz;
			
			sx=sy=sz=0.;
			switch(CENTYPE)
			{
				case CEN_GALCOM:
				#ifdef PERIODIC_BDR
				mx=SubCat.Property[SubCat.GrpOffset_Sub[grpid]].CoM[0];
				my=SubCat.Property[SubCat.GrpOffset_Sub[grpid]].CoM[1];
				mz=SubCat.Property[SubCat.GrpOffset_Sub[grpid]].CoM[2];
				#endif
				for(i=1;i<SubCat.GrpLen_Sub[grpid];i++)
				{
					subid=SubCat.GrpOffset_Sub[grpid]+i;
					#ifdef PERIODIC_BDR
					dx=SubCat.Property[subid].CoM[0]-mx;
					dy=SubCat.Property[subid].CoM[1]-my;
					dz=SubCat.Property[subid].CoM[2]-mz;
					sx+=NEAREST(dx);
					sy+=NEAREST(dy);
					sz+=NEAREST(dz);
					#else
					sx+=SubCat.Property[subid].CoM[0];
					sy+=SubCat.Property[subid].CoM[1];
					sz+=SubCat.Property[subid].CoM[2];
					#endif
				}
				sx/=(SubCat.GrpLen_Sub[grpid]-1);
				sy/=(SubCat.GrpLen_Sub[grpid]-1);
				sz/=(SubCat.GrpLen_Sub[grpid]-1);
				#ifdef PERIODIC_BDR
				sx+=mx;
				sy+=my;
				sz+=mz;
				#endif
				cen[0]=sx;cen[1]=sy;cen[2]=sz;
				break;
				case CEN_DMCOM:
				#ifdef PERIODIC_BDR
				mx=Pdat.Pos[Cat.PIDorIndex[Cat.Offset[grpid]]][0];
				my=Pdat.Pos[Cat.PIDorIndex[Cat.Offset[grpid]]][1];
				mz=Pdat.Pos[Cat.PIDorIndex[Cat.Offset[grpid]]][2];
				#endif
				//~ #pragma omp parallel for private(i,dx,dy,dz) reduction (+:sx,sy,sz)
				for(i=Cat.Offset[grpid];i<Cat.Offset[grpid]+Cat.Len[grpid];i++)
				{
					#ifdef PERIODIC_BDR
					dx=Pdat.Pos[Cat.PIDorIndex[i]][0]-mx;
					dy=Pdat.Pos[Cat.PIDorIndex[i]][1]-my;
					dz=Pdat.Pos[Cat.PIDorIndex[i]][2]-mz;
					sx+=NEAREST(dx);
					sy+=NEAREST(dy);
					sz+=NEAREST(dz);
					#else
					sx+=Pdat.Pos[Cat.PIDorIndex[i]][0];
					sy+=Pdat.Pos[Cat.PIDorIndex[i]][1];
					sz+=Pdat.Pos[Cat.PIDorIndex[i]][2];
					#endif
				}
				sx/=Cat.Len[grpid];
				sy/=Cat.Len[grpid];
				sz/=Cat.Len[grpid];
				#ifdef PERIODIC_BDR
				sx+=mx;
				sy+=my;
				sz+=mz;
				#endif
				cen[0]=sx;cen[1]=sy;cen[2]=sz;
				break;
				case CEN_MINPOT:
				for(j=0;j<3;j++)
				cen[j]=Pdat.Pos[IDpotmin[grpid]][j];
				break;
				case CEN_MAXDEN:
				for(j=0;j<3;j++)
				cen[j]=Pdat.Pos[IDrhomax[grpid]][j];
				break;
				default:
				printf("error: wrong centype\n");
				exit(1);
			}
			
			subid=SubCat.GrpOffset_Sub[grpid];
			if(SubCat.SubLen[subid])
				DCen[grpid]=distance(Pdat.Pos[SubCat.PSubArr[subid][0]],cen);
			else
				DCen[grpid]=BOXSIZE;
				
			dmin=DCen[grpid];imin=0;	
			for(i=1;i<SubCat.GrpLen_Sub[grpid];i++)
			{
				subid=SubCat.GrpOffset_Sub[grpid]+i;
				d=distance(Pdat.Pos[SubCat.PSubArr[subid][0]],cen);
				if(dmin>d)
				{
				dmin=d;
				imin=i;
				}
			}
			DMin[grpid]=dmin;
			RankMin[grpid]=imin;
		}
		
		
		fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
		fwrite(DCen,sizeof(float),SubCat.Ngroups,fp);
		fwrite(DMin,sizeof(float),SubCat.Ngroups,fp);
		for(i=0;i<SubCat.Ngroups;i++)
		fwrite(&Mvir[i][0],sizeof(float),1,fp);
		for(i=0;i<SubCat.Ngroups;i++)
		fwrite(&Rvir[i][0],sizeof(float),1,fp);
		fwrite(RankMin,sizeof(int),SubCat.Ngroups,fp);
		fwrite(SubCat.GrpLen_Sub,sizeof(int),SubCat.Ngroups,fp);
		fwrite(&SubCat.Ngroups,sizeof(int),1,fp);
		fclose(fp);
		
		myfree(RankMin);
		myfree(DMin);
		myfree(DCen);
		myfree(Mvir);
		myfree(Rvir);
		myfree(IDpotmin);
		myfree(IDrhomax);
		erase_sub_catalogue(&SubCat);
		free_catalogue(&Cat);

	
return 0;
}

void load_haloCenID(int *IDrhomax, int *IDpotmin, int Nsnap,int Ngroups)
{  //note these loaded are IDs, need to fresh to index before use.
	int N;
	char buf[1024];
	FILE *fp;
	
	if(Ngroups)
	if(!IDrhomax||!IDpotmin)
	{
		printf("error: allocate IDrhomax and IDpotmin first!\n");
		exit(1);
	}
	sprintf(buf,"%s/profile/haloCenID_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	fread(&N,sizeof(int),1,fp);
	if(N!=Ngroups)
	{
		printf("error loading %s\n 1: %d!=%d\n", buf, N, Ngroups);
		exit(1);
	}
	fread(IDrhomax,sizeof(int),Ngroups,fp);
	fread(IDpotmin,sizeof(int),Ngroups,fp);
	fread(&N,sizeof(int),1,fp);
	if(N!=Ngroups)
	{
		printf("error loading %s\n 2: %d!=%d\n", buf, N, Ngroups);
		exit(2);
	}
	fclose(fp);
}
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
