//to search for overlapped subhalos,i.e,within 4 times smoothing scale
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

struct NgbList
{
	int *list;
	int len;
};
void makell_sub(int np,float **Pos);
void freell_sub();
void collect_subngb(int cenid,struct NgbList *Slist,float rmax,float **SubCen);

#define SUBLEN(i) (i<0?0:SubCat.SubLen[i])
#define SUBPOS(i) Pdat.Pos[SubCat.PSubArr[i][0]]
#define DIST(i,j) ((i<0||j<0||SubCat.SubLen[i]==0||SubCat.SubLen[j]==0)?-1*SofteningHalo:distance(SUBPOS(i),SUBPOS(j)))
int main(int argc, char** argv)
{
	SUBCATALOGUE SubCat;
	float **SubCen;
	int *SubList,(*SubPair)[2];
	int subid,i,j,jmax,Nalive,Npair;
	int Nsnap=0;
	struct NgbList Slist;
	FILE *fp;
	char buf[1024];

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	fresh_ID2Index(&SubCat,-2);
	SubCen=mymalloc(sizeof(float *)*SubCat.Nsubs);
	SubList=mymalloc(sizeof(int)*SubCat.Nsubs);
	SubPair=mymalloc(sizeof(int)*2*SubCat.Nsubs);
	for(i=0,j=0;i<SubCat.Nsubs;i++)  //make clean list of real subhalos (exclude sublen<100)
	{
		if(SubCat.SubLen[i]>100)
		{
			SubList[j]=i;
			SubCen[j]=Pdat.Pos[SubCat.PSubArr[i][0]];
			j++;
		}	
	}
	Nalive=j;Npair=0;jmax=0;
	makell_sub(Nalive,SubCen);
	for(i=0;i<Nalive;i++)
	{
		collect_subngb(i,&Slist,4.*SofteningHalo,SubCen);
		for(j=0;j<Slist.len;j++)
		{
			SubPair[Npair+j][0]=SubList[i];
			SubPair[Npair+j][1]=SubList[Slist.list[j]];
			if(jmax<Slist.len)
			jmax=Slist.len;//maximum number of neighbours for a sub
		}
		Npair+=Slist.len;
		myfree(Slist.list);
	}
	freell_sub();
	printf("Npair %d, <Npair> %g, Nmax %d\n",Npair,Npair*2.0/(float)Nalive,jmax);
	printf("ID1\t\tID2\t\tLen1\t\tLen2\t\tD/soft\n");
	for(i=0;i<Npair;i++)
	{
		printf("%d\t\t%d\t\t%d\t\t%d\t\t%g\n",
		SubPair[i][0],SubPair[i][1],
		SUBLEN(SubPair[i][0]),SUBLEN(SubPair[i][1]),
		distance(SUBPOS(SubPair[i][0]),SUBPOS(SubPair[i][1]))/SofteningHalo);
	}
	
	sprintf(buf,"%s/anal/trace_overlap_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"w");
	for(i=0;i<Npair;i++)//distance
	{
		fprintf(fp,"%g\t",DIST(SubPair[i][0],SubPair[i][1])/SofteningHalo);
	}
	for(i=0;i<Npair;i++)//big mass
	{
		fprintf(fp,"%d\t",SUBLEN(SubPair[i][0]));;
	}
	for(i=0;i<Npair;i++)//small
	{
		fprintf(fp,"%d\t",SUBLEN(SubPair[i][1]));;
	}
	fprintf(fp,"\n");
	
/*	while(Nsnap>0)
	{
	int Npro,Nsplitter,*sp2pro;
	load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);
	int flaglive;
	flaglive=0;
	for(i=0;i<Npair;i++)
		for(j=0;j<2;j++)
		{
			subid=SubPair[i][j];
			if(subid>=0)
			{
				subid=SubCat.HaloChains[subid].ProSubID;
				if(subid>=Npro)
				SubPair[i][j]=sp2pro[subid];
				else
				SubPair[i][j]=subid;
				flaglive=1;
			}
		}
	erase_sub_catalogue(&SubCat);
	free_sp2pro(sp2pro,Npro,Nsplitter);
	if(flaglive==0) break;
	
	Nsnap--;
	printf("Nsnap=%d\n",Nsnap);fflush(stdout);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	fresh_ID2Index(&SubCat,-2);
	
	for(i=0;i<Npair;i++)
	{
		printf("%d\t\t%d\t\t%d\t\t%d\t\t%g\n",
		SubPair[i][0],SubPair[i][1],
		SUBLEN(SubPair[i][0]),SUBLEN(SubPair[i][1]),
		DIST(SubPair[i][0],SubPair[i][1])/SofteningHalo);
	}
	
	for(i=0;i<Npair;i++)//distance
	{
		fprintf(fp,"%g\t",DIST(SubPair[i][0],SubPair[i][1])/SofteningHalo);
	}
	for(i=0;i<Npair;i++)//big mass
	{
		fprintf(fp,"%d\t",SUBLEN(SubPair[i][0]));;
	}
	for(i=0;i<Npair;i++)//small
	{
		fprintf(fp,"%d\t",SUBLEN(SubPair[i][1]));;
	}
	fprintf(fp,"\n");
	
	}
*/	
	
	while(Nsnap<MaxSnap-1)  //downward trace
	{
	int Nsubs,*pro2dest;
	load_pro2dest(Nsnap,&pro2dest,&Nsubs,SUBCAT_DIR);
	int flaglive;
	flaglive=0;
	for(i=0;i<Npair;i++)
		for(j=0;j<2;j++)
		{
			subid=SubPair[i][j];
			if(subid>=0)
			{
				subid=pro2dest[subid];
				SubPair[i][j]=subid;
				flaglive=1;
			}
		}
	erase_sub_catalogue(&SubCat);
	free_pro2dest(pro2dest);
	if(flaglive==0) break;
	
	Nsnap++;
	printf("Nsnap=%d\n",Nsnap);fflush(stdout);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	fresh_ID2Index(&SubCat,-2);
	
	for(i=0;i<Npair;i++)
	{
		printf("%d\t\t%d\t\t%d\t\t%d\t\t%g\n",
		SubPair[i][0],SubPair[i][1],
		SUBLEN(SubPair[i][0]),SUBLEN(SubPair[i][1]),
		DIST(SubPair[i][0],SubPair[i][1])/SofteningHalo);
	}
	
	for(i=0;i<Npair;i++)//distance
	{
		fprintf(fp,"%g\t",DIST(SubPair[i][0],SubPair[i][1])/SofteningHalo);
	}
	for(i=0;i<Npair;i++)//big mass
	{
		fprintf(fp,"%d\t",SUBLEN(SubPair[i][0]));;
	}
	for(i=0;i<Npair;i++)//small
	{
		fprintf(fp,"%d\t",SUBLEN(SubPair[i][1]));;
	}
	fprintf(fp,"\n");
	
	}
	
	fclose(fp);
	
	return 0;
}

#define NDIV 200
static int hoc[NDIV][NDIV][NDIV],*ll;
static float range[3][2], step[3];
void makell_sub(int np,float **Pos)
{
	int i,j,grid[3];
	
	#define POS(i,j) Pos[i][j]
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
	printf("finished linked list.\n");
}
void freell_sub()
{
	myfree(ll);
}
void collect_subngb(int cenid,struct NgbList *Slist,float rmax,float **SubCen)
{
	int i,j,k,pid,subbox_grid[3][2],maxlen,len;
	float *cen,rscale,dr;
	int *Src;
	
	len=0;
	maxlen=10;
	Src=mymalloc(sizeof(float)*maxlen);
	cen=SubCen[cenid];
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
				for(;pid>=0;pid=ll[pid]) //loop in this cell
				{
					if(pid<=cenid) continue;//only count pid>cenid to avoid double counting
					double dx,dy,dz;
					dx=SubCen[pid][0]-cen[0];
					if(dx<-rmax||dx>rmax) continue;
					dy=SubCen[pid][1]-cen[1];
					if(dy<-rmax||dy>rmax) continue;
					dz=SubCen[pid][2]-cen[2];
					if(dz<-rmax||dz>rmax) continue;
					dr=sqrt(dx*dx+dy*dy+dz*dz);
					if(dr<rmax)
					{
						Src[len]=pid;
						len++;
						if(len>=maxlen)
						{
							maxlen*=2;
							Src=realloc(Src,sizeof(int)*maxlen);
						}
					}
				}
			}
			Src=realloc(Src,sizeof(int)*len);
			
	Slist->len=len;
	Slist->list=Src;
}
