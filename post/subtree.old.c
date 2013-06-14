#define OMP
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "globals.c"
#include "load_group.c"

	
struct idsort
	{ 
		int proID;
		int ind;
	};

int compare_id_pro(const void *a, const void *b)//used to sort pro id in ascending order
{
  if(((struct idsort *) a)->proID > ((struct idsort *) b)->proID)
    return +1;

  if(((struct idsort *) a)->proID < ((struct idsort *) b)->proID)
    return -1;

  return 0;
}

int main()
{
	FILE *fp;
	char buf[1024];
	char bufsp[1024]; FILE *fpsp; 
	int i,j,proID,desID,subid,haloid,prohalo,deshalo,proRank,desRank;
	int *pro2halo,*dest2halo, tmp[4],*sp2pro;
	int Nsnap=0;
 	struct idsort *proind;
	SUBCATALOGUE SubCatA,SubCatB;
	char subdir[1024]="/SANdisk5/kambrain/Sim6702/SubCat3";
	
	sprintf(buf, "%s/tree/subtree.dat",subdir);
	if(!(fp= fopen(buf, "w")))
	{
		printf("can't open file `%s'\n", buf);
		exit(1);
	}
	
	sprintf(bufsp, "%s/splitters.log",subdir);
	if(!(fpsp= fopen(bufsp, "r")))
	{
		printf("Error: can't open file `%s'\n", bufsp);
		exit(1);
	}
	fscanf(fpsp,"%d,%d,%d,%d\n",tmp,tmp+1,tmp+2,tmp+3);
	
	load_sub_catalogue(Nsnap,&SubCatA,subdir);
	load_sub_catalogue(Nsnap,&SubCatB,subdir);
	fprintf(fp,"ParentID\tParentHalo\tproinHalo\tDestID\tDestHalo\tdestinHalo\tDestWeight\tDestPortion\tParentSize\tDestSize\n");
	fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",-100,Nsnap,SubCatB.Ngroups,SubCatB.Nsubs,SubCatB.Nbirth,SubCatB.Ndeath,SubCatB.NQuasi,SubCatB.Nsplitter,-100,-100);//-100 indicates the generation division;
	pro2halo=malloc(1);
	dest2halo=mymalloc(sizeof(int)*SubCatB.Nsubs);
	//prepare sub2halo table
#ifdef OMP
	#pragma omp parallel 
	{
	#pragma omp for nowait private(haloid,i) schedule(dynamic)
#endif  
	for(haloid=0;haloid<SubCatB.Ngroups;haloid++)
	{
		for(i=0;i<SubCatB.GrpLen_Sub[haloid];i++)
		{
			dest2halo[SubCatB.GrpOffset_Sub[haloid]+i]=haloid;
		}
	}
	//don't forget the quasihalos
#ifdef OMP
	#pragma omp for
#endif 
	for(subid=SubCatB.Nsubs-SubCatB.NQuasi;subid<SubCatB.Nsubs;subid++)
		dest2halo[subid]=-1;//-1means background halo
#ifdef OMP
	}
#endif
	
	for(Nsnap=1;Nsnap<10;Nsnap++)
	{
		free_sub_catalogue(&SubCatA);
		memcpy(&SubCatA,&SubCatB,sizeof(SUBCATALOGUE));
		free(pro2halo);
		pro2halo=dest2halo;
		
		load_sub_catalogue(Nsnap,&SubCatB,subdir);
		dest2halo=mymalloc(sizeof(int)*SubCatB.Nsubs);
		proind=mymalloc(sizeof(struct idsort)*SubCatB.Nsubs);
		//prepare sub2halo table
#ifdef OMP
	#pragma omp parallel 
	{
	#pragma omp for nowait private(haloid,i,subid) schedule(dynamic)
#endif  
		for(haloid=0;haloid<SubCatB.Ngroups;haloid++)
		{
			for(i=0;i<SubCatB.GrpLen_Sub[haloid];i++)
			{
			subid=SubCatB.GrpOffset_Sub[haloid]+i;
			dest2halo[subid]=haloid;
			proind[subid].proID=SubCatB.ProSubID[subid];
			proind[subid].ind=subid;
			}
		}
		//don't forget the quasihalos
#ifdef OMP
	#pragma omp for
#endif 
		for(subid=SubCatB.Nsubs-SubCatB.NQuasi;subid<SubCatB.Nsubs;subid++)
		{
			dest2halo[subid]=-1;
			proind[subid].proID=SubCatB.ProSubID[subid];
			proind[subid].ind=subid;
		}
#ifdef OMP
	}
#endif
	
	
	//prepare splitter2sub table
	if(SubCatB.Nsplitter)
	{
		if(tmp[1]!=(Nsnap-1))
		{
			printf("error reading splitter info, snap=%d,tmp=%d,%d,%d,%d\n",Nsnap-1,tmp[0],tmp[1],tmp[2],tmp[3]);
			exit(1);
		}
		sp2pro=mymalloc(sizeof(int)*SubCatB.Nsplitter);
		fscanf(fpsp,"%d,%d,%d,%d\n",tmp,tmp+1,tmp+2,tmp+3);
		do
		{
			for(i=0;i<tmp[1]-1;i++)
			{
				sp2pro[tmp[2]+i]=tmp[0];
			}
			fscanf(fpsp,"%d,%d,%d,%d\n",tmp,tmp+1,tmp+2,tmp+3);		
		}while(tmp[0]>=0);
	}
		
	
	qsort(proind,SubCatB.Nsubs,sizeof(struct idsort),compare_id_pro);
		
		j=0;
		//skip baby
		while(proind[j].proID<0)
		{
			j++;
		}
		while(j<SubCatB.Nsubs-SubCatB.Nsplitter)
		{
			desID=proind[j].ind;
			deshalo=dest2halo[desID];
			desRank=((deshalo<0)?(0):(desID-SubCatB.GrpOffset_Sub[deshalo]));
			proID=proind[j].proID;
			prohalo=pro2halo[proID];
			proRank=((prohalo<0)?(0):(proID-SubCatA.GrpOffset_Sub[prohalo]));
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\n",proID,prohalo,proRank,desID,deshalo,desRank,SubCatB.HostWght[desID],SubCatB.HostPortion[desID],SubCatA.SubLen[proID],SubCatB.SubLen[desID]);
			j++;
		}
		if(SubCatB.Nsplitter)//splitters
		{
			sp2pro-=SubCatA.Nsubs;
			while(j<SubCatB.Nsubs)
			{
				desID=proind[j].ind;
				deshalo=dest2halo[desID];
				desRank=((deshalo<0)?(0):(desID-SubCatB.GrpOffset_Sub[deshalo]));
				if((proID=proind[j].proID)<SubCatA.Nsubs)
				{
				printf("error! splitter=%d, Nsubs=%d\n",proID,SubCatA.Nsubs);
				exit(2);
				}
				proID=sp2pro[proID];
				prohalo=pro2halo[proID];
				proRank=((prohalo<0)?(0):(proID-SubCatA.GrpOffset_Sub[prohalo]));
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%d\n",proID,prohalo,proRank,desID,deshalo,desRank,SubCatB.HostWght[desID],SubCatB.HostPortion[desID],SubCatA.SubLen[proID],SubCatB.SubLen[desID]);
				j++;
			}
			free(sp2pro+SubCatA.Nsubs);
		}
		free(proind);
		fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",-100,Nsnap,SubCatB.Ngroups,SubCatB.Nsubs,SubCatB.Nbirth,SubCatB.Ndeath,SubCatB.NQuasi,SubCatB.Nsplitter,-100,-100);//-100 indicates the generation division;
	}
		fclose(fp);
		fclose(fpsp);
	return 0;
}	
