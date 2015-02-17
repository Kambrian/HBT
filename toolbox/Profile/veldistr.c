#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

 int comp_erg(const void *a, const void *b)//used to sort energy in ascending order; note that the most bound will come first( energy <0)
{
  if(((struct Energy *) a)->Erg > ((struct Energy *) b)->Erg)
    return +1;

  if(((struct Energy *) a)->Erg < ((struct Energy *) b)->Erg)
    return -1;

  return 0;
}

int main(int argc,char **argv)
{
	float E,Rvir,bin1[2]={0.4,0.5},bin2[2]={0.9,1.0},dr,E1,E2;
	struct Energy *SubEsort;
	int i,CoreLen,HaloLen,*HaloArr,M1=0,M2=0,Nsnap=99,grpid=0;
	double Hz,sqa,Time,PartMass;
	double vx,vy,vz,sx,sy,sz,dvx,dvy,dvz,dx,dy,dz,*pot;
	CATALOGUE Cat;
	
	char fofdir[512]=GRPCAT_DIR; 
	char snapdir[512]=SNAPSHOT_DIR;
	char outputdir[1024];
	FILE *fp1,*fp2;
	char buf[1024];
	logfile=stdout;
	if(argc!=3){printf("usage: %s [Nsnap] [grpid]\n",argv[0]);exit(1);}
	Nsnap=atoi(argv[1]);
	grpid=atoi(argv[2]);
	sprintf(outputdir,"%s/anal",SUBCAT_DIR);
	sprintf(buf,"%s/veldistr.Sn%02d.G%02d.bin%02d_%02d",outputdir,Nsnap,grpid,(int)(10*bin1[0]),(int)(10*bin1[1]));
	myfopen(fp1,buf,"w");
	sprintf(buf,"%s/veldistr.Sn%02d.G%02d.bin%02d_%02d",outputdir,Nsnap,grpid,(int)(10*bin2[0]),(int)(10*bin2[1]));
	myfopen(fp2,buf,"w");
	
	load_group_catalogue(Nsnap,&Cat,fofdir);
	load_particle_data(Nsnap,snapdir);fresh_ID2Index(&Cat,-1); 
	HaloLen=Cat.Len[grpid];
	HaloArr=Cat.PIDorIndex+Cat.Offset[grpid];
	Time=headerA.time; 
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#end
	Hz=headerA_Hz;	PartMass=headerA.mass[1];
	Rvir=pow(G*HaloLen*PartMass/100/Hz/Hz,1.0/3); 
		
	pot=mymalloc(sizeof(double)*HaloLen);
	tree_tree_allocate(TREE_ALLOC_FACTOR*HaloLen,HaloLen);
	maketree(HaloLen,HaloArr);
	vx=vy=vz=sx=sy=sz=0.;
	SubEsort=mymalloc(sizeof(struct Energy)*HaloLen);	
	
		#pragma omp parallel if(HaloLen>NParaMin)
		{//\\start para
			#pragma omp for //schedule(dynamic)
			for(i=0;i<HaloLen;i++)
			{
				pot[i]=tree_treeevaluate_potential(HaloArr[i],HaloArr);
				pot[i]=(PartMass*pot[i]-PartMass/SofteningHalo)*G/Time;
			}
			#pragma omp for //schedule(dynamic)
			for(i=0;i<HaloLen;i++)
			{
			SubEsort[i].PID=HaloArr[i];
			SubEsort[i].Erg=pot[i];
			}
			#pragma omp single
			{
			qsort(SubEsort,HaloLen,sizeof(struct Energy),comp_erg);
			CoreLen=(int)(HaloLen*CoreFrac0);
			CoreLen=(CoreLen>CoreLenMin)?CoreLen:CoreLenMin;
			}
			#pragma omp for reduction (+:vx,vy,vz,sx,sy,sz) //schedule(dynamic)
			for(i=0;i<CoreLen;i++)
			{
			vx+=Pdat.Vel[SubEsort[i].PID][0];
			vy+=Pdat.Vel[SubEsort[i].PID][1];
			vz+=Pdat.Vel[SubEsort[i].PID][2];
			sx+=Pdat.Pos[SubEsort[i].PID][0];
			sy+=Pdat.Pos[SubEsort[i].PID][1];
			sz+=Pdat.Pos[SubEsort[i].PID][2];
			}
		}//\end para
			vx/=CoreLen;
			vy/=CoreLen;
			vz/=CoreLen;
			sx/=CoreLen;
			sy/=CoreLen;
			sz/=CoreLen;
			free(SubEsort);
			for(i=0;i<HaloLen;i++)
			{
			 dvx=sqa*(Pdat.Vel[HaloArr[i]][0]-vx);//relative vel.
			 dvy=sqa*(Pdat.Vel[HaloArr[i]][1]-vy);
			 dvz=sqa*(Pdat.Vel[HaloArr[i]][2]-vz);
			 dx=Time*(Pdat.Pos[HaloArr[i]][0]-sx);
			 dy=Time*(Pdat.Pos[HaloArr[i]][1]-sy);
			 dz=Time*(Pdat.Pos[HaloArr[i]][2]-sz);
			 dvx+=Hz*dx;//add Hubble flow
			 dvy+=Hz*dy;
			 dvz+=Hz*dz;
			 dr=sqrt(dx*dx+dy*dy+dz*dz)/Rvir;
			 E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
			 if(dr<bin1[0])
			 M1+=1;
			 else if(dr>bin1[0]&&dr<bin1[1])
			 fprintf(fp1,"%g,%g,%g,%g,%g\n",dvx,dvy,dvz,E,dr);
			 if(dr<bin2[0])
			 M2+=1;
			 else if(dr>bin2[0]&&dr<bin2[1])
			 fprintf(fp2,"%g,%g,%g,%g,%g\n",dvx,dvy,dvz,E,dr);
			}
			E1=G*PartMass*M1/(Rvir*bin1[0]);//circular energy
			E2=G*PartMass*M2/(Rvir*bin2[0]);
			fprintf(fp1,"%g,%d,%g,%g,%g\n",Hz*bin1[0]*Rvir,M1,Rvir,E1,(bin1[0]+bin1[1])*0.5);
			fprintf(fp2,"%g,%d,%g,%g,%g\n",Hz*bin2[0]*Rvir,M2,Rvir,E2,(bin2[0]+bin2[1])*0.5);
			fclose(fp1);fclose(fp2);
			return 0;
}
