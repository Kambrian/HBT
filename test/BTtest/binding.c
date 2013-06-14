//applicable to both resimulation and periodical cosmological simulation
//remember to undef HALO_PARA to use parallelization here
#ifdef CoreLenMin
#undef CoreLenMin
#endif
#define CoreLenMin 1

//#define OUTPUT_CENTERS_TAG    //to just output initial CoM and VCoM rather than do the unbinding
#ifdef OUTPUT_CENTERS_TAG
#define OUTPUT_CENTERS {printf("%g,%g,%g,%g,%g,%g\n",sx,sy,sz,sqa*vx,sqa*vy,sqa*vz);exit(0);}
#define OUTPUT_CENTERS_arr {printf("%g,%g,%g,%g,%g,%g\n",s[0],s[1],s[2],sqa*vx,sqa*vy,sqa*vz);exit(0);}
#else
#define OUTPUT_CENTERS
#define OUTPUT_CENTERS_arr
#endif

/*using CoM plus VCoM to calculate kinetic energy,equivalent to unbind_core() with CoreFrac=1 */

 int comp_erg(const void *a, const void *b)//used to sort energy in ascending order; note that the most bound will come first( energy <0)
{
  if(((struct Energy *) a)->Erg > ((struct Energy *) b)->Erg)
    return +1;

  if(((struct Energy *) a)->Erg < ((struct Energy *) b)->Erg)
    return -1;

  return 0;
}

int unbind_mean(int *P2Len,int **P2PIndex, float *CoM,int *P2Len_removed, int **P2PIndex_removed) /*P2Len=&Len, *P2PIndex=PIndex, 
																							*where PIndex[Len] is the array of size Len; 
																							* both will be updated after unbind
																							* *P2PIndex_removed=PIndex_removed need 
																							* not to be allocated as input*/
{
	int * TmpPIndex, *SubPIndex,*SubPIndex_removed,id_minpot;
	float * TmpBindingE, E;
	struct Energy *SubEsort;
	int i,j,Nbound,SubLen,SubLen_removed;
	double Hz,sqa,Time,PartMass;
	double vx,vy,vz,sx,sy,sz,dvx,dvy,dvz,dx,dy,dz, *pot;
	#ifdef PERIODIC_BDR
	double mx,my,mz;//position of minimum potential
	#endif
	//time_t time_start,time_end,t1,t2,t3,t4,tt1,tt2,tt3,tt4,ttt2,ttt3; //for program timing
	//~ FILE *fp;
	//fp=fopen("converg.dat","w");
	
	
	//time_start=time(NULL);
	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	 Hz=header.Hz;
	 PartMass=header.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

	Nbound=*P2Len;
	SubPIndex=*P2PIndex;
	SubLen_removed=0;
	SubPIndex_removed=mymalloc(sizeof(int)*Nbound);
	
	//fprintf(fp,"%d\n",Nbound);
	//if(Nbound>NParaMin) printf("%d\n",Nbound);
	SubLen=Nbound*2+1;//just to start the loop
	TmpBindingE=mymalloc(1);//just to start the loop
	while(Nbound<SubLen*PrecMass)
	{	
		//fprintf(logfile,"Nbound=%d\n",Nbound);fflush(stdout);
		SubLen=Nbound;
		Nbound=0;
		TmpPIndex=mymalloc(sizeof(int)*SubLen);//this has been taken over by SubPIndex
		free(TmpBindingE);TmpBindingE=mymalloc(sizeof(float)*SubLen);
		pot=mymalloc(sizeof(double)*SubLen);//this be freed at the end of each loop
		
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
		maketree(SubLen,SubPIndex,Pdat.Pos);
		vx=vy=vz=0.;sx=sy=sz=0.;
		#pragma omp parallel  if(SubLen>NParaMin) 
		{//\\start para
			#pragma omp for schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{			
				pot[i]=tree_treeevaluate_potential(Pdat.Pos[SubPIndex[i]],SubPIndex,Pdat.Pos);
			    pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;	
			}
			#ifdef PERIODIC_BDR
			#pragma omp single
			{
			id_minpot=SubPIndex[Dmin_of_vec(pot,SubLen)];
			mx=Pdat.Pos[id_minpot][0];
			my=Pdat.Pos[id_minpot][1];
			mz=Pdat.Pos[id_minpot][2];
			}
			#endif
			#pragma omp for private(i,dx,dy,dz) reduction (+:vx,vy,vz,sx,sy,sz) schedule(dynamic,1)
			for(i=0;i<SubLen;i++)
			{
			vx+=Pdat.Vel[SubPIndex[i]][0];
			vy+=Pdat.Vel[SubPIndex[i]][1];
			vz+=Pdat.Vel[SubPIndex[i]][2];
			#ifdef PERIODIC_BDR
			dx=Pdat.Pos[SubPIndex[i]][0]-mx;
			dy=Pdat.Pos[SubPIndex[i]][1]-my;
			dz=Pdat.Pos[SubPIndex[i]][2]-mz;
			sx+=NEAREST(dx);
			sy+=NEAREST(dy);
			sz+=NEAREST(dz);
			//~ #else
			//~ sx+=Pdat.Pos[SubPIndex[i]][0];
			//~ sy+=Pdat.Pos[SubPIndex[i]][1];
			//~ sz+=Pdat.Pos[SubPIndex[i]][2];
			#endif
			}
			#pragma omp single
			{
			vx/=SubLen;
			vy/=SubLen;
			vz/=SubLen;
			#ifdef PERIODIC_BDR
			sx/=SubLen;
			sy/=SubLen;
			sz/=SubLen;
			sx+=mx;
			sy+=my;
			sz+=mz;
			#else
			sx=Nodes_base->way.s[0];//Center of Mass
			sy=Nodes_base->way.s[1];//Center of Mass
			sz=Nodes_base->way.s[2];//Center of Mass
			#endif
			OUTPUT_CENTERS
			}
			#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,E) schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
			 dvx=sqa*(Pdat.Vel[SubPIndex[i]][0]-vx);//relative vel.
			 dvy=sqa*(Pdat.Vel[SubPIndex[i]][1]-vy);
			 dvz=sqa*(Pdat.Vel[SubPIndex[i]][2]-vz);
			 dx=Pdat.Pos[SubPIndex[i]][0]-sx;
			 dy=Pdat.Pos[SubPIndex[i]][1]-sy;
			 dz=Pdat.Pos[SubPIndex[i]][2]-sz;
			 #ifdef PERIODIC_BDR
			 dx=NEAREST(dx);
			 dy=NEAREST(dy);
			 dz=NEAREST(dz);
			 #endif
			 dx*=Time;dy*=Time;dz*=Time;
			 dvx+=Hz*dx;//add Hubble flow
			 dvy+=Hz*dy;
			 dvz+=Hz*dz;
			E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
				#pragma omp critical (unbinding)
				if(E<0)
				{
				TmpPIndex[Nbound]=SubPIndex[i];
				TmpBindingE[Nbound]=E;
				Nbound++;
				}
				else
				{
				SubPIndex_removed[SubLen_removed]=SubPIndex[i];
				SubLen_removed++;	
				}
			}
		}//\\end para
		//fprintf(fp,"%d\n",Nbound);
		
		free(pot);
		tree_tree_free();
		free(SubPIndex);//am i also freeing SubCatTmp->PSubArr[subhaloid]??............................Ok,fine.
		SubPIndex=TmpPIndex;//take over the bound part to begin a new loop 
		if(Nbound<NBOUNDMIN)
		{
			memcpy(SubPIndex_removed+SubLen_removed,TmpPIndex,Nbound*sizeof(int));
			*P2PIndex_removed=SubPIndex_removed;
			*P2Len_removed=SubLen_removed+Nbound;
			free(TmpPIndex);
			free(TmpBindingE);
			*P2Len=0;
			*P2PIndex=NULL;
			CoM[0]=CoM[1]=CoM[2]=0.;
			return 0;//no bound sub found, return 0;
		}			
	}
	//t1=time(NULL);
	
	SubLen=Nbound;//NOTE: this is only necessary when using tolerated (approximated) binding criteria
	SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
	for(i=0;i<SubLen;i++)
	{
		SubEsort[i].PID=SubPIndex[i];
		SubEsort[i].Erg=TmpBindingE[i];
	}
	qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
	for(i=0;i<SubLen;i++)
		SubPIndex[i]=SubEsort[i].PID;//copy the sorted index back

	free(TmpBindingE);
	free(SubEsort);
	*P2Len=SubLen;
	*P2PIndex=SubPIndex;
	CoM[0]=sx;
	CoM[1]=sy;
	CoM[2]=sz; 
	*P2Len_removed=SubLen_removed;
	*P2PIndex_removed=realloc(SubPIndex_removed,sizeof(int)*SubLen_removed);
	//time_end=time(NULL);
	//printf("%ld (%f),%ld (%f) sec\n",t1-time_start,(float)(t1-time_start)/(time_end-time_start),time_end-t1,(float)(time_end-t1)/(time_end-time_start));
	return 1;//found bound structure, return 1;
}

/*using Core-averaged CoM and VCoM for kinetic energy calc*/
int unbind_core(int *P2Len,int **P2PIndex, float *CoM,int *P2Len_removed, int **P2PIndex_removed,float CoreFrac) /*P2Len=&Len, *P2PIndex=PIndex, 
																							*where PIndex[Len] is the array of size Len; 
																							* both will be updated after unbind
																							* *P2PIndex_removed=PIndex_removed need 
																							* not to be allocated as input*/
{
	int * TmpPIndex, *SubPIndex,*SubPIndex_removed;
	float * TmpBindingE, E;
	struct Energy *SubEsort;
	int i,Nbound,SubLen,SubLen_removed,CoreLen;
	double Hz,sqa,Time,PartMass;
	double vx,vy,vz,sx,sy,sz,dvx,dvy,dvz,dx,dy,dz,*pot;
	#ifdef PERIODIC_BDR
	double mx,my,mz;//position of minimum potential
	#endif
	//time_t time_start,time_end,t1,t2,t3,t4,tt1,tt2,tt3,tt4,ttt2,ttt3; //for program timing
	//~ FILE *fp;
	//fp=fopen("converg.dat","w");
	
	
	//time_start=time(NULL);
	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	 Hz=header.Hz;
	 PartMass=header.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

	Nbound=*P2Len;
	SubPIndex=*P2PIndex;
	SubLen_removed=0;
	SubPIndex_removed=mymalloc(sizeof(int)*Nbound);
	
	#ifdef LOG_NBOUND
	fprintf(logfile,"%d\n",Nbound);
	#endif
	//if(Nbound>NParaMin) printf("%d\n",Nbound);
	SubLen=Nbound*2+1;//just to start the loop
	TmpBindingE=mymalloc(1);//just to start the loop
	#ifdef LOG_NBOUND
	while(Nbound<SubLen)
	#else
	while(Nbound<SubLen*PrecMass)
	#endif
	{	
		//~ fprintf(logfile,"Nbound=%d\n",Nbound);fflush(stdout);
		SubLen=Nbound;
		Nbound=0;
		TmpPIndex=mymalloc(sizeof(int)*SubLen);
		TmpBindingE=realloc(TmpBindingE,sizeof(float)*SubLen);
		pot=mymalloc(sizeof(double)*SubLen);
		
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
		maketree(SubLen,SubPIndex,Pdat.Pos);
		vx=vy=vz=sx=sy=sz=0.;
		#pragma omp parallel  if(SubLen>NParaMin) 
		{//\\start para
			#pragma omp for schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
				pot[i]=tree_treeevaluate_potential(Pdat.Pos[SubPIndex[i]],SubPIndex,Pdat.Pos);
				pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;
			}
			#pragma omp single
			{
			SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
			}
			#pragma omp for schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
			SubEsort[i].PID=SubPIndex[i];
			SubEsort[i].Erg=pot[i];
			}
			#pragma omp single
			{
			qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
			CoreLen=(int)(SubLen*CoreFrac);
			CoreLen=(CoreLen>CoreLenMin)?CoreLen:CoreLenMin;
			#ifdef PERIODIC_BDR
			mx=Pdat.Pos[SubEsort[0].PID][0];
			my=Pdat.Pos[SubEsort[0].PID][1];
			mz=Pdat.Pos[SubEsort[0].PID][2];
			#endif
			//~ printf("corelen %d\n",CoreLen);
			}
			#pragma omp for private(i,dx,dy,dz) reduction (+:vx,vy,vz,sx,sy,sz) schedule(dynamic)
			for(i=0;i<CoreLen;i++)
			{
			vx+=Pdat.Vel[SubEsort[i].PID][0];
			vy+=Pdat.Vel[SubEsort[i].PID][1];
			vz+=Pdat.Vel[SubEsort[i].PID][2];
			#ifdef PERIODIC_BDR
			dx=Pdat.Pos[SubEsort[i].PID][0]-mx;
			dy=Pdat.Pos[SubEsort[i].PID][1]-my;
			dz=Pdat.Pos[SubEsort[i].PID][2]-mz;
			sx+=NEAREST(dx);
			sy+=NEAREST(dy);
			sz+=NEAREST(dz);
			#else
			sx+=Pdat.Pos[SubEsort[i].PID][0];
			sy+=Pdat.Pos[SubEsort[i].PID][1];
			sz+=Pdat.Pos[SubEsort[i].PID][2];
			#endif
			}
			#pragma omp single
			{
			vx/=CoreLen;
			vy/=CoreLen;
			vz/=CoreLen;
			sx/=CoreLen;
			sy/=CoreLen;
			sz/=CoreLen;
			free(SubEsort);
			#ifdef PERIODIC_BDR
			sx+=mx;
			sy+=my;
			sz+=mz;
			#endif
			OUTPUT_CENTERS
			}
			#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,E) schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
			 dvx=sqa*(Pdat.Vel[SubPIndex[i]][0]-vx);//relative vel.
			 dvy=sqa*(Pdat.Vel[SubPIndex[i]][1]-vy);
			 dvz=sqa*(Pdat.Vel[SubPIndex[i]][2]-vz);
			 dx=Pdat.Pos[SubPIndex[i]][0]-sx;
			 dy=Pdat.Pos[SubPIndex[i]][1]-sy;
			 dz=Pdat.Pos[SubPIndex[i]][2]-sz;
			 #ifdef PERIODIC_BDR
			 dx=NEAREST(dx);
			 dy=NEAREST(dy);
			 dz=NEAREST(dz);
			 #endif
			 dx*=Time;dy*=Time;dz*=Time;
			 dvx+=Hz*dx;//add Hubble flow
			 dvy+=Hz*dy;
			 dvz+=Hz*dz;
			E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
				#pragma omp critical (unbinding)
				if(E<0)
				{
				TmpPIndex[Nbound]=SubPIndex[i];
				TmpBindingE[Nbound]=E;
				Nbound++;
				}
				else
				{
				SubPIndex_removed[SubLen_removed]=SubPIndex[i];
				SubLen_removed++;	
				}
			}
		}//\\end para
		//~ #undef OMP
	#ifdef LOG_NBOUND
	fprintf(logfile,"%d\n",Nbound);
	#endif
		free(pot);
		tree_tree_free();
		free(SubPIndex);//am i also freeing SubCatTmp->PSubArr[subhaloid]??............................Ok,fine.
		SubPIndex=TmpPIndex;//take over the bound part to begin a new loop 
		if(Nbound<NBOUNDMIN)
		{
			memcpy(SubPIndex_removed+SubLen_removed,TmpPIndex,Nbound*sizeof(int));
			*P2PIndex_removed=SubPIndex_removed;
			*P2Len_removed=SubLen_removed+Nbound;
			free(TmpPIndex);
			free(TmpBindingE);
			*P2Len=0;
			*P2PIndex=NULL;
			CoM[0]=CoM[1]=CoM[2]=0.;
			return 0;//no bound sub found, return 0;
		}			
	}
	//t1=time(NULL);
	
	SubLen=Nbound;//NOTE: this is only necessary when using tolerated (approximated) binding criteria
	SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
	for(i=0;i<SubLen;i++)
	{
		SubEsort[i].PID=SubPIndex[i];
		SubEsort[i].Erg=TmpBindingE[i];
	}
	qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
	for(i=0;i<SubLen;i++)
		SubPIndex[i]=SubEsort[i].PID;//copy the sorted index back

	free(TmpBindingE);
	free(SubEsort);
	*P2Len=SubLen;
	*P2PIndex=SubPIndex;
	CoM[0]=sx;
	CoM[1]=sy;
	CoM[2]=sz; 
	*P2Len_removed=SubLen_removed;
	*P2PIndex_removed=realloc(SubPIndex_removed,sizeof(int)*SubLen_removed);
	//time_end=time(NULL);
	//printf("%ld (%f),%ld (%f) sec\n",t1-time_start,(float)(t1-time_start)/(time_end-time_start),time_end-t1,(float)(time_end-t1)/(time_end-time_start));
	return 1;//found bound structure, return 1;
}

/*using Core-averaged CoM and VCoM for kinetic energy calc
 * also return VCoM than unbind_core*/
int unbind_core_ext(int *P2Len,int **P2PIndex, float *CoM,int *P2Len_removed, int **P2PIndex_removed,float CoreFrac,float *VCoM) /*P2Len=&Len, *P2PIndex=PIndex, 
																							*where PIndex[Len] is the array of size Len; 
																							* both will be updated after unbind
																							* *P2PIndex_removed=PIndex_removed need 
																							* not to be allocated as input*/
{
	int * TmpPIndex, *SubPIndex,*SubPIndex_removed;
	float * TmpBindingE, E;
	struct Energy *SubEsort;
	int i,Nbound,SubLen,SubLen_removed,CoreLen;
	double Hz,sqa,Time,PartMass;
	double vx,vy,vz,sx,sy,sz,dvx,dvy,dvz,dx,dy,dz,*pot;
	#ifdef PERIODIC_BDR
	double mx,my,mz;//position of minimum potential
	#endif
	//time_t time_start,time_end,t1,t2,t3,t4,tt1,tt2,tt3,tt4,ttt2,ttt3; //for program timing
	//~ FILE *fp;
	//fp=fopen("converg.dat","w");
	
	
	//time_start=time(NULL);
	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	 Hz=header.Hz;
	 PartMass=header.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

	Nbound=*P2Len;
	SubPIndex=*P2PIndex;
	SubLen_removed=0;
	SubPIndex_removed=mymalloc(sizeof(int)*Nbound);
	
	//fprintf(fp,"%d\n",Nbound);
	//if(Nbound>NParaMin) printf("%d\n",Nbound);
	SubLen=Nbound*2+1;//just to start the loop
	TmpBindingE=mymalloc(1);//just to start the loop
	while(Nbound<SubLen*PrecMass)
	{	
		//fprintf(logfile,"Nbound=%d\n",Nbound);fflush(stdout);
		SubLen=Nbound;
		Nbound=0;
		TmpPIndex=mymalloc(sizeof(int)*SubLen);
		TmpBindingE=realloc(TmpBindingE,sizeof(float)*SubLen);
		pot=mymalloc(sizeof(double)*SubLen);
		
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
		maketree(SubLen,SubPIndex,Pdat.Pos);
		vx=vy=vz=sx=sy=sz=0.;
		#pragma omp parallel  if(SubLen>NParaMin) 
		{//\\start para
			#pragma omp for schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
				pot[i]=tree_treeevaluate_potential(Pdat.Pos[SubPIndex[i]],SubPIndex,Pdat.Pos);
				pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;
			}
			#pragma omp single
			{
			SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
			}
			#pragma omp for schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
			SubEsort[i].PID=SubPIndex[i];
			SubEsort[i].Erg=pot[i];
			}
			#pragma omp single
			{
			qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
			CoreLen=(int)(SubLen*CoreFrac);
			CoreLen=(CoreLen>CoreLenMin)?CoreLen:CoreLenMin;
			#ifdef PERIODIC_BDR
			mx=Pdat.Pos[SubEsort[0].PID][0];
			my=Pdat.Pos[SubEsort[0].PID][1];
			mz=Pdat.Pos[SubEsort[0].PID][2];
			#endif
			//~ printf("corelen %d\n",CoreLen);
			}
			#pragma omp for private(i,dx,dy,dz) reduction (+:vx,vy,vz,sx,sy,sz) schedule(dynamic)
			for(i=0;i<CoreLen;i++)
			{
			vx+=Pdat.Vel[SubEsort[i].PID][0];
			vy+=Pdat.Vel[SubEsort[i].PID][1];
			vz+=Pdat.Vel[SubEsort[i].PID][2];
			#ifdef PERIODIC_BDR
			dx=Pdat.Pos[SubEsort[i].PID][0]-mx;
			dy=Pdat.Pos[SubEsort[i].PID][1]-my;
			dz=Pdat.Pos[SubEsort[i].PID][2]-mz;
			sx+=NEAREST(dx);
			sy+=NEAREST(dy);
			sz+=NEAREST(dz);
			#else
			sx+=Pdat.Pos[SubEsort[i].PID][0];
			sy+=Pdat.Pos[SubEsort[i].PID][1];
			sz+=Pdat.Pos[SubEsort[i].PID][2];
			#endif
			}
			#pragma omp single
			{
			vx/=CoreLen;
			vy/=CoreLen;
			vz/=CoreLen;
			sx/=CoreLen;
			sy/=CoreLen;
			sz/=CoreLen;
			free(SubEsort);
			#ifdef PERIODIC_BDR
			sx+=mx;
			sy+=my;
			sz+=mz;
			#endif
			}
			#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,E) schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
			 dvx=sqa*(Pdat.Vel[SubPIndex[i]][0]-vx);//relative vel.
			 dvy=sqa*(Pdat.Vel[SubPIndex[i]][1]-vy);
			 dvz=sqa*(Pdat.Vel[SubPIndex[i]][2]-vz);
			 dx=Pdat.Pos[SubPIndex[i]][0]-sx;
			 dy=Pdat.Pos[SubPIndex[i]][1]-sy;
			 dz=Pdat.Pos[SubPIndex[i]][2]-sz;
			 #ifdef PERIODIC_BDR
			 dx=NEAREST(dx);
			 dy=NEAREST(dy);
			 dz=NEAREST(dz);
			 #endif
			 dx*=Time;dy*=Time;dz*=Time;
			 dvx+=Hz*dx;//add Hubble flow
			 dvy+=Hz*dy;
			 dvz+=Hz*dz;
			E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
				#pragma omp critical (unbinding)
				if(E<0)
				{
				TmpPIndex[Nbound]=SubPIndex[i];
				TmpBindingE[Nbound]=E;
				Nbound++;
				}
				else
				{
				SubPIndex_removed[SubLen_removed]=SubPIndex[i];
				SubLen_removed++;	
				}
			}
		}//\\end para
		//~ #undef OMP
		//fprintf(fp,"%d\n",Nbound);
		free(pot);
		tree_tree_free();
		free(SubPIndex);//am i also freeing SubCatTmp->PSubArr[subhaloid]??............................Ok,fine.
		SubPIndex=TmpPIndex;//take over the bound part to begin a new loop 
		if(Nbound<NBOUNDMIN)
		{
			memcpy(SubPIndex_removed+SubLen_removed,TmpPIndex,Nbound*sizeof(int));
			*P2PIndex_removed=SubPIndex_removed;
			*P2Len_removed=SubLen_removed+Nbound;
			free(TmpPIndex);
			free(TmpBindingE);
			*P2Len=0;
			*P2PIndex=NULL;
			CoM[0]=CoM[1]=CoM[2]=0.;
			VCoM[0]=VCoM[1]=VCoM[2]=0.;
			return 0;//no bound sub found, return 0;
		}			
	}
	//t1=time(NULL);
	
	SubLen=Nbound;//NOTE: this is only necessary when using tolerated (approximated) binding criteria
	SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
	for(i=0;i<SubLen;i++)
	{
		SubEsort[i].PID=SubPIndex[i];
		SubEsort[i].Erg=TmpBindingE[i];
	}
	qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
	for(i=0;i<SubLen;i++)
		SubPIndex[i]=SubEsort[i].PID;//copy the sorted index back

	free(TmpBindingE);
	free(SubEsort);
	*P2Len=SubLen;
	*P2PIndex=SubPIndex;
	CoM[0]=sx;
	CoM[1]=sy;
	CoM[2]=sz; 
	VCoM[0]=sqa*vx;
	VCoM[1]=sqa*vy;
	VCoM[2]=sqa*vz;
	*P2Len_removed=SubLen_removed;
	*P2PIndex_removed=realloc(SubPIndex_removed,sizeof(int)*SubLen_removed);
	//time_end=time(NULL);
	//printf("%ld (%f),%ld (%f) sec\n",t1-time_start,(float)(t1-time_start)/(time_end-time_start),time_end-t1,(float)(time_end-t1)/(time_end-time_start));
	return 1;//found bound structure, return 1;
}

 /* Unbind using the particle with minpot as center for Hubble flow,
  * VCoM for residual vel
  */
 
int unbind_minpotH(int *P2Len,int **P2PIndex, float *CoM,int *P2Len_removed, int **P2PIndex_removed) /*P2Len=&Len, *P2PIndex=PIndex, 
																							*where PIndex[Len] is the array of size Len; 
																							* both will be updated after unbind
																							* *P2PIndex_removed=PIndex_removed need 
																							* not to be allocated as input*/
{
	int * TmpPIndex, *SubPIndex,*SubPIndex_removed;
	float * TmpBindingE, E;
	struct Energy *SubEsort;
	int i,Nbound,SubLen,SubLen_removed,id_minpot;
	double Hz,sqa,Time,PartMass;
	double vx,vy,vz,s[3],dvx,dvy,dvz,dx,dy,dz,*pot;
	//time_t time_start,time_end,t1,t2,t3,t4,tt1,tt2,tt3,tt4,ttt2,ttt3; //for program timing
	//~ FILE *fp;
	//fp=fopen("converg.dat","w");
	
	
	//time_start=time(NULL);
	Time=header.time;
	 #ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	 Hz=header.Hz;
	 PartMass=header.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

	Nbound=*P2Len;
	SubPIndex=*P2PIndex;
	SubLen_removed=0;
	SubPIndex_removed=mymalloc(sizeof(int)*Nbound);
	
	//fprintf(fp,"%d\n",Nbound);
	//if(Nbound>NParaMin) printf("%d\n",Nbound);
	SubLen=Nbound*2+1;//just to start the loop
	TmpBindingE=mymalloc(1);//just to start the loop
	while(Nbound<SubLen*PrecMass)
	{	
		//fprintf(logfile,"Nbound=%d\n",Nbound);fflush(stdout);
		SubLen=Nbound;
		Nbound=0;
		TmpPIndex=mymalloc(sizeof(int)*SubLen);
		TmpBindingE=realloc(TmpBindingE,sizeof(float)*SubLen);
		pot=mymalloc(sizeof(double)*SubLen);
		
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
		maketree(SubLen,SubPIndex,Pdat.Pos);
		vx=vy=vz=0.;
		#pragma omp parallel  if(SubLen>NParaMin) 
		{//\\start para
			#pragma omp for reduction (+:vx,vy,vz) schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
				vx+=Pdat.Vel[SubPIndex[i]][0];
				vy+=Pdat.Vel[SubPIndex[i]][1];
				vz+=Pdat.Vel[SubPIndex[i]][2];
				pot[i]=tree_treeevaluate_potential(Pdat.Pos[SubPIndex[i]],SubPIndex,Pdat.Pos);
			    pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;	
			 }
			#pragma omp single
			{
			vx/=SubLen;
			vy/=SubLen;
			vz/=SubLen;
			id_minpot=SubPIndex[Dmin_of_vec(pot,SubLen)];
			s[0]=Pdat.Pos[id_minpot][0];
			s[1]=Pdat.Pos[id_minpot][1];
			s[2]=Pdat.Pos[id_minpot][2];
			//~ printf("MinPot:%g,%g,%g\n",s[0],s[1],s[2]);
			OUTPUT_CENTERS_arr
			}
			#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,E) schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
			 dvx=sqa*(Pdat.Vel[SubPIndex[i]][0]-vx);//relative vel.
			 dvy=sqa*(Pdat.Vel[SubPIndex[i]][1]-vy);
			 dvz=sqa*(Pdat.Vel[SubPIndex[i]][2]-vz);
			 dx=Pdat.Pos[SubPIndex[i]][0]-s[0];
			 dy=Pdat.Pos[SubPIndex[i]][1]-s[1];
			 dz=Pdat.Pos[SubPIndex[i]][2]-s[2];
			 #ifdef PERIODIC_BDR
			 dx=NEAREST(dx);
			 dy=NEAREST(dy);
			 dz=NEAREST(dz);
			 #endif
			 dx*=Time;dy*=Time;dz*=Time;
			 dvx+=Hz*dx;//add Hubble flow
			 dvy+=Hz*dy;
			 dvz+=Hz*dz;
			E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
				#pragma omp critical (unbinding)
				if(E<0)
				{
				TmpPIndex[Nbound]=SubPIndex[i];
				TmpBindingE[Nbound]=E;
				Nbound++;
				}
				else
				{
				SubPIndex_removed[SubLen_removed]=SubPIndex[i];
				SubLen_removed++;	
				}
			}
		}//\\end para
		//fprintf(fp,"%d\n",Nbound);
		free(pot);
		tree_tree_free();
		free(SubPIndex);//am i also freeing SubCatTmp->PSubArr[subhaloid]??............................Ok,fine.
		SubPIndex=TmpPIndex;//take over the bound part to begin a new loop 
		if(Nbound<NBOUNDMIN)
		{
			memcpy(SubPIndex_removed+SubLen_removed,TmpPIndex,Nbound*sizeof(int));
			*P2PIndex_removed=SubPIndex_removed;
			*P2Len_removed=SubLen_removed+Nbound;
			free(TmpPIndex);
			free(TmpBindingE);
			*P2Len=0;
			*P2PIndex=NULL;
			CoM[0]=CoM[1]=CoM[2]=0.;
			return 0;//no bound sub found, return 0;
		}			
	}
	//t1=time(NULL);
	
	SubLen=Nbound;//NOTE: this is only necessary when using tolerated (approximated) binding criteria
	SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
	for(i=0;i<SubLen;i++)
	{
		SubEsort[i].PID=SubPIndex[i];
		SubEsort[i].Erg=TmpBindingE[i];
	}
	qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
	for(i=0;i<SubLen;i++)
		SubPIndex[i]=SubEsort[i].PID;//copy the sorted index back

	free(TmpBindingE);
	free(SubEsort);
	*P2Len=SubLen;
	*P2PIndex=SubPIndex;
	CoM[0]=s[0];
	CoM[1]=s[1];
	CoM[2]=s[2]; 
	*P2Len_removed=SubLen_removed;
	*P2PIndex_removed=realloc(SubPIndex_removed,sizeof(int)*SubLen_removed);
	//time_end=time(NULL);
	//printf("%ld (%f),%ld (%f) sec\n",t1-time_start,(float)(t1-time_start)/(time_end-time_start),time_end-t1,(float)(time_end-t1)/(time_end-time_start));
	return 1;//found bound structure, return 1;
}

 /*using potential weighted CoM plus weighted VCoM to calculate kinetic energy
  * */

int unbind_potW(int *P2Len,int **P2PIndex, float *CoM,int *P2Len_removed, int **P2PIndex_removed) /*P2Len=&Len, *P2PIndex=PIndex, 
																							*where PIndex[Len] is the array of size Len; 
																							* both will be updated after unbind
																							* *P2PIndex_removed=PIndex_removed need 
																							* not to be allocated as input*/
{
	int * TmpPIndex, *SubPIndex,*SubPIndex_removed;
	float * TmpBindingE, E;
	struct Energy *SubEsort;
	int i,Nbound,SubLen,SubLen_removed,id_minpot;
	double Hz,sqa,Time,PartMass;
	double vx,vy,vz,sx,sy,sz,dvx,dvy,dvz,dx,dy,dz,*pot,sumpot;
	double mx,my,mz;
	//time_t time_start,time_end,t1,t2,t3,t4,tt1,tt2,tt3,tt4,ttt2,ttt3; //for program timing
	//~ FILE *fp;
	//fp=fopen("converg.dat","w");
	
	//time_start=time(NULL);
	Time=header.time;
	 #ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	 Hz=header.Hz;
	 PartMass=header.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

	Nbound=*P2Len;
	SubPIndex=*P2PIndex;
	SubLen_removed=0;
	SubPIndex_removed=mymalloc(sizeof(int)*Nbound);
	
	//fprintf(fp,"%d\n",Nbound);
	//if(Nbound>NParaMin) printf("%d\n",Nbound);
	SubLen=Nbound*2+1;//just to start the loop
	TmpBindingE=mymalloc(1);//just to start the loop
	while(Nbound<SubLen*PrecMass)
	{	
		//fprintf(logfile,"Nbound=%d\n",Nbound);fflush(stdout);
		SubLen=Nbound;
		Nbound=0;
		TmpPIndex=mymalloc(sizeof(int)*SubLen);
		TmpBindingE=realloc(TmpBindingE,sizeof(float)*SubLen);
		pot=mymalloc(sizeof(double)*SubLen);
		
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
		maketree(SubLen,SubPIndex,Pdat.Pos);
		vx=vy=vz=sx=sy=sz=0.;sumpot=0.;
		#pragma omp parallel  if(SubLen>NParaMin) 
		{//\\start para
			#pragma omp for schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
				pot[i]=tree_treeevaluate_potential(Pdat.Pos[SubPIndex[i]],SubPIndex,Pdat.Pos);
				pot[i]=(PartMass*pot[i]+PartMass/SofteningHalo)*G/Time;
			}
			#ifdef PERIODIC_BDR
			#pragma omp single
			{
			id_minpot=SubPIndex[Dmin_of_vec(pot,SubLen)];
			mx=Pdat.Pos[id_minpot][0];
			my=Pdat.Pos[id_minpot][1];
			mz=Pdat.Pos[id_minpot][2];
			}
			#endif
			#pragma omp for private(i,dx,dy,dz) reduction (+:vx,vy,vz,sx,sy,sz,sumpot) schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
				if(pot[i]<0)
				{
				vx-=Pdat.Vel[SubPIndex[i]][0]*pot[i];//weighted vel sum
				vy-=Pdat.Vel[SubPIndex[i]][1]*pot[i];
				vz-=Pdat.Vel[SubPIndex[i]][2]*pot[i];
				#ifdef PERIODIC_BDR
				dx=Pdat.Pos[SubPIndex[i]][0]-mx;
				dy=Pdat.Pos[SubPIndex[i]][1]-my;
				dz=Pdat.Pos[SubPIndex[i]][2]-mz;
				sx-=NEAREST(dx)*pot[i];
				sy-=NEAREST(dy)*pot[i];
				sz-=NEAREST(dz)*pot[i];
				#else
				sx-=Pdat.Pos[SubPIndex[i]][0]*pot[i];//weighted pos sum
				sy-=Pdat.Pos[SubPIndex[i]][1]*pot[i];
				sz-=Pdat.Pos[SubPIndex[i]][2]*pot[i];
				#endif
				sumpot-=pot[i];
				}
				else
				{
					printf("error: pot>0? for i=%d, pot=%g\n",i,pot[i]);
					exit(1);
				}
			}
			#pragma omp single
			{
			vx/=sumpot;
			vy/=sumpot;
			vz/=sumpot;
			sx/=sumpot;
			sy/=sumpot;
			sz/=sumpot;	
			#ifdef PERIODIC_BDR
			sx+=mx;
			sy+=my;
			sz+=mz;
			#endif
			//~ printf("%g,%g,%g\n",sumpot,vx,sx);
			OUTPUT_CENTERS
			}
			#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,E) schedule(dynamic)
			for(i=0;i<SubLen;i++)
			{
			 dvx=sqa*(Pdat.Vel[SubPIndex[i]][0]-vx);//relative vel.
			 dvy=sqa*(Pdat.Vel[SubPIndex[i]][1]-vy);
			 dvz=sqa*(Pdat.Vel[SubPIndex[i]][2]-vz);
			 dx=Pdat.Pos[SubPIndex[i]][0]-sx;
			 dy=Pdat.Pos[SubPIndex[i]][1]-sy;
			 dz=Pdat.Pos[SubPIndex[i]][2]-sz;
			 #ifdef PERIODIC_BDR
			 dx=NEAREST(dx);
			 dy=NEAREST(dy);
			 dz=NEAREST(dz);
			 #endif
			 dx*=Time;dy*=Time;dz*=Time;
			 dvx+=Hz*dx;//add Hubble flow
			 dvy+=Hz*dy;
			 dvz+=Hz*dz;
			E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
				#pragma omp critical (unbinding)
				if(E<0)
				{
				TmpPIndex[Nbound]=SubPIndex[i];
				TmpBindingE[Nbound]=E;
				Nbound++;
				}
				else
				{
				SubPIndex_removed[SubLen_removed]=SubPIndex[i];
				SubLen_removed++;	
				}
			}
		}//\\end para
		//fprintf(fp,"%d\n",Nbound);
		free(pot);
		tree_tree_free();
		free(SubPIndex);//am i also freeing SubCatTmp->PSubArr[subhaloid]??............................Ok,fine.
		SubPIndex=TmpPIndex;//take over the bound part to begin a new loop 
		if(Nbound<NBOUNDMIN)
		{
			memcpy(SubPIndex_removed+SubLen_removed,TmpPIndex,Nbound*sizeof(int));
			*P2PIndex_removed=SubPIndex_removed;
			*P2Len_removed=SubLen_removed+Nbound;
			free(TmpPIndex);
			free(TmpBindingE);
			*P2Len=0;
			*P2PIndex=NULL;
			CoM[0]=CoM[1]=CoM[2]=0.;
			return 0;//no bound sub found, return 0;
		}			
	}
	//t1=time(NULL);
	
	SubLen=Nbound;//NOTE: this is only necessary when using tolerated (approximated) binding criteria
	SubEsort=mymalloc(sizeof(struct Energy)*SubLen);
	for(i=0;i<SubLen;i++)
	{
		SubEsort[i].PID=SubPIndex[i];
		SubEsort[i].Erg=TmpBindingE[i];
	}
	qsort(SubEsort,SubLen,sizeof(struct Energy),comp_erg);
	for(i=0;i<SubLen;i++)
		SubPIndex[i]=SubEsort[i].PID;//copy the sorted index back

	free(TmpBindingE);
	free(SubEsort);
	*P2Len=SubLen;
	*P2PIndex=SubPIndex;
	CoM[0]=sx;
	CoM[1]=sy;
	CoM[2]=sz; 
	*P2Len_removed=SubLen_removed;
	*P2PIndex_removed=realloc(SubPIndex_removed,sizeof(int)*SubLen_removed);
	//time_end=time(NULL);
	//printf("%ld (%f),%ld (%f) sec\n",t1-time_start,(float)(t1-time_start)/(time_end-time_start),time_end-t1,(float)(time_end-t1)/(time_end-time_start));
	return 1;//found bound structure, return 1;
}

int unbind_add(int *P2Len,int **P2PIndex, struct SubProperty *Prop, int SubLen, int *SubArr,float SubCoM[3],float SubVCoM[3]) /*P2Len=&Len, *P2PIndex=PIndex, 
																*where PIndex[Len] is the array of size Len; 
																* both will be updated after unbind*/
{
	int * GasArr,*BndGasArr, *BndGasInd;
	float * GasE, E;
	struct Energy *SubEsort;
	int i,j,pid,gaslen,gasbound;
	double Hz,sqa,Time,DMMass,GasMass;
	double vx,vy,vz,sx,sy,sz,dvx,dvy,dvz,dx,dy,dz,*pot;
	double SubPot,SubKin,AMx,AMy,AMz,U;

	if((0==SubLen)||(0==(*P2Len))) //do nothing for null subs
	{
		*P2Len=0;
		*P2PIndex=NULL;
		Prop->CoM[0]=Prop->CoM[1]=Prop->CoM[2]=0.;
		Prop->VCoM[0]=Prop->VCoM[1]=Prop->VCoM[2]=0.;
		Prop->Pot=Prop->Kin=0.;
		Prop->AM[0]=Prop->AM[1]=Prop->AM[2]=0.;
		return 0;
	}
	
	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	Hz=header.Hz;
	GasMass=header.mass[0];
	DMMass=header.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

	gaslen=*P2Len;
	GasArr=*P2PIndex;
	gasbound=0;
	BndGasArr=mymalloc(sizeof(int)*gaslen);
	BndGasInd=mymalloc(sizeof(int)*gaslen);
	GasE=mymalloc(sizeof(float)*gaslen);
	pot=mymalloc(sizeof(double)*gaslen);

	tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
	maketree(SubLen,SubArr,Pdat.Pos);

	sx=sy=sz=vx=vy=vz=0.;
	SubPot=SubKin=AMx=AMy=AMz=0.;
	#pragma omp parallel if(SubLen>NParaMin) 
	{
	#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,E) schedule(dynamic,1) 
	for(i=0;i<gaslen;i++)
	{		
		dvx=sqa*Pdat.Vel[GasArr[i]][0]-SubVCoM[0];//relative vel.Note: SubVCoM is physical Vel
		dvy=sqa*Pdat.Vel[GasArr[i]][1]-SubVCoM[1];
		dvz=sqa*Pdat.Vel[GasArr[i]][2]-SubVCoM[2];
		dx=Pdat.Pos[GasArr[i]][0]-SubCoM[0];
		dy=Pdat.Pos[GasArr[i]][1]-SubCoM[1];
		dz=Pdat.Pos[GasArr[i]][2]-SubCoM[2];
		#ifdef PERIODIC_BDR
		dx=NEAREST(dx);
		dy=NEAREST(dy);
		dz=NEAREST(dz);
		#endif
		dx*=Time;dy*=Time;dz*=Time;
		dvx+=Hz*dx;//add Hubble flow
		dvy+=Hz*dy;
		dvz+=Hz*dz;
		
		pot[i]=tree_treeevaluate_potential(Pdat.Pos[GasArr[i]],SubArr,Pdat.Pos);
		pot[i]=DMMass*pot[i]*G/Time;//no self energy since this is gas in dm pot
		
		E=pot[i]+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);//necessary to add internal energy?

		
		#pragma omp critical (unbinding) 
		if(E<0)
		{
		BndGasInd[gasbound]=i;
		GasE[gasbound]=E;
		gasbound++;
		}
	}
	#pragma omp for private(i,pid,dvx,dvy,dvz,dx,dy,dz) \
							reduction(+:sx,sy,sz,vx,vy,vz,SubPot,SubKin,AMx,AMy,AMz) \
							schedule(dynamic,1) 
	for(i=0;i<gasbound;i++)
	{
		SubPot+=pot[BndGasInd[i]];
		
		pid=GasArr[BndGasInd[i]];
		BndGasArr[i]=pid;
		sx+=Pdat.Pos[pid][0];
		sy+=Pdat.Pos[pid][1];
		sz+=Pdat.Pos[pid][2];
		vx+=Pdat.Vel[pid][0];
		vy+=Pdat.Vel[pid][1];
		vz+=Pdat.Vel[pid][2];
				
		dvx=sqa*Pdat.Vel[pid][0]-SubVCoM[0];//relative vel.Note: VCoM is physical Vel
		dvy=sqa*Pdat.Vel[pid][1]-SubVCoM[1];
		dvz=sqa*Pdat.Vel[pid][2]-SubVCoM[2];
		dx=Pdat.Pos[pid][0]-SubCoM[0];
		dy=Pdat.Pos[pid][1]-SubCoM[1];
		dz=Pdat.Pos[pid][2]-SubCoM[2];
		#ifdef PERIODIC_BDR
		dx=NEAREST(dx);
		dy=NEAREST(dy);
		dz=NEAREST(dz);
		#endif
		dx*=Time;dy*=Time;dz*=Time;
		dvx+=Hz*dx;//add Hubble flow
		dvy+=Hz*dy;
		dvz+=Hz*dz;
		
		AMx+=dy*dvz-dz*dvy;
		AMy+=dx*dvz-dz*dvx;
		AMz+=dx*dvy-dy*dvx;
		SubKin+=0.5*(dvx*dvx+dvy*dvy+dvz*dvz);
	}
	}
	if(gasbound)
	{
	Prop->CoM[0]=sx/gasbound;Prop->CoM[1]=sy/gasbound;Prop->CoM[2]=sz/gasbound;
	Prop->VCoM[0]=sqa*vx/gasbound;Prop->VCoM[1]=sqa*vy/gasbound;Prop->VCoM[2]=sqa*vz/gasbound;
	Prop->Pot=SubPot/gasbound;
	Prop->Kin=SubKin/gasbound;
	Prop->AM[0]=AMx/gasbound;Prop->AM[1]=AMy/gasbound;Prop->AM[2]=AMz/gasbound;
	}
	else
	{
	Prop->CoM[0]=Prop->CoM[1]=Prop->CoM[2]=0.;
	Prop->VCoM[0]=Prop->VCoM[1]=Prop->VCoM[2]=0.;
	Prop->Pot=0.;
	Prop->Kin=0.;
	Prop->AM[0]=Prop->AM[1]=Prop->AM[2]=0.;
	}
	tree_tree_free();
	free(GasArr);
	free(BndGasInd);
	free(pot);
	if(gasbound)
	GasArr=realloc(BndGasArr,sizeof(int)*gasbound);//take over the bound part to begin a new loop 
	else
	GasArr=NULL;
	gaslen=gasbound;
		
	SubEsort=mymalloc(sizeof(struct Energy)*gaslen);
	for(i=0;i<gaslen;i++)
	{
		SubEsort[i].PID=GasArr[i];
		SubEsort[i].Erg=GasE[i];
	}
	qsort(SubEsort,gaslen,sizeof(struct Energy),comp_erg);
	for(i=0;i<gaslen;i++)
		GasArr[i]=SubEsort[i].PID;//copy the sorted index back

	free(GasE);
	free(SubEsort);
	*P2Len=gaslen;
	*P2PIndex=GasArr;
	
	return gaslen;
}
