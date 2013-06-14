#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "gas_stuff.h"
int hoc[NDIV][NDIV][NDIV],ll[NP_gas];
float range[3][2], step[3];

float *Rtidal;
int *GHostSub;//initialized during load_gas_data

 void load_gas_data(int Nsnap, char *SnapPath)/*we do not resolve to PID of each gas particle, 
 																			*rather we use each particle's index in the snapshot directly to refer it, and to get its properties, 
																			*and to store in gascat, so the same particle index in each gascat do not correspond to the same particle*/
{
int i,Ngas,Ndm,Nother,Nmass,dummy,dummy2;
	#ifdef SKIP
		#undef SKIP
		#undef SKIP2
		#undef CHECK
	#endif
	#define SKIP fread(&dummy,sizeof(dummy),1,fp)
	#define SKIP2 fread(&dummy2,sizeof(dummy2),1,fp)
	#define CHECK if(dummy!=dummy2){fprintf(logfile,"error!record brackets not match for Snap%d!\t%d,%d\n",Nsnap,dummy,dummy2);fflush(logfile);exit(1);} 
	FILE *fp;
	char buf[1024];

	sprintf(buf,"%s/snapshot_%03d",SnapPath,Nsnap);
	if((fp=fopen(buf,"r"))==NULL)
	{
		fprintf(logfile,"error: file open failed for %s!\n",buf);
		fflush(logfile);
		exit(1);
	}
	SKIP;
	fread(&headerA,sizeof(headerA),1,fp);
	SKIP2;
	CHECK;

	headerA_Hz=Hubble * sqrt(headerA.Omega0 / (headerA.time * headerA.time * headerA.time) 
			+ (1 - headerA.Omega0 - headerA.OmegaLambda) / (headerA.time * headerA.time)
			+ headerA.OmegaLambda);//Hubble param for the current catalogue;
	
	Ngas=headerA.npart[0];	Ndm=headerA.npart[1];	Nother=headerA.npart[2]+headerA.npart[3]+headerA.npart[4]+headerA.npart[5];
	for(Nmass=0,i=0;i<6;i++)
	{
		if(!headerA.mass[i]) Nmass+=headerA.npart[i];
	}

	SKIP;
	fread(Gdat.Pos,sizeof(float)*3,Ngas,fp);
	fseek(fp,(Ndm+Nother)*sizeof(float)*3,SEEK_CUR);
	SKIP2;
	CHECK;
	
	SKIP;
	fread(Gdat.Vel,sizeof(float)*3,Ngas,fp);
	fseek(fp,(Ndm+Nother)*sizeof(float)*3,SEEK_CUR);
	SKIP2;
	CHECK;
	
	SKIP;
	fread(Gdat.PID,sizeof(int),Ngas,fp);
	fseek(fp,(Ndm+Nother)*sizeof(int),SEEK_CUR);
	SKIP2;
	CHECK;
	
	if(Nmass)
	{
		SKIP;
		fseek(fp,Nmass*sizeof(float),SEEK_CUR);
		SKIP2;
		CHECK;
	}
	
	SKIP;
	fread(Gdat.U,sizeof(float),Ngas,fp);
	SKIP2;
	CHECK;

	SKIP;
	fread(Gdat.Rho,sizeof(float),Ngas,fp);
	SKIP2;
	CHECK;
	
	fclose(fp);
	#undef SKIP
	#undef SKIP2
	#undef CHECK
	GHostSub=mymalloc(sizeof(int)*Ngas);
	for(i=0;i<Ngas;i++)
		GHostSub[i]=-1;
}
	
double tree_treeevaluate_potential_gas(int target, int *PIndex)
{
	//target is PID rather than Pindex
  union NODE *nop = 0;
  int no;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pot, pos_x, pos_y, pos_z;

  pos_x = Gdat.Pos[target][0];
  pos_y = Gdat.Pos[target][1];
  pos_z = Gdat.Pos[target][2];//this is a gas particle

  h = 2.8 * SofteningHalo;
  h_inv = 1.0 / h;

  pot = 0;

  no = NumPart;//start from root node

  while(no >= 0)
    {
      if(no < NumPart)		/* single particle */
	{
	  dx = Pdat.Pos[PIndex[no]][0] - pos_x;
	  dy = Pdat.Pos[PIndex[no]][1] - pos_y;
	  dz = Pdat.Pos[PIndex[no]][2] - pos_z;
	  mass = 1;
	  no = Nextnode[no];
	      r2 = dx * dx + dy * dy + dz * dz;	
	}
      else
	{
	  nop = &Nodes[no];
	  dx = nop->way.s[0] - pos_x;
	  dy = nop->way.s[1] - pos_y;
	  dz = nop->way.s[2] - pos_z;
	  mass = nop->way.mass;
	  r2 = dx * dx + dy * dy + dz * dz;
		/* we have an internal node. Need to check opening criterion */
	  if(nop->way.len * nop->way.len > r2 * ErrTolTheta * ErrTolTheta)
	    {
	      /* open cell */
	      no = nop->way.nextnode;
	      continue;
	    }
	  no = nop->way.sibling;	/* node can be used */
	}
 
      r = sqrt(r2);

      if(r >= h)
	pot -= mass / r;
      else
	{
	  u = r * h_inv;

	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp =
	      -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						   u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

	  pot += mass * h_inv * wp;
	}
    }

  return pot;
}

int collect_gas_particles(int subid,SUBCATALOGUE*SubCat,GASCATALOGUE *GasCat)
{
	int gaslen,i,j,k,pid,subbox_grid[3][2],maxlen;
	float *cen,rscale,dr;
	//~ extern float range[3][2],step[3],hoc[][][],ll[],Rtidal[];
	
	gaslen=0;
	if(SubCat->SubLen[subid])
	{
	maxlen=SubCat->SubLen[subid]*2;
	maxlen=((maxlen<NP_gas)?maxlen:NP_gas);
	GasCat->PSubArr[subid]=mymalloc(sizeof(int)*maxlen);
	cen=Pdat.Pos[SubCat->PSubArr[subid][0]];
	rscale=Rtidal[subid]*GSub_ScaleRelax;
	for(i=0;i<3;i++)
	{
	subbox_grid[i][0]=floor((cen[i]-rscale-range[i][0])/step[i]);
	subbox_grid[i][1]=floor((cen[i]+rscale-range[i][0])/step[i]);
	}
	//~ printf("%d,%d,%d,%d,%d,%d\n",subbox_grid[0][0],subbox_grid[0][1],subbox_grid[1][0],subbox_grid[1][1],subbox_grid[2][0],subbox_grid[2][1]);fflush(stdout);
	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				pid=hoc[i][j][k];
				while(pid>=0)
				{
					if(GHostSub[pid]<0)//have not been captured by other subs 
					{
						dr=distance(Gdat.Pos[pid],cen);
						if(dr<rscale)
						{
							GasCat->PSubArr[subid][gaslen]=pid;
							gaslen++;
							if(gaslen>=maxlen)
							{
								maxlen*=2;
								GasCat->PSubArr[subid]=realloc(GasCat->PSubArr[subid],sizeof(int)*maxlen);
							}
						}
					}
					pid=ll[pid];
				}
			}
			GasCat->SubLen[subid]=gaslen;
			GasCat->PSubArr[subid]=realloc(GasCat->PSubArr[subid],sizeof(int)*gaslen);
	}
	else
	{
		GasCat->SubLen[subid]=gaslen;
		GasCat->PSubArr[subid]=NULL;
	}
			return gaslen;
}
static int comp_erg(const void *a, const void *b)//used to sort energy in ascending order; note that the most bound will come first( energy <0)
{
  if(((struct Energy *) a)->Erg > ((struct Energy *) b)->Erg)
    return +1;

  if(((struct Energy *) a)->Erg < ((struct Energy *) b)->Erg)
    return -1;

  return 0;
}
int unbindgas(int *P2Len,int **P2PIndex, int SubLen, int *SubArr) /*P2Len=&Len, *P2PIndex=PIndex, 
																*where PIndex[Len] is the array of size Len; 
																* both will be updated after unbind*/
{
	int * GasArr,*TmpGasArr;
	float * GasE, E;
	struct Energy *SubEsort;
	int i,j,gaslen,gasbound;
	double Hz,sqa,Time,DMMass,GasMass;
	double vx,vy,vz,s[3],dvx,dvy,dvz,dx,dy,dz,pot;

if((SubLen==0)||((*P2Len)==0)) //do nothing for null subs
	return 0;
	
	Time=headerA.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#end
	 Hz=headerA_Hz;
	 GasMass=headerA.mass[0];
	 DMMass=headerA.mass[1];
	//fprintf(logfile,"unbinding.............\n");fflush(stdout);

		gaslen=*P2Len;
		GasArr=*P2PIndex;
		gasbound=0;
		TmpGasArr=mymalloc(sizeof(int)*gaslen);
		GasE=mymalloc(sizeof(float)*gaslen);
			
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubLen,SubLen);
		maketree(SubLen,SubArr);////////////////////////to be modified here.
		vx=vy=vz=0.;

	#ifdef OMP
		#pragma omp parallel   if(SubLen>NParaMin) 
	#endif
		{//\\start para
		#ifdef OMP
			#pragma omp for nowait 
		#endif
			for(j=0;j<3;j++)
			{
			s[j]=Nodes_base->way.s[j];//Center of Mass
			}
		#ifdef OMP
			#pragma omp for reduction (+:vx,vy,vz) schedule(dynamic)
		#endif		
			for(i=0;i<SubLen;i++)
			{
				vx+=Pdat.Vel[SubArr[i]][0];
				vy+=Pdat.Vel[SubArr[i]][1];
				vz+=Pdat.Vel[SubArr[i]][2];
			}
		#ifdef OMP
			#pragma omp single
		#endif
			{
			vx/=SubLen;
			vy/=SubLen;
			vz/=SubLen;
			}
		#ifdef OMP
			#pragma omp for private(i,dvx,dvy,dvz,dx,dy,dz,pot,E) schedule(dynamic)
		#endif
			for(i=0;i<gaslen;i++)
			{
			 dvx=sqa*(Gdat.Vel[GasArr[i]][0]-vx);//relative vel.
			 dvy=sqa*(Gdat.Vel[GasArr[i]][1]-vy);
			 dvz=sqa*(Gdat.Vel[GasArr[i]][2]-vz);
			 dx=Time*(Gdat.Pos[GasArr[i]][0]-s[0]);
			 dy=Time*(Gdat.Pos[GasArr[i]][1]-s[1]);
			 dz=Time*(Gdat.Pos[GasArr[i]][2]-s[2]);
			 dvx+=Hz*dx;//add Hubble flow
			 dvy+=Hz*dy;
			 dvz+=Hz*dz;
			pot=tree_treeevaluate_potential_gas(GasArr[i],SubArr);////////////////////////to be modified here;
			pot=(DMMass*pot)*G/Time;//gas in DM pot, no self-energy included so no need to correct for the self-energy
			#ifdef THERMAL_BOUND	
			E=pot+0.5*(dvx*dvx+dvy*dvy+dvz*dvz)+Gdat.U[GasArr[i]];//necessary to add internal energy?
			#else
			E=pot+0.5*(dvx*dvx+dvy*dvy+dvz*dvz);//necessary to add internal energy?
			#endif
			
			#ifdef OMP
				#pragma omp critical (unbinding)
			#endif
				if(E<0)
				{
				TmpGasArr[gasbound]=GasArr[i];
				GasE[gasbound]=E;
				gasbound++;
				}
			}
		}//\\end para
		//~ #undef OMP
		
		tree_tree_free();
		free(GasArr);
		GasArr=realloc(TmpGasArr,sizeof(int)*gasbound);//take over the bound part to begin a new loop 
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

void unbind_gas_recursive(int mainsubID,SUBCATALOGUE *SubCat,GASCATALOGUE *GasCat)
{
	int son,sib,i;
	
	if((son=SubCat->sub_hierarchy[mainsubID].sub)>=0)//has sons,unbind them recursively to update its reservior
	{	
		sib=son;
		while(sib>=0)
		{
			unbind_gas_recursive(sib,SubCat,GasCat);
			sib=SubCat->sub_hierarchy[sib].next;
		}
	}
	//now unbind the mainsub itself
	collect_gas_particles(mainsubID,SubCat,GasCat);
	unbindgas(GasCat->SubLen+mainsubID,GasCat->PSubArr+mainsubID,SubCat->SubLen[mainsubID],SubCat->PSubArr[mainsubID]);
	#ifdef EXCLUSIVE_GAS
	for(i=0;i<GasCat->SubLen[mainsubID];i++)		//register particle hosts
		GHostSub[GasCat->PSubArr[mainsubID][i]]=mainsubID;
	#endif
		
}
void create_gas_cat(int Nsubs,GASCATALOGUE *GasCat)
{
	GasCat->Nsubs=Nsubs;
	GasCat->SubLen=mymalloc(sizeof(int)*Nsubs);
	GasCat->SubOffset=mymalloc(sizeof(int)*Nsubs);
	GasCat->PSubArr=mymalloc(sizeof(int *)*Nsubs);
}
void save_gas_cat(int Nsnap,GASCATALOGUE *Cat,char *gasdir)
{	
FILE *fd;
char buf[1024];
int i;

  sprintf(buf, "%s/gascat_%03d", gasdir, Nsnap);
  if(!(fd = fopen(buf, "w")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fwrite(&Cat->Nsubs, sizeof(int), 1, fd);
  fwrite(&Cat->Nids, sizeof(int), 1, fd);
  fwrite(Cat->SubLen, sizeof(int), Cat->Nsubs, fd);
  fwrite(Cat->SubOffset,sizeof(int), Cat->Nsubs, fd);
for(i=0;i<Cat->Nsubs;i++)
  fwrite(Cat->PSubArr[i], sizeof(int), Cat->SubLen[i], fd);

  fclose(fd);
}
void load_gas_cat(int Nsnap,GASCATALOGUE *Cat,char *subdir)
{	
FILE *fd;
char buf[1024];
int i;

  sprintf(buf, "%s/gascat_%03d", subdir, Nsnap);
  if(!(fd = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fread(&Cat->Nsubs, sizeof(int), 1, fd);
  fread(&Cat->Nids, sizeof(int), 1, fd);
  create_gas_cat(Cat->Nsubs,Cat);
  fread(Cat->SubLen, sizeof(int), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(int), Cat->Nsubs, fd);
for(i=0;i<Cat->Nsubs;i++)
{
	Cat->PSubArr[i]=mymalloc(sizeof(int)*Cat->SubLen[i]);
  fread(Cat->PSubArr[i], sizeof(int), Cat->SubLen[i], fd);
}

  fclose(fd);
}

void free_gas_cat(GASCATALOGUE *Cat)
{
	int i;
	free(Cat->SubLen);
	free(Cat->SubOffset);
	for(i=0;i<Cat->Nsubs;i++)
		free(Cat->PSubArr[i]);
	free(Cat->PSubArr);
}

void load_tidal_radius(int Nsnap,float *rtidal, int Nsubs,char*tidaldir)
{
	char buf[1024];
FILE *fp;
	int i;
	
	sprintf(buf,"%s/tidal_radius_%03d",tidaldir,Nsnap);
	if(!(fp=fopen(buf,"r")))
	{
		fprintf(logfile,"Error opening file '%s'\n",buf);
		exit(1);
	}
	for(i=0;i<Nsubs;i++)	
	{
		fscanf(fp,"%e\n",rtidal+i);
		if(ferror(fp))
		{
		fprintf(logfile,"error reading file %s: i=%d,Nsubs=%d\n",buf,i,Nsubs);
		exit(1);
		}
	}
	fclose(fp);
}

void makell(float (*pos)[3],int np)
{
	/*output: hoc[NDIV][NDIV][NDIV],ll[np]*/
	int i,j,grid[3];
	//~ float range[3][2],step[3];
	printf("creating linked list..\n");
	
	/*determining enclosing cube*/
	for(i=0;i<3;i++)
		for(j=0;j<2;j++)
			range[i][j]=pos[0][i];
	for(i=1;i<np;i++)
		for(j=0;j<3;j++)
		{
			if(pos[i][j]<range[j][0])
				range[j][0]=pos[i][j];
			else if(pos[i][j]>range[j][1])
				range[j][1]=pos[i][j];
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
			grid[j]=floor((pos[i][j]-range[j][0])/step[j]);
			if(grid[j]<0) 
				grid[j]=0;
			else if(grid[j]>=NDIV)
				grid[j]=NDIV-1;
		}
		ll[i]=hoc[grid[0]][grid[1]][grid[2]];
		hoc[grid[0]][grid[1]][grid[2]]=i; /*use hoc(floor(xsat)+1) as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
}

