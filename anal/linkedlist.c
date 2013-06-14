//#define NDIV 200
HBTInt hoc[NDIV][NDIV][NDIV],ll[NP_DM];
HBTReal range[3][2],step[3];

#ifndef PERIODIC_BDR
void makell(HBTReal (*pos)[3],HBTInt np,int ndiv)
{
	/*output: hoc[ndiv][ndiv][ndiv],ll[np]*/
	HBTInt i,j,grid[3];
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
		step[j]=(range[j][1]-range[j][0])/ndiv;
	
	/*initialize hoc*/
	HBTInt *phoc=&(hoc[0][0][0]);
	for(i=0;i<ndiv*ndiv*ndiv;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((pos[i][j]-range[j][0])/step[j]);
			if(grid[j]<0) 
				grid[j]=0;
			else if(grid[j]>=ndiv)
				grid[j]=ndiv-1;
		}
		ll[i]=hoc[grid[0]][grid[1]][grid[2]];
		hoc[grid[0]][grid[1]][grid[2]]=i; /*use hoc(floor(xsat)+1) as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
}
#define GRID_IN_BOX(x) ((x)<0?0:((x)>=NDIV?(NDIV-1):(x)))
#else
void makell(HBTReal (*pos)[3],HBTInt np,int ndiv)
{//apply to entire box, and compatible with periodic boundaries

	/*output: hoc[ndiv][ndiv][ndiv],ll[np]*/
	HBTInt i,j,grid[3];
	//~ float range[3][2],step[3];
	printf("creating linked list..\n");
	
	/*determining enclosing cube*/
	for(i=0;i<3;i++)
	{
		range[i][0]=0.;
		range[i][1]=BOXSIZE;
	}
	for(j=0;j<3;j++)
		step[j]=BOXSIZE/ndiv;
	
	/*initialize hoc*/
	HBTInt *phoc=&(hoc[0][0][0]);
	for(i=0;i<ndiv*ndiv*ndiv;i++,phoc++)
		*phoc=-1;
		
	for(i=0;i<np;i++)
	{
		for(j=0;j<3;j++)
		{
			grid[j]=floor((pos[i][j]-range[j][0])/step[j]);
			if(grid[j]<0) 
				grid[j]=0;
			else if(grid[j]>=ndiv)
				grid[j]=ndiv-1;
		}
		ll[i]=hoc[grid[0]][grid[1]][grid[2]];
		hoc[grid[0]][grid[1]][grid[2]]=i; /*use hoc(floor(xsat)+1) as swap varible to temporarily 
			      						store last ll index, and finally the head*/
	}
}
#define GRID_IN_BOX(x) ((x)<0?((x)+NDIV):((x)>=NDIV?((x)-NDIV):(x)))
#endif


