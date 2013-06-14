//#define NDIV 200
HBTInt hoc_DM[NDIV][NDIV][NDIV],ll_DM[NP_DM];
HBTReal range_DM[3][2],step_DM[3];
void makell_DM()
{
	HBTReal (*pos)[3];
	HBTInt np;
	HBTInt i,j,grid[3];
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
	HBTInt *phoc=&(hoc_DM[0][0][0]);
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

