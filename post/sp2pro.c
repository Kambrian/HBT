#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
	
int main()
{
	int i, tmp[4],Nsnap=0;
	int *sp2pro,Npro;
	SUBCATALOGUE SubCatB;
	
	FILE *fp;
	FILE *fpsp; 
	char buf[1024];
	char bufsp[1024]; 
	
	char inputdir[512]=SUBCAT_DIR;
	char outputdir[512];
	logfile=stdout;
	sprintf(outputdir,"%s/splitters",SUBCAT_DIR);
	sprintf(bufsp, "%s/splitters.log",inputdir);
	myfopen(fpsp,bufsp,"r");
	
	Npro=0;
	fscanf(fpsp,"%d,%d,%d,%d\n",tmp,tmp+1,tmp+2,tmp+3);
	
	for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
	{
		load_sub_catalogue(Nsnap,&SubCatB,inputdir);	
		sprintf(buf, "%s/sp2pro_%03d",outputdir,Nsnap);
		myfopen(fp,buf,"w");
		fwrite(&SubCatB.Nsplitter,sizeof(int),1,fp);
		fwrite(&Npro,sizeof(int),1,fp);
		//prepare splitter2sub table
		if(SubCatB.Nsplitter)//subcatLAST splits
		{
			if(tmp[1]!=(Nsnap-1))
			{
				printf("error reading splitter info, snap=%d,%d\n",Nsnap-1,tmp[1]);
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
			fwrite(sp2pro,sizeof(int),SubCatB.Nsplitter,fp);
			free(sp2pro);
		}
		fclose(fp);
		Npro=SubCatB.Nsubs;
		free_sub_catalogue(&SubCatB);
	}
	fclose(fpsp);
return 0;
}
