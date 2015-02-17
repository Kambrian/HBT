//to forge a grpcat for the major merger simulations
#include <stdio.h>
#include <stdlib.h>


int main()
{
	char buf[1024];
	FILE *fp;
	SRCCATALOGUE SrcCat;
	HBTInt *PIDs,i,len0,len1;
	float tmp;
	
	SrcCat.Nsubs=NGRPS;
	create_src_cat(&SrcCat);
	load_particle_data(0,SNAPSHOT_DIR);
	SrcCat.PSubArr[0]=mymalloc(sizeof(HBTInt)*NP_DM);
	SrcCat.PSubArr[1]=mymalloc(sizeof(HBTInt)*NP_DM);
	for(i=0,len0=0,len1=0;i<NP_DM;i++)
	{
		if(Pdat.Pos[i][0]<55.)
		{
		SrcCat.PsubArr[0][len0]=Pdat.PID[i];
		len0++;
		}
		else
		{
		SrcCat.PsubArr[1][len1]=Pdat.PID[i];
		len1++;	
		}
	}
	SrcCat.SubLen[0]=len0;
	SrcCat.SubLen[1]=len1;
	SrcCat.SubLen2[0]=0;
	SrcCat.SubLen2[1]=0;
	SrcCat.Offset[0]=0;
	SrcCat.Offset[1]=len0;
	SrcCat.CoreFrac[0]=CoreFrac0;
	SrcCat.CoreFrac[1]=CoreFrac0;
	SrcCat.Nids=NP_DM;
	SrcCat.NDeathSp=0;
	printf("%d,%d:%d\n",len0,len1,NP_DM);
	
	save_src_catalogue(0,&SrcCat,SUBCAT_DIR);
	
	return 0;
}
