#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "globals.c"
#include "mymath.c"
#include "load_group.c"
#include "tree.c"
#include "binding_minpot.c"

int main(int nargc,char **argv)
{
char subdir[512]="/SANdisk5/kambrain/Sim6702/SubCat5";
char fofdir[512]="/SANdisk5/kambrain/Sim6702/FoFCat"; //"/home/kambrain/fof_hy";
char snapdir[512]="/SANdisk4/data/NewR/SIM6702";
	
int SnapshotNum=0,SnapBind;
CATALOGUE CatB;
SUBCATALOGUE SubCat;
MainIniCat main_sub_ini_cat;//THREE phases of main_sub halos: IN:main_sub_ini_cat; AT: SubCat; OUT: MaskedFoFArr;
int **MaskedFoFArr, *MaskedFoFLen;//to store the masked FoF
int GrpLen, *GrpPIDs,Lmain_removed,*PInd_main_removed;
int subhaloid,fofid,subid,i;
float s[3];
logfile=stdout;	

fofid=atoi(argv[1]);
SnapshotNum=atoi(argv[2]);
SnapBind=atoi(argv[3]);
CoreFrac=atof(argv[4]);

	//===pick up a previous subcat to continue,assign snapnum first==//	

		load_sub_catalogue(SnapshotNum,&SubCat,subdir);
		load_masked_fof(SnapshotNum,&MaskedFoFLen,&MaskedFoFArr,subdir);
		load_main_sub_ini(SnapshotNum,&main_sub_ini_cat,subdir);
		load_group_catalogue(SnapshotNum,&CatB,fofdir);
		load_particle_data(SnapBind,&CatB,&SubCat,snapdir);
		subhaloid=SubCat.GrpOffset_Sub[fofid];
		printf("foflen %d,inilen %d, sublen %d, sublen2 %d, maskedlen %d, nsub %d\n",CatB.Len[fofid],main_sub_ini_cat.Len[fofid],SubCat.SubLen[subhaloid],SubCat.SubLen[subhaloid+1],MaskedFoFLen[fofid],SubCat.GrpLen_Sub[fofid]);
		
		GrpLen=CatB.Len[fofid];
		GrpPIDs=mymalloc(sizeof(int)*GrpLen);
		memcpy(GrpPIDs,CatB.PIDorIndex+CatB.Offset[fofid],sizeof(int)*GrpLen);
		unbind(&GrpLen,&GrpPIDs,s,&Lmain_removed,&PInd_main_removed);	
		printf("after unbinding	len=%d\n",GrpLen);
		free(PInd_main_removed);
		free_catalogue(&CatB);
		free_sub_catalogue(&SubCat);
		
		load_sub_catalogue(SnapshotNum,&SubCat,subdir);
		load_group_catalogue(SnapshotNum,&CatB,fofdir);
		for(i=0;i<SubCat.Ngroups;i++)
		{
		subid=SubCat.GrpOffset_Sub[i];
		free(SubCat.PSubArr[subid]);
		SubCat.PSubArr[subid]=MaskedFoFArr[i];
		SubCat.SubLen[subid]=MaskedFoFLen[i];
		}
		load_particle_data(SnapBind,&CatB,&SubCat,snapdir);
		GrpLen=SubCat.SubLen[subhaloid];
		GrpPIDs=mymalloc(sizeof(int)*GrpLen);
		memcpy(GrpPIDs,SubCat.PSubArr[subhaloid],sizeof(int)*GrpLen);
		unbind(&GrpLen,&GrpPIDs,s,&Lmain_removed,&PInd_main_removed);	
		printf("after unbinding with mask	len=%d\n",GrpLen);
	
		return 0;
}
