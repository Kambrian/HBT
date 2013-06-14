//to forge a srccat for the major merger simulations
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"


//~ #define INPUT_FOF
#ifdef INPUT_FOF
#define NGRPS 2148 
#else
#define NGRPS 2
#endif

int main()
{
	char buf[1024];
	FILE *fp;
	SRCCATALOGUE SrcCat;
	SUBCATALOGUE SubCat;
	HBTInt *PIDs,grplen,grpid,nid;
	float tmp;
	
	logfile=stdout;
	
	SrcCat.Nsubs=NGRPS;
	create_src_cat(&SrcCat);
	
	SubCat.Ngroups=NGRPS;
	SubCat.GrpOffset_Sub=mymalloc(sizeof(HBTInt)*NGRPS);
	SubCat.GrpLen_Sub=mymalloc(sizeof(HBTInt)*NGRPS);
	SubCat.Nsubs=NGRPS;	create_sub_cat(&SubCat);
	
	nid=0;
	for(grpid=0;grpid<NGRPS;grpid++)
	{
	SubCat.GrpLen_Sub[grpid]=1;
	SubCat.GrpOffset_Sub[grpid]=grpid;
	SubCat.SubLen[grpid]=0;
	SubCat.SubOffset[grpid]=0;
	SubCat.SubRank[grpid]=0;
	SubCat.HaloChains[grpid].HostID=grpid;	
	SubCat.HaloChains[grpid].ProSubID=-1;	
	SubCat.sub_hierarchy[grpid].nibs=-1;
	SubCat.sub_hierarchy[grpid].pre=-1;
	SubCat.sub_hierarchy[grpid].next=-1;
	SubCat.sub_hierarchy[grpid].sub=-1;
	SubCat.NQuasi=0;
	SubCat.Nsplitter=0;
	
	PIDs=mymalloc(sizeof(HBTInt)*NP_DM);
	#ifdef INPUT_FOF
	sprintf(buf,"%s/fofs/fof.%d",SNAPSHOT_DIR,grpid);
	#else	
	sprintf(buf,"%s/fofs/grp%d",SNAPSHOT_DIR,grpid);
	#endif
	myfopen(fp,buf,"r");
	grplen=0;
	SrcCat.SubOffset[grpid]=nid;	
	while(1)
	{
		if(fscanf(fp,HBTIFMT" %f %f %f\n", PIDs+grplen,&tmp,&tmp,&tmp)!=EOF)
		grplen++;
		else
		break;
	}	
	SrcCat.SubLen[grpid]=grplen;
	SrcCat.SubLen2[grpid]=0;
	SrcCat.CoreFrac[grpid]=CoreFrac0;
	SrcCat.PSubArr[grpid]=realloc(PIDs,sizeof(HBTInt)*grplen);
	nid+=grplen;
	}
	SrcCat.Nids=nid;
	SrcCat.NDeathSp=0;
	
	save_src_catalogue(0,&SrcCat,SUBCAT_DIR);
	save_sub_catalogue(0,&SubCat,SUBCAT_DIR);
	
			
	return 0;
}
