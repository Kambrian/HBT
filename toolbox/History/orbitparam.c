#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "gas_vars.h"
#include "proto.h"
#include "gas_proto.h"
#include "history_vars.h"
#include "history_proto.h"

int main(int argc,char **argv)
{
SUBCATALOGUE SubCat;
char buf[1024];FILE *fp;
int Nsnap,subid,hsubid,grpid,i;
float r,v,kin,pot,J2;
float pos[3],vel[3],AM[3];
float Time,sqa,Hz,PartMass;
struct orbparam
{
	float r;//physical
	float v;//physical
	float kin;
	float pot;
	float J2;
} *param;
logfile=stdout;

for(Nsnap=0;Nsnap<MaxSnap;Nsnap++)
{
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	load_particle_data(Nsnap,SNAPSHOT_DIR);
	fresh_ID2Index(&SubCat,FRSH_SUBCAT);
	Time=header.time;
	#ifdef VEL_INPUT_PHYSICAL
	sqa=1.0;
	#else
	sqa = sqrt(Time);
	#endif
	 Hz=header.Hz;
	 PartMass=header.mass[1];
	 param=mymalloc(sizeof(struct orbparam)*SubCat.Nsubs);
	for(grpid=0;grpid<SubCat.Ngroups;grpid++)
	{
		hsubid=SubCat.GrpOffset_Sub[grpid];
		if(SubCat.SubLen[hsubid])
		{
		tree_tree_allocate(TREE_ALLOC_FACTOR*SubCat.SubLen[hsubid],SubCat.SubLen[hsubid]);
		maketree(SubCat.SubLen[hsubid],SubCat.PSubArr[hsubid],Pdat.Pos);
		for(subid=hsubid;subid<hsubid+SubCat.GrpLen_Sub[grpid];subid++)
		{
			pot=tree_treeevaluate_potential(SubCat.Property[subid].CoM,SubCat.PSubArr[hsubid],Pdat.Pos);
			pot=PartMass*pot*G/Time;
			for(i=0;i<3;i++)
			{
				pos[i]=SubCat.Property[subid].CoM[i]-SubCat.Property[hsubid].CoM[i];
				pos[i]*=Time;
				vel[i]=SubCat.Property[subid].VCoM[i]-SubCat.Property[hsubid].VCoM[i];//already physical
			}
			r=sqrt(f_prod(pos,pos,3));
			v=sqrt(f_prod(vel,vel,3));
			kin=0.5*f_prod(vel,vel,3);
			f_cross(pos,vel,AM);
			J2=f_prod(AM,AM,3);
			param[subid].r=r;
			param[subid].v=v;
			param[subid].kin=kin;
			param[subid].pot=pot;
			param[subid].J2=J2;
		}
		tree_tree_free();
		}
		else
		{
		param[hsubid].r=0;
		param[hsubid].v=0;
		param[hsubid].kin=0;
		param[hsubid].pot=0;
		param[hsubid].J2=0;
		}
	}
	for(subid=SubCat.Nsubs-SubCat.NQuasi;subid<SubCat.Nsubs;subid++)
	{
		param[subid].r=0;
		param[subid].v=0;
		param[subid].kin=0;
		param[subid].pot=0;
		param[subid].J2=0;
	}
	sprintf(buf,"%s/anal/orbpar/orbparam_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"w");
	fwrite(param,sizeof(struct orbparam),SubCat.Nsubs,fp);
	fclose(fp);
	free(param);
	erase_sub_catalogue(&SubCat);
}
return 0;
}
