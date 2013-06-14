#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "intra_vars.h"
#include "iovars.h"
#include "gas_vars.h"
#include "proto.h"
#include "gas_proto.h"

int main(int nargc,char **argv)
{
char gasdir[1024];
char outputdir[1024];

int SnapLoad,SnapBind;
GASHALOCAT GCat;
SUBCATALOGUE SubCat;
struct GasProperty Prop;
int *HaloDest,*SubDest,*GasLoc,HaloDestLen,SubDestLen,GasLocLen,BndLocLen,Lmain_removed,*PInd_main_removed;
int fofid,subid,i,hostid;
int *Pro2dest,Npro;
float dr,rvir;
float s[3];
  FILE *fp,*fpsub,*fpdes,*fpsrc,*fpdm,*fpcen,*fploc,*fpbnd,*fpsubcen;
  char buf[1024];

logfile=stdout;	
if(nargc!=3)
{
	printf("usage: %s [SnapLoad] [fofid]\n",argv[0]);
	exit(1);
}
SnapLoad=atoi(argv[1]);
fofid=atoi(argv[2]);

  sprintf(gasdir,"%s/gascat",SUBCAT_DIR);
  sprintf(outputdir,"%s/anal/trace/follow",SUBCAT_DIR);
  #ifdef THERMAL_BOUND
  sprintf(buf, "%s/follow_gas_%03d_%d.thermal",outputdir,SnapLoad,fofid);
  #else
    sprintf(buf, "%s/follow_gas_%03d_%d",outputdir,SnapLoad,fofid);
#endif
  if(!(fp= fopen(buf, "w")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }
	
	sprintf(outputdir,"%s/S%03dG%d",outputdir,SnapLoad,fofid);
	mkdir(outputdir,0700);
	sprintf(buf, "%s/hostcen",outputdir);
	myfopen(fpcen,buf,"w");
	sprintf(buf, "%s/subcen",outputdir);
	myfopen(fpsubcen,buf,"w");
	fprintf(fp,"SnapBind,HaloDestLen,SubDestLen,BndLocLen,GasLocLen,DMSubLen,HostSubLen,hostid,subid,subrank,rcen,rcen/rvir\n");

		load_gashalocat(SnapLoad,&GCat,GASCAT_DIR);
		load_sub_catalogue(SnapLoad,&SubCat,SUBCAT_DIR);
		subid=SubCat.GrpOffset_Sub[fofid];
		for(i=0;i<SubCat.Nsubs;i++)
		{
			myfree(SubCat.PSubArr[i]);
		}
		free_sub_catalogue(&SubCat);
		/*==init HaloDest==*/
		HaloDestLen=GCat.Len[fofid];
		HaloDest=mymalloc(sizeof(int)*HaloDestLen);
		memcpy(HaloDest,GCat.PIDorIndex+GCat.Offset[fofid],sizeof(int)*HaloDestLen);
		/*==init SubDest==*/
		SubDestLen=HaloDestLen;
		SubDest=mymalloc(sizeof(int)*SubDestLen);
		memcpy(SubDest,HaloDest,sizeof(int)*HaloDestLen);
		for(SnapBind=SnapLoad;SnapBind<100;SnapBind++)
		{
			sprintf(buf, "%s/gasloc.%03d",outputdir,SnapBind);
			myfopen(fploc,buf,"w");
			sprintf(buf, "%s/gasbnd.%03d",outputdir,SnapBind);
			myfopen(fpbnd,buf,"w");
			sprintf(buf, "%s/gassrc.%03d",outputdir,SnapBind);
			myfopen(fpsrc,buf,"w");
			sprintf(buf, "%s/gasdes.%03d",outputdir,SnapBind);
			myfopen(fpdes,buf,"w");
			sprintf(buf, "%s/gassub.%03d",outputdir,SnapBind);
			myfopen(fpsub,buf,"w");
			sprintf(buf, "%s/dmpos.%03d",outputdir,SnapBind);
			myfopen(fpdm,buf,"w");
			printf("Snap %d\n",SnapBind);fflush(stdout);
			load_particle_data(SnapBind,SNAPSHOT_DIR);
			load_sub_catalogue(SnapBind,&SubCat,SUBCAT_DIR);
			fresh_ID2Index(SubCat.PSubArr[subid],SubCat.SubLen[subid]);
			for(i=0;i<SubCat.SubLen[subid];i++)
			{
				fprintf(fpdm,"%g\t%g\t%g\n",Pdat.Pos[SubCat.PSubArr[subid][i]][0],Pdat.Pos[SubCat.PSubArr[subid][i]][1],Pdat.Pos[SubCat.PSubArr[subid][i]][2]);
			}
			
			load_gas_data(SnapBind,SNAPSHOT_DIR);
			/*==local gas cat==*/
			makell_gas();
			collect_gas_particles(subid,&SubCat,&GasLocLen,&GasLoc);
			for(i=0;i<GasLocLen;i++)
			{
				fprintf(fploc,"%g\t%g\t%g\n",Gdat.Pos[GasLoc[i]][0],Gdat.Pos[GasLoc[i]][1],Gdat.Pos[GasLoc[i]][2]);
			}
			BndLocLen=GasLocLen;
			unbindgas(&BndLocLen,&GasLoc,&Prop,SubCat.SubLen[subid],SubCat.PSubArr[subid],SubCat.Property[subid].CoM,SubCat.Property[subid].VCoM);	
			for(i=0;i<BndLocLen;i++)
			{
				fprintf(fpbnd,"%g\t%g\t%g\n",Gdat.Pos[GasLoc[i]][0],Gdat.Pos[GasLoc[i]][1],Gdat.Pos[GasLoc[i]][2]);
			}
			myfree(GasLoc);
			/*==trace gas cat==*/
			fresh_gasID2Index(HaloDest,HaloDestLen);
			fresh_gasID2Index(SubDest,SubDestLen);
			for(i=0;i<HaloDestLen;i++)
			{
				fprintf(fpsrc,"%g\t%g\t%g\n",Gdat.Pos[HaloDest[i]][0],Gdat.Pos[HaloDest[i]][1],Gdat.Pos[HaloDest[i]][2]);
			}
			unbindgas(&HaloDestLen,&HaloDest,&Prop,SubCat.SubLen[subid],SubCat.PSubArr[subid],SubCat.Property[subid].CoM,SubCat.Property[subid].VCoM);	
			for(i=0;i<HaloDestLen;i++)
			{
				fprintf(fpdes,"%g\t%g\t%g\n",Gdat.Pos[HaloDest[i]][0],Gdat.Pos[HaloDest[i]][1],Gdat.Pos[HaloDest[i]][2]);
			}
			if(SubDestLen)
			{
			unbindgas(&SubDestLen,&SubDest,&Prop,SubCat.SubLen[subid],SubCat.PSubArr[subid],SubCat.Property[subid].CoM,SubCat.Property[subid].VCoM);	
			}
			for(i=0;i<SubDestLen;i++)
			{
				fprintf(fpsub,"%g\t%g\t%g\n",Gdat.Pos[SubDest[i]][0],Gdat.Pos[SubDest[i]][1],Gdat.Pos[SubDest[i]][2]);
			}
			hostid=SubCat.HaloChains[subid].HostID;
			dr=distance(SubCat.Property[SubCat.GrpOffset_Sub[hostid]].CoM,SubCat.Property[subid].CoM);
			rvir=comoving_virial_radius(SubCat.SubLen[SubCat.GrpOffset_Sub[hostid]]);
			fprintf(fp,"%d,%d,%d,%d,%d,%d,%d,",SnapBind,HaloDestLen,SubDestLen,BndLocLen,GasLocLen,SubCat.SubLen[subid],SubCat.SubLen[SubCat.GrpOffset_Sub[hostid]]);
			fprintf(fp,"%d,%d,%d,%g,%g\n",hostid,subid,SubCat.SubRank[subid],dr,dr/rvir);
			fflush(fp);
			fprintf(fpcen,"%g\t%g\t%g\n",SubCat.Property[SubCat.GrpOffset_Sub[hostid]].CoM[0],SubCat.Property[SubCat.GrpOffset_Sub[hostid]].CoM[1],SubCat.Property[SubCat.GrpOffset_Sub[hostid]].CoM[2]);
			fprintf(fpsubcen,"%g\t%g\t%g\n",SubCat.Property[subid].CoM[0],SubCat.Property[subid].CoM[1],SubCat.Property[subid].CoM[2]);
			/*==restore HaloDest==*/
			HaloDestLen=GCat.Len[fofid];
			free(HaloDest);
			HaloDest=mymalloc(sizeof(int)*HaloDestLen);
			memcpy(HaloDest,GCat.PIDorIndex+GCat.Offset[fofid],sizeof(int)*HaloDestLen);
			/*==fresh_index2ID==*/
			for(i=0;i<SubDestLen;i++)
				SubDest[i]=Gdat.PID[SubDest[i]];
			/*==free_sub_cat==*/
			for(i=0;i<SubCat.Nsubs;i++)
			{
				myfree(SubCat.PSubArr[i]);
			}
			free_sub_catalogue(&SubCat);
			/*==update subid for next snapshot==*/
			if(SnapBind<99)
			{
			load_pro2dest(SnapBind,&Pro2dest,&Npro,SUBCAT_DIR);
			subid=Pro2dest[subid];
			free_pro2dest(Pro2dest);
			if(subid<0) break;
			}
			fclose(fpdm);
			fclose(fpsub);
			fclose(fpdes);
			fclose(fpsrc);
			fclose(fpbnd);
			fclose(fploc);
		}
		fclose(fp);	
		fclose(fpcen);
		fclose(fpsubcen);
		return 0;
}
