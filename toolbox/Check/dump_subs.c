#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define SUBFIND_DIR "/gpfs/data/jch/MilliMillennium/Snapshots/"
extern void load_subfind_catalogue(int Nsnap,SUBCATALOGUE *SubCat,char *inputdir);	
#define load_sub_catalogue load_subfind_catalogue

int main(int argc, char** argv)
{
	CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	HBTInt Nsubs,Nhalo;
	
	HBTInt Nsnap=0;
	HBTInt grpid,subid,pid;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=2)
	{printf("usage: %s [Nsnap], otherwise Nsnap=%d\n",argv[0],Nsnap);fflush(stdout);}
	else
	Nsnap=atoi(argv[1]);
	
char buf[1024];
#ifdef SUBFIND_DIR	
 load_sub_catalogue(Nsnap,&SubCat,SUBFIND_DIR);
 sprintf(buf,"%s/anal/subsnap_%d.subfind",SUBCAT_DIR,Nsnap);	
#else
 load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
 sprintf(buf,"%s/anal/subsnap_%d",SUBCAT_DIR,Nsnap);	
#endif
 /*
 load_particle_data(Nsnap, SNAPSHOT_DIR);
 fill_PIDHash();
 fresh_ID2Index(&SubCat, FRSH_SUBCAT);
 free_PIDHash();
 */
  FILE *fp;
  myfopen(fp,buf,"w");
  fprintf(fp,"Rank\tNBound\tX\tY\tZ\tHostId\n");
  HBTxyz NullPos={0,0,0};
  for(subid=0;subid<SubCat.Nsubs;subid++)
  {
	HBTReal *pos, *vel;
	if(SubCat.SubLen[subid]) 
        {
// 	  pos=Pdat.Pos[SubCat.PSubArr[subid][0]];
	  pos=SubCat.Property[subid].CoM;
          vel=SubCat.Property[subid].VCoM;
        }
	else
        {
	  pos=NullPos;
          vel=NullPos;
        }
	fprintf(fp,"%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n",SubCat.SubRank[subid],SubCat.SubLen[subid],
		pos[0],pos[1],pos[2], vel[0], vel[1], vel[2], SubCat.HaloChains[subid].HostID);
  }
  fclose(fp);
 
  sprintf(buf, "%s/anal/halosize_%d.subfind", SUBCAT_DIR, Nsnap);
  myfopen(fp, buf, "w")
  fprintf(fp, "CenSubId\t Mtophat\t Rtophat\t Mcrit\t Rcit\t Mmean\t Rmean\n");
  float (*Mvir)[3]=mymalloc(sizeof(float)*SubCat.Ngroups*3);
  float (*Rvir)[3]=mymalloc(sizeof(float)*SubCat.Ngroups*3);
  printf("%d\n", SubCat.Ngroups);
  load_subfind_halo_size(Mvir, Rvir, 0., SubCat.Ngroups, Nsnap);
  for(grpid=0;grpid<SubCat.Ngroups;grpid++)
  {
      if(SubCat.GrpLen_Sub[grpid])
      fprintf(fp, "%d\t %g\t %g\t %g\t %g\t %g\t %g\n", SubCat.GrpOffset_Sub[grpid], Mvir[grpid][0], Rvir[grpid][0], Mvir[grpid][1], Rvir[grpid][1], Mvir[grpid][2], Rvir[grpid][2]);
      else
      fprintf(fp, "0\t0\t0\t0\t0\t0\t0\n");
  }
  myfree(Mvir);
  myfree(Rvir);
  
	free_sub_catalogue(&SubCat);
// 	free_particle_data();
	return 0;
}
