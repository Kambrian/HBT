//last infall
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
#include "BranchFuncs.h"

extern void parse_snap_args(HBTInt SnapRange[2], int argc, char **argv);
int main(int argc,char **argv)
{
	CrossSection SecLast, Sec;
	HBTInt Nsnap,SnapRange[2];
        char olddir[1024], newdir[1024]="/mnt/ddnfs/jxhan/6610BranchTable";
        sprintf(olddir, "%s/BranchTable", SUBCAT_DIR);
	
	logfile=stdout;
	parse_snap_args(SnapRange, argc, argv);	
	Nsnap=SnapRange[0]-1;
	load_cross_section(Nsnap, &SecLast, newdir);//load previous cross-section
	
	for(Nsnap=SnapRange[0];Nsnap<=SnapRange[1];Nsnap++)
	{
          load_cross_section(Nsnap, &Sec, olddir);
	  load_particle_header(Nsnap, SNAPSHOT_DIR);
          HBTInt NodeId;
#pragma omp parallel for
          for(NodeId=0;NodeId<SecLast.NumNode;NodeId++)
          {
              BranchNode *node=Sec.Node+NodeId;
              node->Vmax/=header.time;//to fix the bug
              node->VmaxPeak=SecLast.Node[NodeId].VmaxPeak;
              node->SnapNumVpeak=SecLast.Node[NodeId].SnapNumVpeak;
              if(node->VmaxPeak<node->Vmax)
              {
                  node->VmaxPeak=node->Vmax;
                  node->SnapNumVpeak=Nsnap;
              }
          }
#pragma omp parallel for
          for(NodeId=SecLast.NumNode;NodeId<Sec.NumNode;NodeId++)
          {
              BranchNode *node=Sec.Node+NodeId;
              node->Vmax/=header.time;
              node->VmaxPeak=node->Vmax;
          }
              
	  save_cross_section(Nsnap, &Sec, newdir);
          section_free(&SecLast);
          SecLast=Sec;
	}
	section_free(&Sec);
	return 0;
}