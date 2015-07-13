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
	SUBCATALOGUE SubCat;
	CATALOGUE Cat;
	CrossSection Sec;
	HBTInt *pro2dest, Npro,NsubLast;
	HBTInt Nsnap,SnapRange[2];
	
	parse_snap_args(SnapRange, argc, argv);	
	Nsnap=SnapRange[0]-1;
	load_cross_section(Nsnap, &Sec);//load previous cross-section
	if(Nsnap<IniSnap)
	  NsubLast=0;
	else
	{
	  load_sub_table(SnapRange[0]-1, &SubCat, SUBCAT_DIR);
	  NsubLast=SubCat.Nsubs;
	  free_sub_table(&SubCat);
	}
	
	for(Nsnap=SnapRange[0];Nsnap<SnapRange[1];Nsnap++)
	{
	  load_sub_catalogue(Nsnap, &SubCat, SUBCAT_DIR);
	  load_group_catalogue(Nsnap, &Cat, GRPCAT_DIR);
	  load_pro2dest(Nsnap-1, &pro2dest, &Npro, SUBCAT_DIR);
	  load_particle_data(Nsnap, SNAPSHOT_DIR);
	  fill_PIDHash();
	  fresh_ID2Index(&Cat,FRSH_GRPCAT); 
	  fresh_ID2Index(&SubCat,FRSH_SUBCAT);//now CatB.PIDorIndex has been occupied with Index to the freshly loaded Pdat
								//also refresh subPindex according to the new Pdat.
	  section_ID2Index(&Sec);
	  free_PIDHash();
	  prepare_ind2halo(&Cat);
	  
	  HBTInt NodeID,SubID;
	  #pragma omp parallel for
	  for(NodeID=0;NodeID<Sec.NumNode;NodeID++)
		node_update(Sec.Node+NodeID, pro2dest, &SubCat, &Cat);
	  NodeID=Sec.NumNode-1;//last node
	  #pragma omp parallel
	  #pragma omp single
	  for(SubID=0;SubID<SubCat.Nsubs;SubID++)
	  {
		if(SubCat.HaloChains[SubID].ProSubID<0||SubCat.HaloChains[SubID].ProSubID>=NsubLast) /*no progenitor or splitter;
																							*TODO: consider treating splitters differently */
		{
		  NodeID++;
		  if(NodeID==Sec.NumNodeAlloc)
		  {
			Sec.NumNodeAlloc*=2;
			Sec.Node=realloc(Sec.Node, sizeof(BranchNode)*Sec.NumNodeAlloc);
		  }//this is finished before any of the tasks below is executed, so they all see the updated Sec.
		  #pragma omp task firstprivate(NodeID,SubID)//each task will have its own NodeID and SubID at the time of creation
		  section_fill_new_node(&Sec, NodeID, SubID, &SubCat);
		}
	  }
	  NsubLast=SubCat.Nsubs;
	  section_Ind2ID(&Sec);
	  save_cross_section(Nsnap, &Sec);
	  free_catalogue(&Cat);
	  free_sub_catalogue(&SubCat);
	  free_pro2dest(pro2dest);
	  free_particle_data();
	}
	section_free(&Sec);
	return 0;
}