#include "datatypes.h"

typedef struct
{
  HBTInt BranchID; //id of branch, starting from 0.
  HBTInt SubID; 
  HBTInt HostID;
  HBTInt SubRank; //0: central; >0: satellite.
  HBTInt NpBnd; //number of bound particles
  HBTInt NpBndPeak; //maximum value of NBound in all previous snapshots
  HBTInt SnapNumPeak; //snapshot number at which NBoundPeak is found.
  HBTInt MstBndID; //ID of most bound particle
  HBTxyz MstBndPos; //comoving coordinate of MostBound particle
  HBTxyz MstBndVel; //physical velocity of MostBound particle
} BranchNode;

typedef struct
{
  HBTInt NumNode;
  HBTInt NumNodeAlloc; //allocated number of nodes.
  BranchNode *Node;
} CrossSection;

#define NumNodeAllocMin 1024
extern void section_ID2Index(CrossSection *sec);
extern void section_Ind2ID(CrossSection *sec);
extern void section_init(CrossSection * sec);
extern void section_free(CrossSection * sec);
extern void save_cross_section(HBTInt Nsnap, CrossSection *sec);
extern void load_cross_section(HBTInt Nsnap, CrossSection * sec);
extern void section_fill_new_node(CrossSection *sec, HBTInt NodeID, HBTInt SubID, SUBCATALOGUE *SubCat);
extern void node_update(BranchNode *node, HBTInt *pro2dest, SUBCATALOGUE *SubCat, CATALOGUE *Cat);
extern void parse_snap_args(HBTInt SnapRange[2], int argc, char **argv);