#include <stdio.h>
#include "datatypes.h"
#include "intra_vars.h"

FILE *logfile;// logfile=stdout; if needed.

struct PIDHASH PIDHash;

/*===this part used for tree algorithm and binding energy calculation and sorting===*/
/*===they are applied as single halo var==*/
union NODE *Nodes_base,                    /*!< points to the actual memory allocted for the nodes */
  			   *Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				                                                     gives the first allocated node */

//float MaxNodeFilledFraction=0.;
size_t MaxNodes,MaxNids;
HBTInt NumPart;
HBTInt * Nextnode;	/*this is for particle node, for each single halo*/
//~ #pragma omp threadprivate(Nodes_base,Nodes,MaxNodes,MaxNids,NumPart,Nextnode)
