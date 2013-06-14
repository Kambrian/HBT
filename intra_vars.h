/*hash table data for ID2Index*/
struct ID2Ind
{
	HBTInt PID;
	HBTInt PInd;  //address in Pdat.PID
};
struct PIDHASH
{
HBTInt np;	
struct ID2Ind *table; //general ID2Ind pairs in case hash is really needed
HBTInt *PIndex;  //a simple Index array is sufficient if max(PID) is small enough
HBTInt PInd_offset;
};

struct Chain_data 
{
HBTInt ProSubID;
HBTInt HostID;
};

struct LinkInfo
{
	HBTInt *ProCount;
	HBTInt *ChainOffset;
	HBTInt *pro2dest;
};
struct Hierarchy
{
HBTInt nibs;	//its host sub
HBTInt pre;	//its previous sibling
HBTInt next;	//its next sibling
HBTInt sub;	//its sub-in-sub
} ;
struct SubProperty
{
	HBTReal CoM[3];//center of mass position (comoving)
	HBTReal VCoM[3];//CoM velocity (physical)
	HBTReal Pot;//average gavitational potential per unit mass per particle (GM/R_physical)
	HBTReal Kin;//average Kinetic energy per unit mass per particle (0.5*V_physical^2), or 0.5*VelocityDispersion
	HBTReal AM[3];//average angular momentum per unit mass per particle (R_physical x V_physical)
};
typedef struct
{
	HBTInt Ngroups;
	HBTInt Nsubs;
	HBTInt Nids;
	HBTInt *GrpLen_Sub;//number of subs in a group
	HBTInt *GrpOffset_Sub;//sub index offset of a group
	HBTInt *SubLen;
	HBTInt *SubOffset;

	HBTInt *SubRank;
	struct Chain_data *HaloChains;
	
	struct Hierarchy *sub_hierarchy;
	
	struct SubProperty *Property;   
	
	HBTInt **PSubArr;//PSubArr[Nsubs] as a pointer will point to PIDs block (or PIndices during unbind()) of each sub;i.e.,PSubArr[subhaloid][particle] will give PIDs(or PInd);
	
	HBTInt Nbirth;
	HBTInt NQuasi;
	HBTInt Ndeath;
	HBTInt Nsplitter;//the splitted-out subs from last snap
} SUBCATALOGUE;

typedef struct
{
	HBTInt Nsubs;
	HBTInt Nids;
	HBTInt *SubLen;
	HBTInt *SubLen2;
	HBTReal *CoreFrac;
	HBTInt *SubOffset;//? necessary?
	HBTInt **PSubArr;//PSubArr[Nsubs] as a pointer will point to PIDs block (or PIndices during unbind()) of each sub;i.e.,PSubArr[subhaloid][particle] will give PIDs(or PInd);
	HBTInt **PSubArr2;
	struct Chain_data *HaloChains;//this is only create during splitt_srccat() and then passed away to subcat during make_srcsub(); at any other times it remains NULL;
																			//add it here only to make it act like a global var; 
	HBTInt NDeathSp;//splitted to death
}SRCCATALOGUE;

struct cand_data
{
	HBTInt *SubArr;
	HBTInt SubLen;
	HBTInt desID;
};
struct Energy
  {
	  HBTInt PID;
	  HBTReal Erg;
  } ;	//energy array used for sorting Energy;
union NODE
{
	HBTInt sons[8];		/*!< temporary pointers to daughter nodes */
	struct
	{
	HBTReal s[3];               /*!< center of mass of node */
	HBTReal len;		/*!< sidelength of treenode */
	HBTInt mass;            /*!< mass of node */
	//     HBTInt cost;            /*!< counts the number of interactions in which this node is used */
	HBTInt sibling;         /*!< this gives the next node in the walk in case the current node can be used */
	HBTInt nextnode;        /*!< this gives the next node in case the current node needs to be opened */
	}
	way;
};

extern FILE *logfile;// logfile=stdout; if needed.

extern struct PIDHASH PIDHash;

/*===this part used for tree algorithm and binding energy calculation and sorting===*/
/*===they are applied as single halo var==*/
extern union NODE *Nodes_base,                    /*!< points to the actual memory allocted for the nodes */
  							*Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				  gives the first allocated node */

//HBTReal MaxNodeFilledFraction=0.;
extern size_t MaxNodes,MaxNids;
extern HBTInt NumPart;
extern HBTInt * Nextnode;	/*this is for particle node, for each single halo*/
#ifdef HALO_PARA
#pragma omp threadprivate(Nodes_base,Nodes,MaxNodes,MaxNids,NumPart,Nextnode)
#endif

//data structure for group-finding
struct ParticleGroup
{
	HBTInt PIndex;
	HBTInt GrpID;
};
struct GroupData
{
	HBTInt Np;
	HBTInt Ngrp;
	HBTInt *GrpLen;
	struct ParticleGroup *GrpTags;
};

typedef HBTReal * AccessPosFunc(HBTInt pid,void *PosData);
typedef struct 
{
	HBTInt ndiv;
	HBTInt np;
	HBTInt UseFullBox; //whether to use a local enclosing-box or full simulation box (latter for periodic boudary)
	HBTInt *hoc;
	HBTInt *list;
	HBTReal range[3][2];
	HBTReal step[3];
	void *PosData; //to be passed to GetPos() to access particle positions
	AccessPosFunc *GetPos;
} LINKLIST;


