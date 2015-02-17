typedef struct
{
	HBTInt Mdm;
	HBTInt HostID;
	HBTInt SubRank;
	HBTInt ProID;
	HBTInt DesID;
} SubData;
typedef struct
{
	HBTInt Nsnap;
	HBTInt Ngroups;
	HBTInt Nsubs;
	SubData *SubHalo;
} BRFCAT;

typedef struct 
{
	HBTInt Mdm;//sub DM mass
	HBTInt SubID;
	HBTInt HostID;//host haloid
	HBTInt SubRank;
} SubNodePre;
typedef struct
{
	HBTInt SnapBirth;
	HBTInt SnapDeath;//when subhalo becomes under resolution
	HBTInt SnapEnter;
	HBTInt ProHistID;// -1 means the same as itself; positive integer means this is a splitter, 
					//so its progenitor is in a different history,given by ProHistID
	SubNodePre *Member;
} HISTORY_Pre;
typedef struct
{
	HBTInt NNode;
	HBTInt NHist;
	HBTInt *HistLen;
	HBTInt *HistOffset;
	HISTORY_Pre *History;
}EVOLUTIONCAT_Pre;

typedef struct 
{
	HBTInt Mdm;//sub DM mass
	HBTInt Mhost;//virial mass if possible,otherwise fof-mass
	HBTReal Chost;//-1 if not NFW-fittable
	HBTInt SubID;
	HBTInt HostID;//host haloid,now dynamically extended in EvoCatRev
	HBTInt SubRank;
} SubNode;
typedef struct
{
	HBTInt SnapBirth;
	HBTInt SnapDeath;//when subhalo becomes under resolution
	HBTInt SnapEnter;
	HBTInt ProHistID;// -1 means the same as itself; positive integer means this is a splitter, 
					//so its progenitor is in a different history,given by ProHistID
	SubNode *Member;
} HISTORY;
typedef struct
{
	HBTReal PartMass;
	HBTInt NNode;
	HBTInt NHist;
	HBTInt *HistLen;
	HBTInt *HistOffset;
	HISTORY *History;
}EVOLUTIONCAT;

typedef struct
{
	HBTInt HostID;
	HBTInt subid;
	HBTInt Mdm;
	HBTInt Mhost;//host halo mass
	HBTInt Mcen;//host subhalo mass
	HBTInt Msub;//subhalo total dm mass (including sub-in-subs);
#ifdef GAS_TRACE	
	HBTInt Mgas;
	HBTInt Mhostgas;//host halo gas mass
	HBTReal Umsat;//sat average thermal energy (Temperature)
#endif
	HBTReal CoM[3];//relative position,comoving
	HBTReal VCoM[3];//relative Vel,physical
	HBTReal Chost;//host halo concentration
}SubNodeExp;

//typedef HBTReal TYPE_XYZ[3]; //replaced by HBTxyz in datatypes.h
typedef struct
{
	int mass;//fof mass
	int Mvir[3];
	float Rvir[3];//[tophat,c200,b200],comoving
	int flag_badvir[3];
	int flag_fakehalo;//set to 1 when halo is not self-bound,to 0 otherwise.
} HALOSIZE;
struct ShardParam
{
	HBTInt SnapBirth;
	HBTInt SnapTidal;
	HBTInt SnapRvir; //the first snapshot when D<Rvir; -1 means no virial-crossing
	HBTReal Mrate[2];//Msat/Mhost, interpolated at D=RTidal and D=Rvir
	HBTReal Mhost[2];//host virial mass, interpolated at RTidal and Rvir
	HBTReal Rhost[2];
	HBTReal Kappa[2]; //energy parameter, V^2/Vcirc^2, interpolated at RTidal and Rvir
	HBTReal j2[2];  //angular momentum parameter, circularity^2
	HBTReal ErrK[2];//error estimation due to discrete output
	HBTReal Errj2[2];
	HBTReal Chost[2];//host concentration
	HBTReal Csat[2];//sat concentration
};
struct HistoryShards
{
	HBTInt NumHist;
	HBTInt NumShards;//number of mergers (infalls) ever happened
	HBTInt *NBirth;
	HBTInt *HistoryOffset;
	struct ShardParam **Par;
};


typedef struct //most-bound particle 
{
	HBTInt MBD_PID;
	HBTInt HistID; //unique history-id in EvoCat
	HBTInt HostID; //host haloid
	HBTReal Pos[3];//position for MstbndID,comoving; 
	HBTReal Vel[3];//velocity for MstbndID,physical
} MBDNode;

typedef struct //most-bound particle catalogue
{
	HBTInt NSubs;  //only living subs (same as in subcat); total: Nsubs+NOrphans.
	HBTInt NOrphans;  //only orphans  
	HBTInt Ngroups;
	HBTInt NQuasi; //subhalos without a host
	HBTInt NOrphanQuasi;  //quasi without a host
	HBTInt *GrpLen_Sub;     //same as in subcat
	HBTInt *GrpOffset_Sub;  //same as in subcat
	HBTInt *GrpLen_Orphan; //number of orphans for each group
	HBTInt *GrpOffset_Orphan; //the starting position for the orphans of each group, in *Nodes
	MBDNode *Nodes;
}MBDCATALOGUE;
