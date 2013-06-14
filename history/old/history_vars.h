typedef struct 
{
	int Mdm;//sub DM mass
	int Mgas;//bound gas mass
	int Mhost;//host halo DM mass
	float Chost;//host halo concentration
	int SubID;
	int HostID;//host haloid
	int SubRank;
} SubNode;
typedef struct
{
	int SnapBirth;
	int SnapEnter;
	int ProHistID;// -1 means the same as itself; positive integer means this is a splitter, 
					//so its progenitor is in a different history,given by ProHistID
	SubNode Member[MaxSnap];
 }SubHist;
