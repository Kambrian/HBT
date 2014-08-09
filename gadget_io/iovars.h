typedef	struct 
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned npartTotal[6];  //differ from standard. to be able to hold large integers
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  double Hz;//current Hubble param in internal units, BT extension
  int Nsnap; //current snapshot
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8 -8 - 4];  /* fills to 256 Bytes */
} IO_HEADER;

typedef struct 
{
HBTInt Ngroups;
HBTInt Nids;
HBTInt *Len;
HBTInt *Offset;
#ifdef GRP_JINGFORMAT
HBTReal *HaloCen[3];//HaloCen[3][Ngroups], center of halos
#endif
HBTInt *PIDorIndex; //stores PID when loading FOF then changes to PIndex after loading particles
char *HaloMask; //HaloMask[NP_DM],HaloMask==1 means the particle does not belong to any sub,i.e, it's free.
char *HaloMaskSrc;
HBTInt *ID2Halo;//for index2halo, this is ID2halo[NP_DM];
}CATALOGUE;

struct ParticleData
{
#ifndef PID_ORDERED
HBTInt *PID; //since Jing's snap data is stored according to PID, so this array is not necessary here
#endif
HBTxyz *Pos;
HBTxyz *Vel;
HBTInt Nsnap; //this exists to help check whether the Pdat is loaded for the current snapshot
};	

extern IO_HEADER header;//this header filled during load_particle_data();
extern struct ParticleData Pdat;

#ifdef SNAPLIST
extern HBTInt	snaplist[MaxSnap];
#endif
