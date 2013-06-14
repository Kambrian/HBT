typedef struct
{
   int Np;   //number of particles in the simulation
   int ips;  //time step number, or snapshot num
   float ztp; // current redshift
   float Omega0;  //omega_m0
   float OmegaLambda; //OmegaLambda0;   
   float Omegat;  // current omega_m, mass density
   float Lambdat; // current Omega_lambda, dark energy density
   float rLbox;  // boxsize in unit of Mpc/h
   float xscale;
   float vscale;
   float Hz; //Hubble Param at ztp
   float vunit; //velocity unit to get physical peculiar velocity in km/s
   /*==extension needed for calculating binding energy:==*/
   float time;//current reduced scale factor 
   float mass[2];//Gas (mass[0]) and DM (mass[1]) particle masses, in units of 10^10Msun/h
   int Nsnap; //current snapnum
}IO_HEADER; //header of jing's data structure

typedef struct 
{
int Ngroups;
int Nids;
int *Len;
int *Offset;
float *HaloCen[3];//HaloCen[3][Ngroups], center of halos
int *PIDorIndex; //stores PID when loading FOF then changes to PIndex after loading particles
short *HaloMask; //HaloMask[NP_DM],HaloMask==1 means the particle does not belong to any sub,i.e, it's free.
short *HaloMaskSrc;
int *ID2Halo;//for index2halo, this is ID2halo[NP_DM];
}CATALOGUE;

struct ParticleData
{
//#ifndef PID_ORDERED
int *PID; //since Jing's snap data is stored according to PID, so this array is not necessary here
//#endif
float (*Pos)[3];
float (*Vel)[3];
int Nsnap; //this exists to help check whether the Pdat is loaded for the current snapshot
};	

extern float PMass;
extern IO_HEADER header;//this header filled during load_particle_data();
extern struct ParticleData Pdat;
extern int snaplist[MaxSnap];

