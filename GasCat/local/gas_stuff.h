#ifndef GAS_STUFF_INCLUDE

#define EXCLUSIVE_GAS
 //~ #define THERMAL_BOUND
 #ifndef SofteningGas
 	#define SofteningGas 1.5
#endif
#define GSub_ScaleRelax 1.5
#define NDIV 200

typedef struct
{
	int Nsubs;
	int Nids;
	int *SubLen;
	int *SubOffset;
	int **PSubArr;//PSubArr[Nsubs] as a pointer will point to PIDs block (or PIndices during unbind()) of each sub;i.e.,PSubArr[subhaloid][particle] will give PIDs(or PInd);
} GASCATALOGUE;

 struct GasData
{
int PID[NP_gas];
float Pos[NP_gas][3];
float Vel[NP_gas][3];
float U[NP_gas];
float Rho[NP_gas];
}Gdat;

extern int hoc[NDIV][NDIV][NDIV],ll[NP_gas];
extern float range[3][2], step[3];

extern float *Rtidal;
extern int *GHostSub;//initialized during load_gas_data
extern void load_gas_data(int Nsnap, char *SnapPath);
//we do not resolve to PID of each gas particle, rather we use each particle's index in the snapshot directly to refer it, and to get its properties, 
//and to store in gascat, so the same particle index in each gascat do not correspond to the same particle
extern double tree_treeevaluate_potential_gas(int target, int *PIndex);	
extern int collect_gas_particles(int subid,SUBCATALOGUE*SubCat,GASCATALOGUE *GasCat);
extern int unbindgas(int *P2Len,int **P2PIndex, int SubLen, int *SubArr); /*P2Len=&Len, *P2PIndex=PIndex, 
																							*where PIndex[Len] is the array of size Len; 
																							* both will be updated after unbind*/
extern void unbind_gas_recursive(int mainsubID,SUBCATALOGUE *SubCat,GASCATALOGUE *GasCat);
extern void create_gas_cat(int Nsubs,GASCATALOGUE *GasCat);
extern void save_gas_cat(int Nsnap,GASCATALOGUE *Cat,char *gasdir);
extern void load_gas_cat(int Nsnap,GASCATALOGUE *Cat,char *subdir);
extern void free_gas_cat(GASCATALOGUE *Cat);
extern void load_tidal_radius(int Nsnap,float *rtidal, int Nsubs,char*tidaldir);
extern void makell(float (*pos)[3],int np);

#define GAS_STUFF_INCLUDE
#endif
