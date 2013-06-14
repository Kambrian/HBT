struct GasData
{
HBTInt *PID;
HBTReal (*Pos)[3];
HBTReal (*Vel)[3];
HBTReal *U;
HBTReal *Rho;
};
struct GasProperty
{
	HBTReal CoM[3];//center of mass position (comoving)
	HBTReal VCoM[3];//CoM velocity (physical)
	HBTReal Pot;//average gavitational potential per unit mass per particle (GM/R_physical)
	HBTReal Kin;//average Kinetic energy per unit mass per particle (0.5*V_physical^2),w.r.t. DM center
	HBTReal AM[3];//average angular momentum per unit mass per particle (R_physical x V_physical),w.r.t. DM center
	HBTReal U;//average thermal energy per mass per particle
};
typedef CATALOGUE GASHALOCAT;
typedef struct
{
	HBTInt Nsubs;
	HBTInt Nids;
	HBTInt *SubLen;
	HBTInt *SubOffset;
	HBTInt **PSubArr;//PSubArr[Nsubs] as a pointer will point to PIDs block (or PIndices during unbind()) of each sub;i.e.,PSubArr[subhaloid][particle] will give PIDs(or PInd);
} GASSRCCAT;
typedef struct
{
	HBTInt Nsubs;
	HBTInt Nids;
	HBTInt *SubLen;
	HBTInt *SubOffset;
	HBTInt **PSubArr;//PSubArr[Nsubs] as a pointer will point to PIDs block (or PIndices during unbind()) of each sub;i.e.,PSubArr[subhaloid][particle] will give PIDs(or PInd);
	struct GasProperty *Property;
} GASSUBCAT;
extern struct GasData Gdat;
#define GasCollect_ScaleRelax 3
