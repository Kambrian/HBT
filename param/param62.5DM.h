#ifndef PARAM_FILE_INCLUDED	//to avoid multiple inclusion

//SubLenMin=20
	//#define USE_AHF_IO //to use ahf_io in the analysis routines
	#define PERIODIC_BDR
#ifndef DISABLE_HALO_PARA
	#define HALO_PARA
#endif
	/*=========program IO params==========*/
	#define ROOT_DIR "/gpfs/data/jvbq85/HBT/data/62.5DM/"
//	#define ROOT_DIR "/home/sussing/Working/HBT/DATASET_I/"
//	#define ROOT_DIR "/work/Projects/HBT/code/sussing2013/DATA/DATASET_I/"			//the output directory for gashalocat,gassubcat and gassrccat, this must be an existing directory
//	#define SUBCAT_DIR  ROOT_DIR "/subcatV3"			//the output directory for subcatalogues and srccatalogues, this must be an existing directory
	#define SUBCAT_DIR ROOT_DIR "subcatRaw"
	//necessary sub-dirs: splitters,pro2dest,history
	#define GRPCAT_DIR  ROOT_DIR "fof"				//the input directory for GrpCatalogues
	#define  SNAPSHOT_DIR  ROOT_DIR "simu"					//the input directory for simulation snapshots
	#define LOGFILE_NAME  "logfile"																			//the name of program logfile, set to "stdout" to use stdout

	/*======simulation params===========*/
	#define NP_SIM	19683000LL //total Number of all kinds of particles
	#define NP_GAS 	0
	#define NP_DM  NP_SIM
	#define BOXSIZE 	62.5e3  //kpc/h
	#define BOXHALF     31.25e3
	#define MP_DM  9.36395e-2
	#define MP_GAS   0.
	#define OMEGA0 0.272
	#define OMEGAL0 0.728

	/*=======Tree algorithm params========*/
	#define TREE_ALLOC_FACTOR 4//  MaxNodes = ((maxnodes>500)?maxnodes:500);
	#define ErrTolTheta 0.45
	#define SofteningHalo 5  
	#define NodeResolution 0.5  //0.1*Softening
	#define NodeReso_half 0.25
	
	/*=======unbinding algorithm params=====*/
	#define SAT_ACCR_ON      //define this to enable hierarchical accretion of satellite subs; comment this out to disable it
	#define PrecMass 0.995 //relative converge error for unbind (ErrorTolMass=1-PrecMass)
	#define NBOUNDMIN 20  //nbound<2 means only 1 particle bound, i.e, nothing left.
	#define NSRCMIN NBOUNDMIN  //min len of src_sub (used for splitting)
	#define NParaMin 100 //the lower bound of particle numbers for parallelization of unbind()
	#define MassRelax_Input 2/*control the dynamic range of the unbinding procedure; 
															* it is the maximum allowed ratio that a sattellite sub is allowed to reaccrete/grow, 
															* its square (MassRelax_Input*MassRelax_Input) also "roughly"(maybe either bigger or smaller) gives the maximum safe-unbinding  factor 
															* by which  a sub is stripped  ,i.e. if a sub reduces its mass by a factor of N>MRelax*MRelax between two subsequent snapshot, 
															* then the unbinding  procedure may not find the bound part safely (since CoreFrac=1.0/MRelax^2 for algorithm efficiency reason)
															* so it seems the bigger the better. 
															* However, bigger MRelax means smaller CoreFrac, then small CoreLen if the sub in hand is small,
															* thus may also undermine the unbinding process for small subs
															* But for big enough subs, bigger MRelax (thus smaller CoreFrac) usually makes unbinding more robust to disturbing particles */
	#define CoreFrac0 0.25 //1.0/MassRelax_Input/MassRelax_Input;
	#define CoreLenMin NBOUNDMIN//CoreLenMin<=NBOUNDMIN

	/*if GADGET units used, this gives the corresponding Gravitational constant; 
	 * otherwise calculate it in your own units and some modification to the unbinding procedure may be needed*/
	#define CONVERT_LENGTH_MPC_KPC //this tells to convert unit from Mpc to kpc, 
	                   //because this is usually the only difference between standard gadget unit and user unit
	#define G 43007.1
	#define HUBBLE0 0.1    //H_0 in internal units
	#define IniSnap 0
	#define MaxSnap 62//total number of snapshot outputs
	#define SNAPFILE_BASE "62.5_dm"
	#define NFILES_GRP 32 //number of group files per snapshot
	//#define GRP_HBTFORMAT  //fof produced from HBT
	#define NFILES 16
	#define HBT_INT8
	#define HBT_REAL8
	#define GRP_V3FORMAT
//	#define GRPINPUT_INDEX
//	#define VEL_INPUT_PHYSICAL  //this macro allows for the input velocity from load_particle_data() to be physical rather than in GADGET manner,i.e, do not need to multiply sqrt(a)
	#define IDOFFSET_UNIT 1000000000000LL  //to deal with AHF haloid
	
//	#define CLEAR_DUPLICATE_AHF_PARTICLES

#define PARAM_FILE_INCLUDED
#endif

