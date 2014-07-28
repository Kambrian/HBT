#ifndef PARAM_FILE_INCLUDED	//to avoid multiple inclusion
//use merger_io!
	
	/*=========program IO params==========*/
    #define ROOT_DIR "/gpfs/data/jvbq85/HBT/data/majormerger/static"
	#define GASCAT_DIR  ROOT_DIR "/gascat"			//the output directory for gashalocat,gassubcat and gassrccat, this must be an existing directory
	#define SUBCAT_DIR  ROOT_DIR "/subcat"			//the output directory for subcatalogues and srccatalogues, this must be an existing directory
	//necessary sub-dirs: splitters,pro2dest,history
	#define GRPCAT_DIR  ROOT_DIR "/fof"				//the input directory for GrpCatalogues
	#define  SNAPSHOT_DIR  ROOT_DIR "/simu"					//the input directory for simulation snapshots
	#define LOGFILE_NAME  "logfile"																			//the name of program logfile, set to "stdout" to use stdout

	/*======simulation params===========*/
	#define NP_SIM	3106774 //total Number of all kinds of particles
	#define NP_GAS 	0
	#define NP_DM 	3106774
	#define BOXSIZE 250.0
	//#define partmass=0.01
	#define MP_DM   0.01
	#define MP_GAS   0.  //10^10Msun/h; 
	
	/*=======Tree algorithm params========*/
	#define TREE_ALLOC_FACTOR 2//  MaxNodes = ((maxnodes>500)?maxnodes:500);
	#define ErrTolTheta 0.45
	#define SofteningHalo 2.1e-3
	#define NodeResolution 2.1e-4  //0.1*Softening
	#define NodeReso_half 1.05e-4
	
	/*=======unbinding algorithm params=====*/
	#define SAT_ACCR_ON      //define this to enable hierarchical accretion of satellite subs; comment this out to disable it
	#define PrecMass 0.995 //relative converge error for unbind (ErrorTolMass=1-PrecMass)
	#define NBOUNDMIN 20  //nbound<2 means only 1 particle bound, i.e, nothing left.
	#define NSRCMIN 20  //min len of src_sub (used for splitting)
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
	#define CoreLenMin 20//CoreLenMin<=NBOUNDMIN

	/*if GADGET units used, this gives the corresponding Gravitational constant; 
	 * otherwise calculate it in your own units and some modification to the unbinding procedure may be needed*/
//	#define CONVERT_LENGTH_MPC_KPC //this tells to convert unit from Mpc to kpc, 
	                   //because this is usually the only difference between standard gadget unit and user unit
	#define G 43.0071
	#define HUBBLE0 100.    //H_0 in internal units
	#define MaxSnap 12//total number of snapshot outputs
	#define SNAPFILE_BASE "snap"
//	#define SHIFT_PID_0_1  //if the PID in snapshot and groupcat starts from 0, define this to convert them to 1+
	#define NFILES 1   //number of files per snapshot
	#define NFILES_GRP 1 //number of group files per snapshot
	#define LOAD_SNAP_AS_GRP
	
#define PARAM_FILE_INCLUDED
#endif

