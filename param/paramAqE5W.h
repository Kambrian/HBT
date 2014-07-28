#ifndef PARAM_FILE_INCLUDED	//to avoid multiple inclusion
		//once named W200
	/*=========program IO params==========*/
	#define GASCAT_DIR  "/gpfs/data/jvbq85/HBT/data/AqE5W/gascat"			//the output directory for gashalocat,gassubcat and gassrccat, this must be an existing directory
	#define SUBCAT_DIR  "/gpfs/data/jvbq85/HBT/data/AqE5W/subcat"			//the output directory for subcatalogues and srccatalogues, this must be an existing directory
	//necessary sub-dirs: splitters,pro2dest,history
	#define GRPCAT_DIR  "/gpfs/data/jvbq85/HBT/data/AqE5W/fof"				//the input directory for GrpCatalogues
	#define SNAPSHOT_DIR  "/gpfs/data/jvbq85/HBT/data/AqE5W/simu"					//the input directory for simulation snapshots
	#define LOGFILE_NAME  "logfile"																			//the name of program logfile, set to "stdout" to use stdout

	/*======simulation params===========*/
	#define NP_SIM	2951478 //total Number of all kinds of particles
	#define NP_gas 	0
	#define NP_DM 	2317091
	#define BOXSIZE 	100000.0  //kpc/h
	//#define partmass=2.29431e6  Msun/h
	#define MP_DM   2.29431e-4
	#define MP_GAS   0.  //10^10Msun/h; 

	/*=======Tree algorithm params========*/
	#define TREE_ALLOC_FACTOR 2//  MaxNodes = ((maxnodes>500)?maxnodes:500);
	#define ErrTolTheta 0.45
	#define SofteningHalo 0.85
	#define NodeResolution 0.085  //0.1*Softening
	#define NodeReso_half 0.0425
	
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
	#define CONVERT_LENGTH_MPC_KPC //this tells to convert unit from Mpc to kpc, 
	                   //because this is usually the only difference between standard gadget unit and user unit
	#define G 43007.1
	#define HUBBLE0 0.1    //H_0 in internal units
	#define MaxSnap 128//total number of snapshot outputs
	#define SNAPFILE_BASE "snap"
	#define NFILES 1   //number of files per snapshot
	#define NFILES_GRP 1 //number of group files per snapshot
//	#define INPUT_DOUBLEPRECISION   //datatype for input in double precision
//	#define OUTPUT_DOUBLEPRECISION   //datatype for output in double precision
//	#define LONG_ID                //long int for particle id
//	#define GRP_V3FORMAT           //default using PGADGET-3's group(subfind) format
#define PARAM_FILE_INCLUDED
#endif

