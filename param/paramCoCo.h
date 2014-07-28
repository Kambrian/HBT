#ifndef PARAM_FILE_INCLUDED	//to avoid multiple inclusion
	
	#define PERIODIC_BDR
#ifndef DISABLE_HALO_PARA
	#define HALO_PARA
#endif
	/*=========program IO params==========*/
#define ROOTDIR "gpfs/data/jvbq85/HBT/data/CoCo"
	#define SUBCAT_DIR  ROOTDIR"/subcat/"			//the output directory for subcatalogues and srccatalogues, this must be an existing directory
	//necessary sub-dirs: splitters,pro2dest,history
	#define GRPCAT_DIR  ROOTDIR"/simu/"				//the input directory for GrpCatalogues
	#define  SNAPSHOT_DIR  ROOTDIR"/simu/"					//the input directory for simulation snapshots
	#define LOGFILE_NAME  "logfile"																			//the name of program logfile, set to "stdout" to use stdout

	/*======simulation params===========*/
	#define NP_SIM	4794310656LL //total Number of all kinds of particles
	#define NP_GAS 	0
	#define NP_DM 	4286872576LL
	#define BOXSIZE  70.4  //Mpc/h
	#define BOXHALF 35.2
	#define MP_DM    1.135728e-05
	#define MP_GAS   0.

	/*=======Tree algorithm params========*/
	#define TREE_ALLOC_FACTOR 1//  MaxNodes = ((maxnodes>500)?maxnodes:500);
	#define ErrTolTheta 0.45
	#define SofteningHalo  2.3e-4 //Mpc/h
	#define NodeResolution 2.3e-5 //0.1*Softening
	#define NodeReso_half 1.15e-5
	
	/*=======unbinding algorithm params=====*/
	//~ #define E_Relax 9  //(E_Relax+1)*pot+K<0
	#define SAT_ACCR_ON      //define this to enable hierarchical accretion of satellite subs; comment this out to disable it
	#define PrecMass 0.99 //relative converge error for unbind (ErrorTolMass=1-PrecMass); 0.995 by default.
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
	#define G 43.0071 //Mpc/h
	#define HUBBLE0 100    //H_0 in internal units
	#define MaxSnap 160//total number of snapshot outputs
	#define SNAPFILE_BASE "snap_coco_big_run_1-groupordered"
//	#define SNAPLIST 1131 //snaplist identifier to enable the corresponding snaplist in iovars.c
	#define NFILES 512   //number of files per snapshot
	#define NFILES_GRP 512 //number of group files per snapshot
// 	#define INPUT_REAL8   //datatype for input in double precision
	#define INPUT_INT8    //datatype for input ID
// 	#define HBT_REAL8    //datatype for HBT calculation and output
	#define HBT_INT8      //datatype of HBT integers, must be able to hold all PIDs
	//~ #define HBTPID_RANKSTYLE //replace PIDs with their ranks at input(in the range [0~NP_DM-1])
                             //in this case HBT_INT can be smaller than INPUT_INT if the original pids are not on ground state
							//these ground state pids will also be saved into subcat files
	#define PID_NEED_HASH //only needed if index table is not efficient, e.g., when PIDs are very large or discontinous, 
							//so that the index table consumes incredible amount of memory
	#define GRP_V4FORMAT  
	
#if defined(GRP_V4FORMAT)&&!defined(GRPINPUT_INDEX)
  #define GRPINPUT_INDEX
#endif

#define PARAM_FILE_INCLUDED
#endif

