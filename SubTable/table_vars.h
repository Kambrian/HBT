typedef MBDCAT //most-bound particle catalogue
{
	HBTInt MBD_PID;
	HBTInt HistID; //unique history-id in EvoCat
	HBTInt HostSubID;  //subid for the central-subhalo of its host halo
	HBTReal Pos[3];//position for MstbndID,comoving; 
	HBTReal Vel[3];//velocity for MstbndID,physical
}
