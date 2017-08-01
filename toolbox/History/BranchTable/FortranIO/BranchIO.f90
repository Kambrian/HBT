! For reading BranchTables
! =====Definitions/Data Organizations=======
! Branch: The evolution history of each subhalo, extending from birth time to z=0, is called a branch.
! 		if the subhalo is disrupted, then the most-bound particle is recorded.
!       each branch is labelled by a BranchID, from 0 to NBranch-1 where NBranch is the number of branches at a given snapshot.
! BranchNode: the subhalo information at a given snapshot inside the branch.
! BranchTable: the collection of all the branches.
! CrossSection: all the BranchNodes at given single snapshot in the BranchTable
!
! =====CAUTION: ============
!    whether the integers and reals are 4 or 8 bytes depends on the 
!    parameters for the simulation (in param/paramRUNNAME.h). need to manually 
!    adjust this for each run. check param/paramRUNNAME.h to see if 
!    HBT_INT8 and HBT_REAL8 are defined. 
! 	
!    The default integer types can be changed at compile time, with something like:
!    
!     gfortran -c -fdefault-integer-8 -fdefault-real-8 BranchIO.f90
!    or
!     ifort -c -i8 -r8 BranchIO.f90
!
!!

module BranchIO
implicit none

type BranchNode
  integer BranchID !ID of Branch, from 0 to NumNode-1. nodes with the same BranchID from different snapshot form a Branch.
  integer SubID !ID of subhalo, 0 to Nsub-1 at each snapshot; special case: SubID=-1 means disrupted subhalo.
  integer HostID !Host FoF ID, 0 to Ngroups; special case: HostID=-1 means no fof-host. (found in the background).
  integer SubRank !Rank of subhalo inside its host halo. 0: central; >0: satellite; -1: disrupted subhalo (only most-bound particle available)
  integer NpBnd !number of bound particles. =0 if unbound.
  integer NpBndPeak ! maximum value of NpBnd in all previous snapshots.
  integer SnapNumPeak ! the snapshot number (0 to MaxNumSnap) at NpBndPeak is found. 
  integer MstBndID ! most bound particle ID (as in the snapshot file; this is ID, not index.); special case: MstBndID=-1 if the subhalo has not been bound since birth upto the current snapshot, which means this is a fake object since the beginning.
  real MstBndPos(3) ! comoving coordinate of mostbound particle, kpc/h; =[0,0,0] if MstBndID=-1
  real MstBndVel(3) ! physical peculiar velocity of mostbound particle, km/s; =[0,0,0] if MstBndID=-1
  real Vmax !physical vmax, km/s;
  real VmaxPeak !peak of Vmax in all previous snapshots
  integer SnapNumVpeak ! the snapshot number at VmaxPeak
end type BranchNode

type CrossSection
  integer NumNode !number of nodes (Branches) in this CrossSection
  type(BranchNode),allocatable:: Node(:) ! list of nodes
end type CrossSection

contains
subroutine load_cross_section(Nsnap, Sec, SubCatPath)
  integer Nsnap
  type(CrossSection) Sec
  
  character(Len=*) SubCatPath
  character(Len=1024):: infile,snum
  integer NumNodeCheck
  
  write(snum,'(I3.3)') Nsnap
  infile=trim(SubCatPath)//'/CrossSection_'//trim(snum)
  open(unit=20,file=infile,form='unformatted',status='old',action='read',access='stream')
	
  read(20) Sec%NumNode
  allocate(Sec%Node(0:Sec%NumNode-1))
  read(20) Sec%Node
  read(20) NumNodeCheck !should be the same as Sec%NumNode.
  
  if(NumNodeCheck.ne.Sec%NumNode) then
	write(*,*) 'Error reading file '//infile
	write(*,*) 'File corruption or wrong datatype.'
	write(*,*) Sec%NumNode, NumNodeCheck
	close(20)
	stop 1
  end if
  
  close(20)
end subroutine load_cross_section

subroutine free_cross_section(Sec)
  type(CrossSection) Sec
  
  Sec%NumNode=0
  deallocate(Sec%Node)
end subroutine free_cross_section

end module BranchIO
