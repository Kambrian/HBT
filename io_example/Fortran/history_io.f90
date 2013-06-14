! For reading Evolution catalogue, which divides all the subhalos ever alive into many histories 
! each history consists of all the nodes of a subhalo at different snapshots
! each node record the properties of that subhalo at that snapshot 
! there're two types of evolution catalogue: preliminary and revised
! * the preliminary evolution catalogue (EVOLUTIONCAT_PRE)
!   is constructed using the fof-halo hosting information,as directly given in the subhalo catalogue.
!   better use this catalogue for getting halo formation time.
! * the revised evolution catalogue  (EVOLUTIONCAT)
!   #extends the host halos of a satellite subhalo to the time before infall and after crossing-out,
!    i.e., a subhalo's host is always set to be the one it falls into even at the time it's still outside the fof region of that host.
!   #also the host halos are switched according to R/Rvir when a subhalo is ejected from one halo and then fall into another
!   #each history is break into several pieces(Shards) according to host-crossing
!   #the merger parameters of each shard of history are given in ShardParam
!!

module history_io
implicit none
type SubNodePre
	integer*4 Mdm  !sub DM mass
	integer*4 SubID
	integer*4 HostID !host haloid
	integer*4 SubRank
end type SubNodePre;   ! primordial subhalo properties

type EVOLUTIONCAT_Pre
	integer*4 NNode   !number of Nodes (subhalos)
	integer*4 NHist   !number of Histories
	integer*4,allocatable::HistLen(:)    !Length of each history
	integer*4,allocatable::HistOffset(:)  !offset of each history
	integer*4,allocatable::SnapBirth(:)   !the first snapshot when this subhalo exist
	integer*4,allocatable::SnapDeath(:)   !the first snapshot when this subhalo does not exist
	integer*4,allocatable::SnapEnter(:)   !the snapshot right after which the subhalo becomes a satellite-sub,i.e.,SubRank>0
	integer*4,allocatable::ProHistID(:)   ! the mother history id for splitters, set to -1 for normal (non-splitter) subhalo's histories
	type (SubNodePre),allocatable::Node(:)  ! nodes of subhalo properties for all the histories
end type EVOLUTIONCAT_Pre

type SubNode
	integer*4 Mdm  !sub DM mass
	integer*4 Mhost !virial mass if possible,otherwise fof-mass
	real*4 Chost     !-1 if not NFW-fittable
	integer*4 SubID
	integer*4 HostID !host haloid,now dynamically extended in EvoCat
	integer*4 SubRank
end type SubNode   !subhalo mass after host-halo extension


type EVOLUTIONCAT
	real*4 PartMass   !particle mass, in units of 10^10Msun/h
	integer*4 NNode
	integer*4 NHist
	integer*4,allocatable::HistLen(:)    
	integer*4,allocatable::HistOffset(:)  
	integer*4,allocatable::SnapBirth(:)   
	integer*4,allocatable::SnapDeath(:)   
	integer*4,allocatable::SnapEnter(:)   
	integer*4,allocatable::ProHistID(:)  
	type (SubNode),allocatable::Node(:)
end type EVOLUTIONCAT

type ShardParam
	integer*4 SnapBirth  ! the starting snapshot of this shard
	integer*4 SnapTidal   ! the first snapshot when tidal stripping dominates over accretion
	integer*4 SnapRvir  ! the first snapshot when R<Rvir
	real*4 Mrate(2) !m/M at Rtidal, Rvir
	real*4 Mhost(2) ! host virial mass at Rtidal, Rvir
	real*4 Rhost(2)  !host virial radius
	real*4 Kappa(2) ! V^2/Vcirc^2
	real*4 j2(2)    ! circularity^2
	real*4 ErrK(2) !error estimation due to discrete output
	real*4 Errj2(2)
	real*4 Chost(2) !host concentration
	real*4 Csat(2) !satellite concentration
end type

type HistoryShards
	integer*4 NumHist
	integer*4 NumShards   
	integer*4,allocatable::NBirth(:) !number of shards for each history
	integer*4,allocatable:: HistOffset(:) ! starting shard position of each history in Par(:)
	type(ShardParam),allocatable::Par(:)  !parameters of each shard
end type

contains

subroutine load_sub2hist(Nsnap,sub2hist,Nsubs,subcatdir)
! load the table to convert subhaloid to history_id
integer*4:: Nsnap,i
integer*4,intent(out)::Nsubs
integer*4,allocatable,intent(out)::sub2hist(:)
character(*) subcatdir
character(1024) filename,snum

	write(snum,'(I3)') Nsnap
	do i=1,3 
	 if (snum(i:i) .eq. ' ') then 
		snum(i:i)='0'
	 end if
	end do
	filename=trim(subcatdir)//'/history/sub2hist_'//trim(snum)
	open(unit=20,file=filename,form='unformatted',status='old',action='read',access='stream')

	read(20) Nsubs
	allocate(sub2hist(0:Nsubs-1))
	read(20) sub2hist
	close(20)
end subroutine

subroutine load_evocat_pre(EvoCat,subcatdir)
!to load the primordial histories
type(EVOLUTIONCAT_Pre) EvoCat
character(*) subcatdir
character(1024) filename

filename=trim(subcatdir)//'/history/EvoCat_pre'
open(unit=20,file=filename,form='unformatted',status='old',action='read',access='stream')

read(20) EvoCat%NNode,EvoCat%NHist
allocate(EvoCat%HistLen(0:EvoCat%NHist-1))
allocate(EvoCat%HistOffset(0:EvoCat%NHist-1))
allocate(EvoCat%SnapBirth(0:EvoCat%NHist-1))
allocate(EvoCat%SnapDeath(0:EvoCat%NHist-1))
allocate(EvoCat%SnapEnter(0:EvoCat%NHist-1))
allocate(EvoCat%ProHistID(0:EvoCat%NHist-1))
allocate(EvoCat%Node(0:EvoCat%NNode-1))
read(20) EvoCat%HistLen
read(20) EvoCat%HistOffset
read(20) EvoCat%SnapBirth
read(20) EvoCat%SnapDeath
read(20) EvoCat%SnapEnter
read(20) EvoCat%ProHistID
read(20) EvoCat%Node

close(20)
end subroutine

subroutine free_evocat_pre(EvoCat)
type(EVOLUTIONCAT_PRE) EvoCat

deallocate(EvoCat%HistLen)
deallocate(EvoCat%HistOffset)
deallocate(EvoCat%SnapBirth)
deallocate(EvoCat%SnapDeath)
deallocate(EvoCat%SnapEnter)
deallocate(EvoCat%ProHistID)
deallocate(EvoCat%Node)
end subroutine

subroutine load_evocat(EvoCat,subcatdir)
!to load the host-halo extended histories
type(EVOLUTIONCAT) EvoCat
character(*) subcatdir
character(1024) filename

filename=trim(subcatdir)//'/history/EvoCat_rev'
open(unit=20,file=filename,form='unformatted',status='old',action='read',access='stream')

read(20) EvoCat%PartMass,EvoCat%NNode,EvoCat%NHist
allocate(EvoCat%HistLen(0:EvoCat%NHist-1))
allocate(EvoCat%HistOffset(0:EvoCat%NHist-1))
allocate(EvoCat%SnapBirth(0:EvoCat%NHist-1))
allocate(EvoCat%SnapDeath(0:EvoCat%NHist-1))
allocate(EvoCat%SnapEnter(0:EvoCat%NHist-1))
allocate(EvoCat%ProHistID(0:EvoCat%NHist-1))
allocate(EvoCat%Node(0:EvoCat%NNode-1))
read(20) EvoCat%HistLen
read(20) EvoCat%HistOffset
read(20) EvoCat%SnapBirth
read(20) EvoCat%SnapDeath
read(20) EvoCat%SnapEnter
read(20) EvoCat%ProHistID
read(20) EvoCat%Node

close(20)
end subroutine

subroutine free_evocat(EvoCat)
type(EVOLUTIONCAT) EvoCat

deallocate(EvoCat%HistLen)
deallocate(EvoCat%HistOffset)
deallocate(EvoCat%SnapBirth)
deallocate(EvoCat%SnapDeath)
deallocate(EvoCat%SnapEnter)
deallocate(EvoCat%ProHistID)
deallocate(EvoCat%Node)
end subroutine

recursive function GetNode(EvoCat,HistID,Nsnap)
! to access the node for history HistID at snapshot Nsnap
type(EVOLUTIONCAT) EvoCat
integer*4::HistID,Nsnap
type(SubNode) GetNode
integer*4 NodeID,SnapBirth,SnapDeath,ProHistID

SnapBirth=EvoCat%SnapBirth(HistID)
SnapDeath=EvoCat%SnapDeath(HistID)

if (Nsnap<SnapBirth) then
	ProHistID=EvoCat%ProHistID(HistID)
	if (ProHistID>0) then
	GetNode=GetNode(EvoCat,ProHistID,Nsnap)
	else
	GetNode=SubNode(0,0,0.0,-1,-2,-1)
	return
	end if
end if

if (Nsnap>=SnapDeath) then
	GetNode=SubNode(0,0,0.0,-1,-2,-1) 
	return
endif

NodeID=EvoCat%HistOffset(HistID)-SnapBirth+Nsnap
GetNode=EvoCat%Node(NodeID)
end function

subroutine load_historyshards(Shard,subcatdir)
!load shard parameters
type(HistoryShards) Shard
integer NumHist
character(*) subcatdir
character(1024) filename

filename=trim(subcatdir)//'/history/HistoryShards'
open(unit=20,file=filename,form='unformatted',status='old',action='read',access='stream')

read(20) NumHist
Shard%NumHist=NumHist
allocate(Shard%NBirth(0:NumHist-1))
allocate(Shard%HistOffset(0:NumHist-1))

read(20) Shard%NumShards
allocate(Shard%Par(0:Shard%NumShards-1))

read(20) Shard%NBirth
read(20) Shard%HistOffset
read(20) Shard%Par

close(20)
end subroutine	

end module history_io
