program read_CrossSection

use BranchIO
implicit none

type(CrossSection) sec
type(BranchNode) node
character(len=1024):: subdir
integer Nsnap,branchID

Nsnap=10
branchID=20000
!subdir='/mnt/ddnfs/jxhan/uv35/6610/subcat/BranchTable' 
subdir='/mnt/ddnfs/jxhan/6610BranchTable'
call load_cross_section(Nsnap,sec,subdir)

write(*,*)  'Snapshot',Nsnap
write(*,*) 'Total Number of Branches=',sec%NumNode
write(*,*) 'printing BranchID=',branchID
write(*,*) '============================================'
node=sec%Node(branchID)
write(*,*) 'HostID=',node%HostID
write(*,*) 'SubID=', node%SubID
write(*,*) 'Rank=', node%SubRank
write(*,*) 'Most Bound ParticleID=', node%MstBndID
write(*,*) 'Bound Mass=', node%NpBnd, 'particles'
write(*,*) 'Peak Bound Mass=', node%NpBndPeak
write(*,*) 'Peak Snapshot=', node%SnapNumPeak
write(*,*) 'Vmax[km/s]=', node%Vmax
write(*,*) 'Peak Vmax=', node%VmaxPeak
write(*,*) 'Peak Vmax Snapshot=', node%SnapNumVpeak
write(*,*) 'Current Pos[kpc/h]=', node%MstBndPos
write(*,*) 'Vel[km/s]=', node%MstBndVel

call free_cross_section(sec)

end program read_CrossSection
