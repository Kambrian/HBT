program read_history
use history_io
implicit none
type (evolutioncat) evocat
type (historyshards) shard
integer*4,allocatable::sub2hist(:)
character(len=1024):: SubDir
integer*4 Nsnap,Nsubs,subid,histid
type(subnode) node

Nsnap=15
subid=11
SubDir='/home/kambrain/data/6113/subcat' 
!subhalo catalogues for different simulations all located under /home/kambrain/data/***/subcat
! currently the evolution catalogue has only been produced for 6113

call load_evocat(evocat,subdir)
call load_historyshards(shard,subdir)
call load_sub2hist(Nsnap,sub2hist,Nsubs,subdir)
histid=sub2hist(subid)
write(*,*)  'Snapshot',Nsnap
write(*,*) 'Number of histories=',evocat%NHist,'; total number of subhalos=',evocat%NNode
write(*,*) 'subid=',subid,'corresponds to histid=',histid
write(*,*) 'it lives',evocat%histlen(subid),'snapshots,birth at Snap=',evocat%snapbirth(subid),', death at Snap=',evocat%snapdeath(subid)
node=getnode(evocat,histid,nsnap)
write(*,*)  'its mass=',node%Mdm,', host virial mass=',node%Mhost, ', host concentration=',node%chost
write(*,*)  'its subid=',node%subid, ', hostid=',node%hostid
write(*,*) 'this history has', shard%nbirth(histid),'merger(s)'
write(*,*) 'the first merger crosses virial radius at snapshot',shard%par(shard%histoffset(histid))%SnapRvir
write(*,*) 'with circularity',sqrt(shard%par(shard%histoffset(histid))%j2(2))

call free_evocat(evocat)
deallocate(sub2hist)
end program read_history
