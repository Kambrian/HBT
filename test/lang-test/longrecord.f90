program longrecord
integer,parameter::nid=1024*1024*1024
integer:: i,j
!~ integer,allocatable:: id(:,:)
real*4:: f(2)

f(1)=0.3
f(2)=0.4

!~ Allocate(ID(3,nid))
!~ forall(i=1:3,j=1:nid)
!~ id(i,j)=0
!~ end forall

open(10,file='longrecord.dat',status='new',form='unformatted')
write(10) ((i,i=0,nid-1),j=1,3)     !((ID(i,j),i=1,3),j=1,nid)
write(10) f
close(10)
end program longrecord
