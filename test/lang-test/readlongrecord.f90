program readlongrecord
integer,parameter::nid=1024*1024*1024
integer:: i,j
integer,allocatable::id(:,:)
real*4:: f(2)

allocate(ID(3,nid))
ID=1
open(10,file='longrecord.dat',status='old',form='unformatted')
read(10) ID     !((ID(i,j),i=1,3),j=1,nid)
read(10) f
close(10)
do i=1,3
do j=1,nid
if(ID(i,j).NE.0) then
	write(*,*) i,j,ID(i,j)
endif
enddo
enddo
write(*,*) f

end program readlongrecord