program read
implicit none
real :: start, finish
call cpu_time(start)
!$omp parallel
! !$omp single
! !$omp task
!$omp sections
!$omp section
call read_id()
! !$omp end task
! !$omp task
!$omp section
call read_pos()
! !$omp end task
! !$omp task
!$omp section
call read_vel()
!$omp end sections
! !$omp end task
! !$omp end single 
!$omp end parallel
call cpu_time(finish)
print '("total Time = ",E8.2," seconds.")',finish-start
end program

subroutine read_id()
INTEGER*8 ::np,ips,i
real*4::ztp,omgt,lbdt,boxsize,xscale,vscale
integer*8, allocatable:: id(:)
real :: start, finish
call cpu_time(start)
open(10,file='/mnt/ddnfs/jxhan/uv35/6610/simu/id6610.0247.01',form='unformatted',status='old')
! read(10) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale
np=3072
np=np*np*np
allocate(id(np/20))
read(10) id !(id(i),i=1,1000)
close(10)
call cpu_time(finish)
print '("id Time = ",E8.2," seconds.")',finish-start
return
end

subroutine read_pos()
INTEGER*8 ::np,ips,i
real*4::ztp,omgt,lbdt,boxsize,xscale,vscale
real*4, allocatable:: pos(:,:)
real :: start, finish
call cpu_time(start)
open(11,file='/mnt/ddnfs/jxhan/uv35/6610/simu/pos6610.0247.01',form='unformatted',status='old')
read(11) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale
allocate(pos(3,np/40))
read(11) pos ! (pos(1, i),i=1,1000)
close(11)
call cpu_time(finish)
print '("pos Time = ",E8.2," seconds.")',finish-start
return
end

subroutine read_vel()
INTEGER*8 ::np,ips,i
real*4::ztp,omgt,lbdt,boxsize,xscale,vscale
real*4, allocatable:: vel(:,:)
real :: start, finish
call cpu_time(start)
open(12,file='/mnt/ddnfs/jxhan/uv35/6610/simu/vel6610.0247.01',form='unformatted',status='old')
read(12) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale
allocate(vel(3,np/40))
read(12) vel !(vel(1, i),i=1,1000)
close(12)
call cpu_time(finish)
print '("vel Time = ",E8.2," seconds.")',finish-start
return
end