subroutine write_fortran_file(filename,fileno)
implicit none
character(1024) filename
integer*4:: fileno,error_stat
!$omp critical (fortran_open)
open(unit=fileno,FILE=filename,form='unformatted',status='new',iostat=error_stat)
!$omp end critical (fortran_open)
write(fileno) fileno
close(fileno)
end subroutine

subroutine read_fortran_file(filename,fileno, i)
implicit none
character(1024) filename
integer*4:: fileno,error_stat,i
!$omp critical (fortran_open)
write(*,*) fileno,filename
open(unit=fileno,file=filename,form='unformatted',status='old',iostat=error_stat)
!$omp end critical (fortran_open)
read(fileno) i
! !$omp critical
! write(*,*) i,fileno,filename
close(fileno)
end subroutine