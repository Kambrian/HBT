function alloc_file_unit()
integer alloc_file_unit
integer,save::fileno=10
fileno=fileno+1
if(fileno.eq.100) then !skip 100-102, reserved for stdin,out,err
  fileno=103
endif
alloc_file_unit=fileno  
return
end function alloc_file_unit

subroutine open_fortran_file(filename,fileno,bigendian_flag,error_stat)
implicit none
integer alloc_file_unit
character(1024) filename
integer*4:: fileno,bigendian_flag,error_stat
!$omp critical (fortran_open) !this is necessary for parallel fortran read. no effect is outside parallel region.
fileno=alloc_file_unit()
if(bigendian_flag.NE.0) then
open(unit=fileno,FILE=filename,form='unformatted',status='old',iostat=error_stat,convert='big_endian')
else
open(unit=fileno,FILE=filename,form='unformatted',status='old',iostat=error_stat)
endif
!$omp end critical (fortran_open)
end subroutine

subroutine close_fortran_file(fileno)
implicit none
integer*4 fileno
!$omp critical (fortran_open)
close(fileno)
!$omp end critical (fortran_open)
end subroutine 

subroutine read_fortran_record1(arr,arr_len,fileno)
implicit none
integer*8 arr_len
integer*4 fileno,base_size
integer*1 arr(arr_len)
read(fileno) arr
end subroutine 
subroutine read_fortran_record2(arr,arr_len,fileno)
implicit none
integer*8 arr_len
integer*4 fileno,base_size
integer*2 arr(arr_len)
read(fileno) arr
end subroutine 
subroutine read_fortran_record4(arr,arr_len,fileno)
implicit none
integer*8 arr_len
integer*4 fileno,base_size
integer*4 arr(arr_len)
read(fileno) arr
end subroutine 
subroutine read_fortran_record8(arr,arr_len,fileno)
implicit none
integer*8 arr_len
integer*4 fileno,base_size
integer*8 arr(arr_len)
read(fileno) arr
end subroutine 
