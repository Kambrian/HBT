subroutine open_fortran_file(filename,fileno,bigendian_flag,error_stat)
implicit none
character (len=1024) filename
integer*4 fileno,bigendian_flag,error_stat
if(bigendian_flag.NE.0) then
open(unit=fileno,FILE=filename,form='unformatted',status='old',iostat=error_stat,convert='big_endian')
else
open(unit=fileno,FILE=filename,form='unformatted',status='old',iostat=error_stat)
endif
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
subroutine close_fortran_file(fileno)
implicit none
integer*4 fileno
close(fileno)
end subroutine 
