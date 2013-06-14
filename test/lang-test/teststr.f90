program main
implicit none
character(50) str
str="hello world!"
write(*,'(2A)') str,str
write(*,'(2A)') trim(str),str
write(*,'(A)') trim(str)//str,str//str
write(*,'(A,I3,I,A,I2)') trim(str),17,5,trim(str),17
write(*,'(I)') 3,10
end program main
