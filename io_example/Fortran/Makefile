SRC=sub_io.f90 read_example.f90
FC=ifort

read_example: $(SRC)
	$(FC) -o read_example $(SRC)

HSRC=history_io.f90 read_history.f90	
read_history: $(HSRC)
	$(FC) -o read_history $(HSRC)	

clean:
	rm -f read_example read_history
	rm -f *.mod *.o
