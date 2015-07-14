SRC=sub_io.f90 read_example.f90
FC=ifort

read_example: $(SRC)
	$(FC) -o read_example $(SRC)

HSRC=history_io.f90 read_history.f90	
read_history: $(HSRC)
	$(FC) -o read_history $(HSRC)	

read_CrossSection: BranchIO.f90 read_CrossSection.f90
	$(FC) $^ -o $@ #-fdefault-integer-8 -fdefault-real-8

clean:
	rm -f read_example read_history read_CrossSection
	rm -f *.mod *.o
