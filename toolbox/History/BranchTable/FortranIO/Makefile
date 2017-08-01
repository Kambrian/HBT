FC=ifort
FFLAGS=-i8 #-r8

#FC=gfortran
#FFLAGS=-fdefault-integer-8 #-fdefault-real-8

read_CrossSection: BranchIO.f90 read_CrossSection.f90
	$(FC) $^ -o $@ $(FFLAGS)

clean:
	rm -f read_CrossSection
	rm -f *.mod *.o
