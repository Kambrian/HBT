FTNLIB=-L/opt/intel/fc/10.1.018/lib -lifport -lifcore 
#~ FTNLIB= -lgfortran
OBJS=testfortread.o fortread.o
CC=icc
FC=ifort
#~ CFLAGS=-g

ttfread: $(OBJS) Makefile
		$(CC)  -o ttfread $(OBJS) $(FTNLIB) 

testfortread.o:testfortread.c Makefile

fortread.o:fortread.f90 Makefile
		$(FC) fortread.f90 -c 
		
clean:
	rm $(OBJS) ttfread -rf
