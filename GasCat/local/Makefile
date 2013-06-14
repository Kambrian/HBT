BTDIR=../bt
PARAM=param6702.h

CC			=icc
CFLAGS=-I$(BTDIR) -include $(PARAM) -openmp
#~ LDFLAGS= -L$(BTDIR) -lbtio 
gaslib: gas_stuff.o tree.o
		ar rc libgas.a gas_stuff.o tree.o
		
gas_stuff.o: gas_stuff.c gas_stuff.h  Makefile
				$(CC) -c gas_stuff.c $(CFLAGS) -g

tree.o: $(BTDIR)/tree.c	Makefile
				$(CC) -c $(BTDIR)/tree.c $(CFLAGS) -g
clean:
	rm -rf *.o gaslib 
