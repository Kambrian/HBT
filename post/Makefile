include ../Makefile.runs
BTDIR=..

USRINCL=-I /home/kambrain/opt/include
USRLIB=-L /home/kambrain/opt/lib
CFLAGS=-I$(BTDIR) -I$(BTDIR)/$(IODIR) -include $(PARAM) -openmp $(USRINCL) 
LDFLAGS= -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB) -openmp $(USRLIB)

sp2pro:sp2pro.c Makefile
	$(CC) -o sp2pro sp2pro.c $(CFLAGS) $(LDFLAGS) 

sp2pro.patch:sp2pro.patch.c Makefile
	$(CC) -o sp2pro.patch sp2pro.patch.c $(CFLAGS) $(LDFLAGS) 
		
	
pro2dest:pro2dest.o  Makefile
	$(CC) -o pro2dest pro2dest.o $(CFLAGS) $(LDFLAGS) 
				
traceup:traceup.c Makefile
	$(CC) -o traceup traceup.c $(CFLAGS) $(LDFLAGS) 

DeathMassfun:DeathMassfun.c Makefile
	$(CC) -o DeathMassfun DeathMassfun.c $(CFLAGS) $(LDFLAGS)	
				
tidal_tree: tidal_tree.o treesearch_binlin.o tree.o Makefile
		$(CC) -o tidal_tree tidal_tree.o treesearch_binlin.o tree.o  $(LDFLAGS) 

treesearch_binlin.o: treesearch.c Makefile
		$(CC) -c -o treesearch_binlin.o treesearch.c $(CFLAGS) -DBIN_LIN_0

tree.o: $(BTDIR)/tree.c	Makefile
				$(CC) -c $(BTDIR)/tree.c $(CFLAGS) 
				
hierarchy.o: $(BTDIR)/hierarchy.c	Makefile
				$(CC) -c $(BTDIR)/hierarchy.c $(CFLAGS) 
				
clean:
	rm -f *.o tidal_tree  sp2pro traceup pro2dest sp2pro.patch DeathMassfun
