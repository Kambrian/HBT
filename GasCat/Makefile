include ../Makefile.runs
BTDIR=..

EXEC	=	gashalo.$(RUN_NUM)
INCL	= gas_vars.h gas_proto.h
OBJS_GAS=gas_io.o gas_binding.o gas_vars.o
OBJS_TREE=treesearch.o tree.o
CFLAGS=  -openmp -include $(PARAM) -I $(BTDIR) -I $(BTDIR)/$(IODIR) -I.
LDFLAGS= -openmp -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io

$(EXEC): gashalo.o $(OBJS_GAS) $(OBJS_TREE) Makefile
	$(CC) gashalo.o $(OBJS_GAS) $(OBJS_TREE) $(LDFLAGS) -o  $(EXEC) 

gashalo.o: $(INCL) Makefile
$(OBJS_GAS):$(INCL) Makefile

GT: GT.o $(OBJS_GAS) $(OBJS_TREE)
	$(CC) GT.o $(OBJS_GAS) $(OBJS_TREE) $(LDFLAGS) -o GT.$(RUN_NUM)

gas_check: gas_check.o gas_io.o gas_vars.o
	$(CC) gas_check.o gas_io.o gas_vars.o $(LDFLAGS) -o gas_check	

gas_check.o:$(INCL)
	$(CC) gas_check.c -c -g $(CFLAGS)
	 
follow_gashalo:follow_gashalo.o tree.o $(OBJS_GAS)
	$(CC) follow_gashalo.o tree.o $(OBJS_GAS) $(LDFLAGS) -o follow_gashalo
	

treesearch.o:$(BTDIR)/treesearch.c	Makefile
	$(CC) -c $(BTDIR)/treesearch.c $(CFLAGS)	
	
tree.o: $(BTDIR)/tree.c	Makefile
	$(CC) -c $(BTDIR)/tree.c $(CFLAGS)	

OBJS_IO=gas_io.o gas_vars.o
gasio: $(OBJS_IO)	
		ar -r libgasio.a $(OBJS_IO) $(OBJS_FTN)
clean:
	rm -f *.o $(EXEC) follow_gashalo GT.$(RUN_NUM) gas_check
	rm -f core.*
	rm -f libgasio.a
