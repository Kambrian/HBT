include ../Makefile.runs
BTDIR=..
GASDIR=../GasCat

EXEC	=	history.$(RUN_NUM)
INCL	= history_vars.h history_proto.h
OBJS  = history.o history_io.o
CFLAGS=  -openmp -include $(PARAM) -I $(BTDIR) -I $(GASDIR) -I $(BTDIR)/$(IODIR) -I.
LDFLAGS= -openmp -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB)
GASLIB=-L$(GASDIR) -lgasio

$(EXEC): $(OBJS) Makefile
	$(CC) $(OBJS) $(LDFLAGS) $(GASLIB) -o  $(EXEC) 
history_extend: history_extend.c histio
	$(CC) history_extend.c $(CFLAGS) -L. -lhistio $(LDFLAGS)  -o history_extend
history_rev: history_revised.c histio
	$(CC) history_revised.c $(CFLAGS) -L. -lhistio $(LDFLAGS)  -o history_rev

history_brf: history_brief.c 
	$(CC) history_brief.c $(CFLAGS) $(LDFLAGS)  -o history_brf

evo-cat: evo-cat.c 
	$(CC) evo-cat.c $(CFLAGS) $(LDFLAGS)  -o evo-cat

evo-cat-rev: evo-cat-rev.c 
	$(CC) evo-cat-rev.c $(CFLAGS) $(LDFLAGS)  -o evo-cat-rev	

hist_diff: hist_diff.c 
	$(CC) hist_diff.c $(CFLAGS) $(LDFLAGS)  -o hist_diff	-g		
hist_dbg: hist_dbg.c histio
	$(CC) hist_dbg.c $(CFLAGS) -L. -lhistio $(LDFLAGS)   -o hist_dbg	-g
histio:history_io.o
	ar -r libhistio.a history_io.o
	
clean:
	rm -f *.o $(EXEC) history_extend hist_dbg history_rev history_brf
	rm -f core.* 
	rm -f libhistio.a
