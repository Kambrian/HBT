include ../Makefile.runs
BTDIR=..
GASDIR=../GasCat

EXEC	 =evocat_pre evocat_rev MBDCat\
		  trace_bin_rev unevo_mass pick_history 
MAIN_SRC =$(addsuffix .c, $(basename $(EXEC)))
OBJS_MAIN=$(MAIN_SRC:%.c=%.o)
INCL	 = history_vars.h history_proto.h
SRC		 = history_io.c
OBJS     = history_io.o
EXE_ALL  :=$(basename $(wildcard *.c))


CFLAGS=  -openmp -include $(PARAM) -I $(BTDIR) -I $(GASDIR) -I $(BTDIR)/$(IODIR) -I.
#GASLIB=-L$(GASDIR) -lgasio
LDFLAGS= -openmp -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB) $(GASLIB) -limf

% : %.c   #cancel this implicit rule first; otherwise it will override the following rule

% : %.c $(OBJS)      #this rule does not override the above implicit rule because of extra pre-requisites
	$(CC) -o $@ $^ $(CFLAGS) $(OBJS) $(LDFLAGS) 

all: $(EXEC)
	
#$(EXEC): $(OBJS)
#	$(CC) -o $@ $(addsuffix .c, $(basename $@)) $^ $(CFLAGS) $(LDFLAGS)
	
hist_diff: CFLAGS+= -g

histio: libhistio.a
		
libhistio.a:history_io.o
	ar -r libhistio.a history_io.o
	
depend:
	makedepend --$(CFLAGS)-- -Y $(SRC) $(MAIN_SRC)	


clean:
	\rm -f *.o $(EXE_ALL)
	\rm -f core.*
	\rm -f libhistio.a
