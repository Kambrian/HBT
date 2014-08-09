include ../Makefile.runs

BTDIR=..

CFLAGS=-I$(BTDIR) -I$(BTDIR)/$(IODIR) -I. -include $(PARAM) -openmp -g 
LDFLAGS= -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB) -limf -openmp -g

#this is only necessary when code+data size>2Gb
#CFLAGS+=-mcmodel=medium -share-intel

VPATH=$(BTDIR)

SOURCES = $(wildcard *.c) $(wildcard **/*.c) $(wildcard **/*/*.c)
EXEC =$(basename $(notdir $(SOURCES)))

#this rule automatically adds any additional dependences which are specified later to the recipe.
#so do not need to specify $(LDFLAGS) anymore
% : %.o
	$(CC) $^ $(LDFLAGS) -o $@

%.o : %.c
	$(CC) $< $(CFLAGS) -c -o $@
	      
#this is not so good if your .c file depend on another .c file but that file cannot be compiled alone (e.g, linkedlist.c)
% : %.c
	$(CC) $^ $(CFLAGS) $(LDFLAGS) -o $@
		     
Vmax:tree.o

clean:
	rm -f *.o $(EXEC)
