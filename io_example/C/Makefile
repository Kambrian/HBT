include library/Makefile.runs

#Sample Programs==
#loadview.c: loads all the subhalo data for a given snapshot
#load_virial.c: loads halo virial properties
#load_history.c: loads subhalo merger histories
#
# before you can compile these programs, you need to compile the libraries first, by:
#    cd library
#    make
#
# once you've made the libraries, come back to this directory and simply call "make" will make all the executables
#    make
# or specify the executable name to only make that program, e.g,
#    make load_virial

#simply add your own file to this SRC list,as long as they are single file
SRC=loadview.c load_virial.c load_history.c
EXEC=$(basename $(SRC))
OBJS = $(SRC:%.c=%.o)

CFLAGS=-include $(PARAM) -I library -I library/$(IODIR)  -openmp -g
LDFLAGS=-openmp $(FTNLIB) -L library -l$(RUN_NUM)io 

all:$(EXEC)

# a static pattern rule, to compile any member x of $(EXEC) using
#  $(CC) x.o $(LDFLAGS) -o x
$(EXEC):%:%.o
	$(CC) $^ $(LDFLAGS) -o $@

load_history: CFLAGS+=-I library/history
		LDFLAGS+=-lhistio
# the above rule is equivalent to listing them implicitly below
#loadview: loadview.o

#load_virial: load_virial.o

#load_history: load_history.o
	
clean:
	rm -f $(EXEC) $(OBJS) core.* core	
