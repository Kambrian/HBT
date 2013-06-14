include Makefile.runs

OBJS_MAIN= HBT.o FoF.o FoF_mpi.o hiercheck.o 
OBJS_COMM= intra_vars.o $(IODIR)/iovars.o mymath.o sub_IO.o $(IODIR)/user_IO.o hierarchy.o treesearch.o tree.o binding_minpot.o $(OBJS_FTN)
INCL	 = datatypes.h intra_vars.h $(IODIR)/iovars.h proto.h Makefile $(PARAM)
EXEC     = $(OBJS_MAIN:%.o=%.$(RUN_NUM))
targets  = $(basename $(OBJS_MAIN))

ifeq ($(CC), icc)
CFLAGS=  -openmp -include $(PARAM) -I. -I $(IODIR) -g  
LDFLAGS= -openmp  $(FTNLIB) -limf
else 												#using gcc
CFLAGS =-fopenmp -include $(PARAM) -I. -I $(IODIR) -g
LDFLAGS=-lm -lgfortran -fopenmp 
endif

#this is used only when code+data size >2GB
#CFLAGS+=-mcmodel=medium -shared-intel


default: HBT
all:$(targets)
#the implicit rule search is too picky; let's define it static.
$(targets):% : %.$(RUN_NUM) ;
%.$(RUN_NUM): %.o $(OBJS_COMM)
	$(CC) $(LDFLAGS) $^ -o $@
	
#FoF_mpi.$(RUN_NUM): LDFLAGS+=-lmpi
FoF_mpi.$(RUN_NUM): CC=mpicc
$(OBJS_MAIN) $(OBJS_COMM): $(INCL)

$(OBJS_FTN): $(IODIR)/fortread.f90
		$(FC) -c $(IODIR)/fortread.f90 -o $(IODIR)/fortread.o

iolib:$(IODIR)/lib$(RUN_NUM)io.a ;
$(IODIR)/lib$(RUN_NUM)io.a: sub_IO.o $(IODIR)/user_IO.o mymath.o intra_vars.o $(IODIR)/iovars.o $(OBJS_FTN)
	ar -r $(IODIR)/lib$(RUN_NUM)io.a $^

sync4700:
	rsync -avz $(shell pwd) kambrain@a4700:data/

syncuv:
	rsync -e "ssh -p 4702" -avz $(shell pwd) jxhan@localhost:data/HBT/

clean:
	rm -f $(OBJS_COMM) $(OBJS_MAIN) $(HCHECK) $(IODIR)/lib*io.a $(EXEC)
	rm -f core.* *.o
	
.PHONY: $(targets) iolib clean

