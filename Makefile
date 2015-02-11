include Makefile.runs

OBJS_MAIN= $(addsuffix .o, $(targets))
OBJS_COMM= intra_vars.o $(IODIR)/iovars.o mymath.o sub_IO.o $(IODIR)/user_IO.o hierarchy.o treesearch.o tree.o binding_minpot.o $(OBJS_FTN)
INCL	 = datatypes.h intra_vars.h $(IODIR)/iovars.h proto.h Makefile $(PARAM) 
targets  = HBT FoF FoF_mpi #hiercheck.o 
EXEC     = $(addsuffix .$(RUN_NUM), $(targets))

ANAL     = haloprof snap2 massfun_plot NFW_fit snaps #analysis routines
targets += $(ANAL)

SRC =$(OBJS_COMM:%.o=%.c) $(addsuffix .c, $(targets))

VPATH=anal/Image anal/MassFunction anal/Profile \
      anal/Profile/concentration   
      
ifeq ($(CC), icc)
CFLAGS=  -openmp -include $(PARAM) -I. -I $(IODIR) -g  
LDFLAGS= -openmp  $(FTNLIB) -limf
else 												#using gcc
CFLAGS =-fopenmp -include $(PARAM) -I. -I $(IODIR) -g
LDFLAGS=-lm -lgfortran -fopenmp 
endif

#this is used only when code+data size >2GB
#CFLAGS+=-mcmodel=medium -shared-intel
NFW_fit: CFLAGS+=$(GSLINC)
NFW_fit: LDFLAGS+=$(GSLLIB)
snap2:   CFLAGS+=-DDISABLE_HALO_PARA

default: HBT
all:$(targets) ReBuild
daily:HBT $(ANAL) ReBuild
ReBuild:
	@$(MAKE) snap2 -B
	
#the implicit rule search is too picky; let's define it static.
$(targets):% : %.$(RUN_NUM) ;
%.$(RUN_NUM): %.o $(OBJS_COMM)
	$(CC) $(LDFLAGS) $^ -o $@
	
#FoF_mpi.$(RUN_NUM): LDFLAGS+=-lmpi
FoF_mpi.$(RUN_NUM): CC=mpicc
$(OBJS_MAIN) $(OBJS_COMM): $(INCL) Makefile

$(OBJS_FTN): $(IODIR)/fortread.f90
		$(FC) -c $(IODIR)/fortread.f90 -o $(IODIR)/fortread.o

iolib:$(IODIR)/lib$(RUN_NUM)io.a ;
$(IODIR)/lib$(RUN_NUM)io.a: sub_IO.o $(IODIR)/user_IO.o mymath.o intra_vars.o $(IODIR)/iovars.o $(OBJS_FTN)
	ar -r $(IODIR)/lib$(RUN_NUM)io.a $^
	
synccosma: clean
# 	rsync -avzL $(shell pwd) -e "ssh -p 4800" jvbq85@localhost:data/HBT/code/
	rsync -avzL $(shell pwd) jvbq85@cosma-c:data/HBT/code/
	
sync4700:
	rsync -avz $(shell pwd) kambrain@a4700:data/

syncuv:
	rsync -e "ssh -p 4702" -avz $(shell pwd) jxhan@localhost:data/HBT/

syncuvj:
	rsync -e "ssh -p 4719" -avz $(shell pwd) jxhan@localhost:HBT/

sync: clean
	rsync -avzL $(shell pwd) sussing@suss:Working/HBT/
syncback: clean
	rsync -avz sussing@suss:Working/HBT/v8.7c/ $(shell pwd)
	
echo:
	@echo $(EXEC)
clean:
	rm -f $(OBJS_COMM) $(OBJS_MAIN) $(HCHECK) $(IODIR)/lib*io.a 
	rm -f core.* *.o
distclean: clean
	rm -f $(EXEC)
depend:
	makedepend --$(CFLAGS)-- -Y $(SRC_COMM) $(SRC_MAIN) $(SRC)
	
	
.PHONY: $(targets) iolib clean distclean echo depend ReBuild

# DO NOT DELETE

intra_vars.o: param/paramHY300.h datatypes.h intra_vars.h
gadget_io/iovars.o: param/paramHY300.h datatypes.h intra_vars.h
gadget_io/iovars.o: gadget_io/iovars.h
mymath.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h
mymath.o: proto.h
sub_IO.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h
sub_IO.o: proto.h
gadget_io/user_IO.o: param/paramHY300.h datatypes.h intra_vars.h
gadget_io/user_IO.o: gadget_io/iovars.h proto.h
hierarchy.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h
hierarchy.o: proto.h
treesearch.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h
treesearch.o: proto.h
tree.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h
tree.o: proto.h
binding_minpot.o: param/paramHY300.h datatypes.h intra_vars.h
binding_minpot.o: gadget_io/iovars.h proto.h
HBT.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h proto.h
FoF.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h proto.h
FoF_mpi.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h
FoF_mpi.o: proto.h
HBT.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h proto.h
FoF.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h proto.h
FoF_mpi.o: param/paramHY300.h datatypes.h intra_vars.h gadget_io/iovars.h
FoF_mpi.o: proto.h
haloprof.o: param/paramHY300.h
snap2.o: param/paramHY300.h
massfun_plot.o: param/paramHY300.h
NFW_fit.o: param/paramHY300.h
