include Makefile.runs

CFLAGS=  -openmp -include $(PARAM) -I. -I $(IODIR)
LDFLAGS= -openmp  $(FTNLIB)

INCL	= intra_vars.h $(IODIR)/iovars.h proto.h Makefile $(PARAM)
OBJS_IO=sub_IO.o $(IODIR)/user_IO.o mymath.o intra_vars.o $(IODIR)/iovars.o 
OBJS_HIST=history/history_io.o

all: iolib histio

iolib:$(OBJS_IO) $(OBJS_FTN)
		ar -r lib$(RUN_NUM)io.a $(OBJS_IO) $(OBJS_FTN)
	
histio: $(OBJS_HIST) iolib
	ar -r libhistio.a $(OBJS_HIST)

$(OBJS_IO): $(INCL)

$(OBJS_FTN): $(IODIR)/fortread.f90
		$(FC) -c $(IODIR)/fortread.f90 -o $(IODIR)/fortread.o

$(OBJS_HIST):history/history_proto.h history/history_vars.h

clean:
	rm -f $(OBJS_FTN) lib*io.a *_io/*.o
	rm -f core.* *.o
