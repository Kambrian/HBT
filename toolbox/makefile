include ../Makefile.runs

BTDIR=..
GASDIR=../GasCat
HISTDIR=../history

USRINCL=-I /home/kambrain/opt/include
USRLIB=-L /home/kambrain/opt/lib
CFLAGS=-I$(BTDIR) -I$(BTDIR)/$(IODIR) -I$(GASDIR) -I$(HISTDIR) -I. -include $(PARAM) $(USRINCL) -openmp -g -I/usr/include/mpi
LDFLAGS= -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB) $(USRLIB) -limf -openmp -g
GASLIB= -L$(GASDIR) -lgasio
HISTLIB= -L$(HISTDIR) -lhistio
# GSLINC=-I/home/qypeng/bin/gsl-1.15/include
# GSLLIB=-L/home/qypeng/bin/gsl-1.15/lib -lgsl -lgslcblas
GSLINC=-I/home/qypeng/softwares/gsl-1.15/include
GSLLIB=-L/home/qypeng/softwares/gsl-1.15/lib -lgsl -lgslcblas
#GSLINC=-I/gpfs/COSMA/gsl/intel_11.1/1.14/include
#GSLLIB=-L/gpfs/COSMA/gsl/intel_11.1/1.14/lib -lgsl -lgslcblas

#this is only necessary when code+data size>2Gb
#CFLAGS+=-mcmodel=medium -share-intel

VPATH=Check History Image MassFunction Profile \
      Profile/concentration  Profile/geometrical  Profile/sub  Profile/xray \
      $(BTDIR) $(BTDIR)/history

SOURCES = $(wildcard *.c) $(wildcard **/*.c) $(wildcard **/*/*.c)
EXEC =$(basename $(notdir $(SOURCES)))

echo:
	#$(EXEC)
	
#this rule automatically adds any additional dependences which are specified later to the recipe.
#so do not need to specify $(LDFLAGS) anymore
% : %.o
	$(CC) $^ $(LDFLAGS) -o $@

%.o : %.c
	$(CC) $< $(CFLAGS) -c -o $@
	      
#this is not so good if your .c file depend on another .c file but that file cannot be compiled alone (e.g, linkedlist.c)
% : %.c
	$(CC) $^ $(CFLAGS) $(LDFLAGS) -o $@
		     
#history	
NEED_GAS+=pick_subsnap trace_bin_ext trace_bin_rev trace_bin
NEED_HIST+=pick_subsnap trace_bin_ext trace_bin_rev trace_bin quasiMmax spMmax wrongSpMain bigSplinters
NEED_TREE+=check_centers
p: pick_Chunyan
	mv pick_Chunyan pick_Chunyan4.${RUN_NUM}
	
#massfun
NEED_HIST+=pick_halosample
NEED_SUBFIND+=massfun_cosm
	
massfun_plotN massfun_pickN abundanceN: %N: %.c
	$(CC) -o $@ $^ -DNORM $(CFLAGS) $(LDFLAGS)  								
			 								
#image
#NEED_SUBFIND+=subparticles subcen snaps snap2 snap3 snaps_phase snaps_phase_diff
NEED_TREE+=snap2 snaps_phase snaps_reacc snaps_energy snaps_energy2
NEED_TREESEARCH+=snap2
NEED_GSL+=subparticles NFW_fit NFW_halo
		
snap_lowres:CFLAGS=-I$(BTDIR) -include $(PARAM_LOW)		

#profile
#haloprof.o: subfindread.o
gashaloprof dynamprof dynamprof_sph:LDFLAGS+=$(GASLIB)	
veldistr: tree.o 
halocenter: treesearch.o tree.o 
haloprof_nosub: CFLAGS+=-DWITHOUT_SUBCAT
haloprof_nosub: haloprof
	mv haloprof haloprof_nosub
haloprof.$(RUN_NUM): haloprof
	mv haloprof haloprof.$(RUN_NUM)

#Profile/concentration/
NFW_fit.o NFW_halo.o: NFW_fun.c lnNFW_fun.c
#NFW_fit NFW_halo: CFLAGS+=$(GSLINC) LDFLAGS+=$(GSLLIB)

#Profile/geometrical							
tidalradius tidalradius_log tidalradius_aqua: subfindread.o
fofpot:  tree.o  
fofpot: LDFLAGS+=tree.o
subshape: subfindread.o CFLAGS+=$(GSLINC) LDFLAGS+=$(GSLLIB)				
notts_shape: CFLAGS+=$(GSLINC) LDFLAGS+=$(GSLLIB)
 				
#Profile/sub				
subprofile: CFLAGS+=-USUBFIND_DIR
subprofileF: subprofile.c subfindread.o
	$(CC) -o $@ $< subfindread.o $(CFLAGS) $(LDFLAGS) 				
Vmax:tree.o 
Vmax:CFLAGS+=-USUBFIND_DIR		

plotSubInSub:CFLAGS+=-UDYNAMICAL
plotSubInSub.dyn:plotSubInSub.c
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS) -DDYNAMICAL								
										
#Profile/xray/				
Lxsub gasprofile gasprofile2: LDFLAGS+=-lgsl -lgslcblas
				
#Check/				
			
loadview:CFLAGS+=-g
match_subcat back_binding:subfindread.o			

notts: saveNotts Vmax
				./Vmax
				./saveNotts

subfindread: subfindread.c Makefile
	$(CC) -o subfindread subfindread.c $(CFLAGS) $(LDFLAGS) -DSTANDALONE

$(NEED_GAS): LDFLAGS+=$(GASLIB) 
$(NEED_HIST):LDFLAGS+=$(HISTLIB)
$(NEED_GSL): CFLAGS+=$(GSLINC) 
$(NEED_GSL): LDFLAGS+=$(GSLLIB)
   
$(NEED_SUBFIND): subfindread.o
#$(NEED_SUBFIND): LDFLAGS+=subfindread.o

$(NEED_TREE): tree.o 
#$(NEED_TREE): LDFLAGS+=tree.o

$(NEED_TREESEARCH): treesearch.o 
#$(NEED_TREESEARCH): LDFLAGS+=tree.o

hdf_util.o: 

StarPot: hdf_util.o tree.o
StarPot : LDFLAGS+=-lhdf5_hl -lhdf5

saveHDFhalo saveHDFhalo_addstrm saveHDFhalo_debug: hdf_util.o
saveHDFhalo saveHDFhalo_addstrm saveHDFhalo_debug: LDFLAGS+=-lhdf5_hl -lhdf5
saveHDFhalo saveHDFhalo_debug: subfindread.o

collect_mostbounds: hdf_util.o history_io.o
collect_mostbounds: LDFLAGS+=-lhdf5_hl -lhdf5 

collect_mostbounds_virial: hdf_util.o history_io.o
collect_mostbounds_virial: LDFLAGS+=-lhdf5_hl -lhdf5 

clean:
	rm -f *.o $(EXEC)
