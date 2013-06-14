include ../Makefile.runs

BTDIR=..
GASDIR=../GasCat
HISTDIR=../history

USRINCL=-I /home/kambrain/opt/include
USRLIB=-L /home/kambrain/opt/lib
CFLAGS= -I$(BTDIR) -I$(BTDIR)/$(IODIR) -I$(GASDIR) -I$(HISTDIR) -I. -include $(PARAM) $(USRINCL) -openmp -g #-Wbrief
LDFLAGS= -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB) $(USRLIB) -openmp -g
GASLIB= -L$(GASDIR) -lgasio
HISTLIB= -L$(HISTDIR) -lhistio
GSLINC=-I/gpfs/COSMA/gsl/intel_11.1/1.14/include
GSLLIB=-L/gpfs/COSMA/gsl/intel_11.1/1.14/lib -lgsl -lgslcblas

#history
pick_MADorbit: History/pick_MADorbit.c Makefile
	$(CC) -o pick_MADorbit History/pick_MADorbit.c $(CFLAGS) $(LDFLAGS) 
orbitparam: History/orbitparam.c tree.o Makefile
	$(CC) -o orbitparam History/orbitparam.c tree.o $(CFLAGS) $(LDFLAGS) 
pick_subsnap: History/pick_subsnap.c Makefile
	$(CC) -o pick_subsnap History/pick_subsnap.c $(CFLAGS) $(LDFLAGS) $(GASLIB) $(HISTLIB)
trace_bin_ext:History/trace_bin_ext.c Makefile
	$(CC) -o trace_bin_ext History/trace_bin_ext.c $(CFLAGS) $(LDFLAGS) $(GASLIB) $(HISTLIB)	
trace_bin_rev:History/trace_bin_rev.c Makefile
	$(CC) -o trace_bin_rev History/trace_bin_rev.c $(CFLAGS) $(LDFLAGS) $(GASLIB) $(HISTLIB)		
trace_bin:History/trace_bin.c Makefile
	$(CC) -o trace_bin History/trace_bin.c $(CFLAGS) $(LDFLAGS) $(GASLIB) $(HISTLIB)
steller: History/steller.c Makefile
	$(CC) -o steller History/steller.c $(CFLAGS) $(LDFLAGS)
steller2: History/steller2.c Makefile
	$(CC) -o steller2 History/steller2.c $(CFLAGS) $(LDFLAGS)
Decay: History/Decay.c Makefile
	$(CC) -o Decay History/Decay.c $(CFLAGS) $(LDFLAGS)	
check_centers: History/check_centers.c tree.o Makefile
	$(CC) -o check_centers History/check_centers.c tree.o $(CFLAGS) $(LDFLAGS)
		
stellerDeath: History/stellerDeath.c Makefile
	$(CC) -o stellerDeath History/stellerDeath.c $(CFLAGS) $(LDFLAGS)	
lostsub: History/lostsub.c Makefile
	$(CC) -o lostsub History/lostsub.c	$(CFLAGS) $(LDFLAGS)			
HostHist: History/HostHist.c Makefile
	$(CC) -o HostHist History/HostHist.c	$(CFLAGS) $(LDFLAGS)		
SubHist: History/SubHist.c Makefile
	$(CC) -o SubHist History/SubHist.c	$(CFLAGS) $(LDFLAGS)		
quasiMmax:History/quasiMmax.c Makefile
	$(CC) -o quasiMmax History/quasiMmax.c $(CFLAGS) $(HISTLIB) $(LDFLAGS) 
spMmax:History/spMmax.c Makefile
	$(CC) -o spMmax History/spMmax.c $(CFLAGS) $(HISTLIB) $(LDFLAGS) 
wrongSpMain:History/spMmax.c Makefile
	$(CC) -o wrongSpMain History/wrongSpMain.c $(CFLAGS) $(HISTLIB) $(LDFLAGS) 	
bigSplinter:History/bigSplinters.c Makefile
	$(CC) -o bigSplinter History/bigSplinters.c $(CFLAGS) $(HISTLIB) $(LDFLAGS) 		
#MassFunction		
biggestsub:MassFunction/biggestsub.c Makefile
				$(CC) -o biggestsub MassFunction/biggestsub.c $(CFLAGS) $(LDFLAGS)  
					
massfun:MassFunction/massfun.c Makefile
				$(CC) -o massfun MassFunction/massfun.c $(CFLAGS) $(LDFLAGS)  

massfun_all:MassFunction/massfun_all.c Makefile
				$(CC) -o massfun_all MassFunction/massfun_all.c $(CFLAGS) $(LDFLAGS)  				

massfun_cosm:MassFunction/massfun_cosm.c  subfindread.o Makefile
				$(CC) -o massfun_cosm MassFunction/massfun_cosm.c subfindread.o $(CFLAGS) $(LDFLAGS) # $(GASLIB)
					
massfun_rem:MassFunction/massfun_rem.c Makefile
				$(CC) -o massfun_rem MassFunction/massfun_rem.c $(CFLAGS) $(LDFLAGS) 					

massfun_plot:MassFunction/massfun_plot.c Makefile
				$(CC) -o massfun_plot MassFunction/massfun_plot.c $(CFLAGS) $(LDFLAGS)  				
massfunN_plot:MassFunction/massfun_plot.c Makefile
				$(CC) -o massfunN_plot MassFunction/massfun_plot.c -DNORM $(CFLAGS) $(LDFLAGS)  								

pick_halosample:MassFunction/pick_halosample.c Makefile
				$(CC) -o pick_halosample MassFunction/pick_halosample.c $(CFLAGS) $(HISTLIB) $(LDFLAGS)  
massfun_pick:MassFunction/massfun_pick.c Makefile
				$(CC) -o massfun_pick MassFunction/massfun_pick.c $(CFLAGS) $(LDFLAGS)  				
massfunN_pick:MassFunction/massfun_pick.c Makefile
				$(CC) -o massfunN_pick MassFunction/massfun_pick.c -DNORM $(CFLAGS) $(LDFLAGS)  
abundanceN:MassFunction/abundance.c Makefile
				$(CC) -o abundanceN MassFunction/abundance.c -DNORM $(CFLAGS) $(LDFLAGS)  				
abundance:MassFunction/abundance.c Makefile
				$(CC) -o abundance MassFunction/abundance.c $(CFLAGS) $(LDFLAGS)  								
#image
subpart:Image/subparticles.c subfindread.o Makefile
				$(CC) -o subpart Image/subparticles.c subfindread.o $(CFLAGS) $(LDFLAGS)  
subcen:Image/subcen.c subfindread.o Makefile
				$(CC) -o subcen Image/subcen.c subfindread.o $(CFLAGS) $(LDFLAGS)  				
snap_low:Image/snap_lowres.c Makefile
				$(CC) -o snap_low Image/snap_lowres.c -I$(BTDIR) -include $(PARAM_LOW) $(LDFLAGS)	
sub_movie:Image/sub_movie.c Makefile
		$(CC) -o sub_movie Image/sub_movie.c $(CFLAGS) $(LDFLAGS)  						
				
snap_reacc:Image/snaps_reacc.c tree.o Makefile
				$(CC) -o snap_reacc Image/snaps_reacc.c tree.o $(CFLAGS) $(LDFLAGS)					
snap_energy:Image/snaps_energy.c tree.o Makefile
				$(CC) -o snap_energy Image/snaps_energy.c tree.o $(CFLAGS) $(LDFLAGS)		
snap_energy2:Image/snaps_energy2.c tree.o Makefile
				$(CC) -o snap_energy2 Image/snaps_energy2.c tree.o $(CFLAGS) $(LDFLAGS)													
snap:Image/snaps.c subfindread.o Makefile
				$(CC) -o snap Image/snaps.c  subfindread.o $(CFLAGS) $(LDFLAGS)			
snap2:Image/snap2.c subfindread.o treesearch.o tree.o Makefile
				$(CC) -o snap2 Image/snap2.c  subfindread.o treesearch.o tree.o $(CFLAGS) $(LDFLAGS)
snap3:Image/snap3.c subfindread.o Makefile
				$(CC) -o snap3 Image/snap3.c  subfindread.o $(CFLAGS) $(LDFLAGS)  												
snap_phase:Image/snaps_phase.c tree.o subfindread.o Makefile
				$(CC) -o snap_phase Image/snaps_phase.c tree.o subfindread.o $(CFLAGS) $(LDFLAGS)	
snap_phase_diff:Image/snaps_phase_diff.c subfindread.o Makefile
				$(CC) -o snap_phase_diff Image/snaps_phase_diff.c subfindread.o $(CFLAGS) $(LDFLAGS)		
snap_big:Image/snaps_big.c Makefile
				$(CC) -o snap_big Image/snaps_big.c  $(CFLAGS) $(LDFLAGS)				
snap_conf:Image/snaps_conf.c Makefile
				$(CC) -o snap_conf Image/snaps_conf.c $(CFLAGS) $(LDFLAGS)
#profile
M_c:Profile/M_c.c Makefile
				$(CC) -o M_c Profile/M_c.c $(CFLAGS) $(LDFLAGS)	
haloprof:Profile/haloprof.c linkedlist.c  Makefile #subfindread.o
				$(CC) -o haloprof Profile/haloprof.c  $(CFLAGS) $(LDFLAGS)	#subfindread.o		
gashaloprof:Profile/gashaloprof.c linkedlist.c Makefile
				$(CC) -o gashaloprof Profile/gashaloprof.c $(CFLAGS) $(LDFLAGS)	$(GASLIB)	
dynamprof:Profile/dynamprof.c linkedlist.c Makefile 
				$(CC) -o dynamprof Profile/dynamprof.c $(CFLAGS) $(LDFLAGS)	$(GASLIB)	
dynamprof.sph:Profile/dynamprof_sph.c linkedlist.c Makefile 
				$(CC) -o dynamprof.sph Profile/dynamprof_sph.c $(CFLAGS) $(LDFLAGS)	$(GASLIB)	
vel_tan:Profile/vel_tan.c  Makefile
				$(CC) -o vel_tan Profile/vel_tan.c $(CFLAGS) $(LDFLAGS) 
vdistr: Profile/veldistr.c tree.o  Makefile
				$(CC) -o vdistr Profile/veldistr.c tree.o $(CFLAGS) $(LDFLAGS) 
halocenter: Profile/halocenter.c treesearch.o tree.o Makefile
				$(CC) -o halocenter Profile/halocenter.c  treesearch.o tree.o $(CFLAGS) $(LDFLAGS)
CentralSub: Profile/CentralSub.c  Makefile
				$(CC) -o CentralSub Profile/CentralSub.c  $(CFLAGS) $(LDFLAGS)
#Profile/concentration/
NFW_fit: Profile/concentration/NFW_fit.c Profile/concentration/NFW_fun.c Profile/concentration/lnNFW_fun.c Makefile
	$(CC) -o NFW_fit Profile/concentration/NFW_fit.c $(CFLAGS) $(LDFLAGS) -lgsl -lgslcblas
NFW_halo: Profile/concentration/NFW_halo.c Profile/concentration/NFW_fun.c Profile/concentration/lnNFW_fun.c Makefile
	$(CC) -o NFW_halo Profile/concentration/NFW_halo.c $(CFLAGS) $(LDFLAGS) -lgsl -lgslcblas					
#Profile/geometrical							
tidalradius:Profile/geometrical/tidalradius.c linkedlist.c subfindread.o Makefile
				$(CC) -o tidalradius Profile/geometrical/tidalradius.c subfindread.o $(CFLAGS) $(LDFLAGS)
tidalradius_log:Profile/geometrical/tidalradius_log.c subfindread.o linkedlist.c Makefile
				$(CC) -o tidalradius_log Profile/geometrical/tidalradius_log.c subfindread.o $(CFLAGS) $(LDFLAGS)
tidalradius_aqua:Profile/geometrical/tidalradius_aqua.c linkedlist.c subfindread.o Makefile
				$(CC) -o tidalradius_aqua Profile/geometrical/tidalradius_aqua.c subfindread.o $(CFLAGS) $(LDFLAGS)
fofpot: Profile/geometrical/fofpot.c tree.o  Makefile
				$(CC) -o fofpot Profile/geometrical/fofpot.c tree.o $(CFLAGS) $(LDFLAGS) 		
grpcen: Profile/geometrical/grpcen.c Makefile
				$(CC) -o grpcen Profile/geometrical/grpcen.c $(CFLAGS) $(LDFLAGS) 					
subshape: Profile/geometrical/subshape.c subfindread.o Makefile
				$(CC) -o subshape Profile/geometrical/subshape.c subfindread.o $(CFLAGS)  $(GSLINC) $(LDFLAGS) $(GSLLIB)				
notts_shape: Profile/geometrical/notts_shape.c Makefile
				$(CC) -o notts_shape Profile/geometrical/notts_shape.c $(CFLAGS) $(GSLINC) $(LDFLAGS) $(GSLLIB) 				
#Profile/sub				
subprof:Profile/sub/subprofile.c linkedlist_DM.c Makefile
				$(CC) -o subprof Profile/sub/subprofile.c $(CFLAGS) $(LDFLAGS) -USUBFIND_DIR
subfindprof:Profile/sub/subprofile.c linkedlist_DM.c subfindread.o Makefile
				$(CC) -o subfindprof Profile/sub/subprofile.c subfindread.o $(CFLAGS) $(LDFLAGS) 				
Vmax:Profile/sub/Vmax.c tree.o Makefile
				$(CC) -o Vmax Profile/sub/Vmax.c $(CFLAGS) $(LDFLAGS) tree.o -USUBFIND_DIR		
SubInSub:Profile/sub/SubInSub.c Makefile
				$(CC) -o SubInSub Profile/sub/SubInSub.c $(CFLAGS) $(LDFLAGS) 
plotSubInSub:Profile/sub/plotSubInSub.c Makefile
				$(CC) -o plotSubInSub Profile/sub/plotSubInSub.c $(CFLAGS) $(LDFLAGS) -UDYNAMICAL
plotSubInSub.dyn:Profile/sub/plotSubInSub.c Makefile
				$(CC) -o plotSubInSub.dyn Profile/sub/plotSubInSub.c $(CFLAGS) $(LDFLAGS) -DDYNAMICAL								
										
#Profile/xray/				
Lxsub: Profile/xray/Lxsub.c linkedlist.c Makefile
				$(CC) -o Lxsub Profile/xray/Lxsub.c $(CFLAGS) $(LDFLAGS) -lgsl -lgslcblas
gasprof: Profile/xray/gasprofile.c linkedlist.c Makefile
				$(CC) -o gasprof Profile/xray/gasprofile.c $(CFLAGS) $(LDFLAGS) -lgsl -lgslcblas				
gasprof2: Profile/xray/gasprofile2.c linkedlist.c Makefile
				$(CC) -o gasprof2 Profile/xray/gasprofile2.c $(CFLAGS) $(LDFLAGS) -lgsl -lgslcblas
gprf_mkl:Profile/xray/gasprofile4mekl.c linkedlist.c Makefile
				$(CC) -o gprf_mkl Profile/xray/gasprofile4mekl.c $(CFLAGS) $(LDFLAGS) 
				
#Check/				
pop:Check/populations.c Makefile
	$(CC) -o pop Check/populations.c $(CFLAGS) $(LDFLAGS)
dest:Check/dest.c Makefile
				$(CC) -o dest Check/dest.c $(CFLAGS) $(LDFLAGS)				
follow_group:Check/follow_group.c Makefile
	$(CC) -o follow_group Check/follow_group.c $(CFLAGS) $(LDFLAGS)				
loadview:Check/loadview.c Makefile
				$(CC) -o loadview Check/loadview.c $(CFLAGS) $(LDFLAGS)	-g
forge_srccat:Check/forge_srccat.c Makefile
				$(CC) -o forge_srccat Check/forge_srccat.c $(CFLAGS) $(LDFLAGS)	-g	
forge_srccat_static:Check/forge_srccat_static.c Makefile
				$(CC) -o forge_srccat_static Check/forge_srccat_static.c $(CFLAGS)$(LDFLAGS) -g					
allsubs:Check/allsubs.c Makefile
				$(CC) -o allsubs Check/allsubs.c $(CFLAGS) $(LDFLAGS)	-g
redshift:Check/redshift.c Makefile
				$(CC) -o redshift Check/redshift.c $(CFLAGS) $(LDFLAGS)	
byteorders:Check/byteorders.c Makefile
				$(CC) -o byteorders Check/byteorders.c $(CFLAGS) $(LDFLAGS)																
testprop:Check/testprop.c Makefile
				$(CC) -o testprop Check/testprop.c $(CFLAGS) $(LDFLAGS)	
loadhead:Check/loadhead.c Makefile
				$(CC) -o loadhead Check/loadhead.c $(CFLAGS) $(LDFLAGS)	
overlap:Check/overlap.c Makefile
				$(CC) -o overlap Check/overlap.c $(CFLAGS) $(LDFLAGS)	
find_a_halo:Check/find_a_halo.c Makefile
				$(CC) -o find_a_halo Check/find_a_halo.c $(CFLAGS) $(LDFLAGS)
match_subcat:Check/match_subcat.c subfindread.o Makefile
				$(CC) -o match_subcat Check/match_subcat.c subfindread.o $(CFLAGS) $(LDFLAGS)				
back_binding:Check/back_binding.c subfindread.o Makefile
				$(CC) -o back_binding Check/back_binding.c subfindread.o $(CFLAGS) $(LDFLAGS)
saveNotts:Check/saveNotts.c Makefile
				$(CC) -o saveNotts Check/saveNotts.c $(CFLAGS) $(LDFLAGS) 
				
saveNotts_src:Check/saveNotts_src.c Makefile
				$(CC) -o saveNotts_src Check/saveNotts_src.c $(CFLAGS) $(LDFLAGS) 

notts: saveNotts Vmax
				./Vmax
				./saveNotts

subfindread: subfindread.c Makefile
			$(CC) -o subfindread subfindread.c $(CFLAGS) $(LDFLAGS) -DSTANDALONE

subfindread.o: subfindread.c Makefile

tree.o: $(BTDIR)/tree.c	Makefile
				$(CC) -c $(BTDIR)/tree.c $(CFLAGS) 
treesearch.o: $(BTDIR)/treesearch.c	Makefile
				$(CC) -c $(BTDIR)/treesearch.c $(CFLAGS) 
								
clean:
	rm -f *.o a.out massfun massfun_cosm massfun_plot subpart subprof haloprof gashaloprof dynamprof dynamprof.sph
	rm -f  tidalradius tidalradius_log tidalradius_aqua subshape trace_bin pick_subsnap orbitparam trace_bin_ext
	rm -f steller steller2 Lxsub gasprof gprf_mkl dest loadview loadhead snap_low snap vdistr fofpot grpcen
	rm -f trace_bin_rev subfindread.o redshift byteorders biggestsub match_subcat back_binding overlap check_centers
	rm -f NFW_fit NFW_halo quasiMmax pop pick_halosample massfun_pick abundanceN snap_phase snap_phase_diff
	rm -f snap2 snap3 Vmax halocenter subfindread SubInSub plotSubInSub plotSubInSub.dyn saveNotts SubHist Decay
