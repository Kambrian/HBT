include ../../Makefile.runs
BTDIR=../..

# you need to link against a library compiled with HALO_PARA undefined for cosmological simu
CFLAGS=  -openmp -include $(PARAM) -I $(BTDIR) -I $(BTDIR)/$(IODIR) -I. 
LDFLAGS= -openmp -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB) tree.o
XLDFLAGS= -openmp -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB)
DEPS= binding.c tree.o

follow_fof:follow_fof.mean follow_fof.minpotH follow_fof.potW follow_fof.core follow_fof.adpt follow_fof.local
	  

follow_fof.mean: follow_fof.c $(DEPS)
	$(CC) follow_fof.c -DFOLLOW_MEAN $(CFLAGS) $(LDFLAGS) -o follow_fof.mean
follow_fof.minpotH: follow_fof.c $(DEPS)
	$(CC) follow_fof.c -DFOLLOW_MINPOTH $(CFLAGS) $(LDFLAGS) -o follow_fof.minpotH
follow_fof.potW: follow_fof.c $(DEPS)
	$(CC) follow_fof.c -DFOLLOW_POTW $(CFLAGS) $(LDFLAGS) -o follow_fof.potW
follow_fof.core: follow_fof.c $(DEPS)
	$(CC) follow_fof.c -DFOLLOW_CORE $(CFLAGS) $(LDFLAGS) -o follow_fof.core
	
follow_main.mean: follow_fof.c $(DEPS)
	$(CC) follow_fof.c -DFOLLOW_MEAN -DSUBID_INPUT $(CFLAGS) $(LDFLAGS) -o follow_main.mean
follow_main.minpotH: follow_fof.c $(DEPS)
	$(CC) follow_fof.c -DFOLLOW_MINPOTH -DSUBID_INPUT $(CFLAGS) $(LDFLAGS) -o follow_main.minpotH
follow_main.potW: follow_fof.c $(DEPS)
	$(CC) follow_fof.c -DFOLLOW_POTW -DSUBID_INPUT $(CFLAGS) $(LDFLAGS) -o follow_main.potW
follow_main.core: follow_fof.c $(DEPS)
	$(CC) follow_fof.c -DFOLLOW_CORE -DSUBID_INPUT $(CFLAGS) $(LDFLAGS) -o follow_main.core	

follow_fof.adpt:follow_fof_adptive.c $(DEPS)	
	$(CC) follow_fof_adptive.c  $(CFLAGS) $(LDFLAGS) -o follow_fof.adpt

follow_fof.stat:follow_fof_stat.c $(DEPS)	
	$(CC) follow_fof_stat.c  $(CFLAGS) $(LDFLAGS) -o follow_fof.stat
	
follow_fof.local:follow_fof_local.c $(DEPS)	
	$(CC) follow_fof_local.c  $(CFLAGS) $(LDFLAGS) -o follow_fof.local
	
follow_sub.local:follow_sub_local.c $(DEPS)	
	$(CC) follow_sub_local.c  $(CFLAGS) $(LDFLAGS) -o follow_sub.local

follow_sub.pos:follow_sub_pos.c $(DEPS)	
	$(CC) follow_sub_pos.c  $(CFLAGS) $(LDFLAGS) -o follow_sub.pos	

unbinding_curve:unbinding_curve.c $(DEPS)
	$(CC) unbinding_curve.c  $(CFLAGS) $(LDFLAGS) -o unbinding_curve
	
checkcenter:checkcenter.c 	
	$(CC) checkcenter.c  $(CFLAGS) $(XLDFLAGS) -o checkcenter
	
tree.o: $(BTDIR)/tree.c	Makefile
				$(CC) -c $(BTDIR)/tree.c $(CFLAGS) 
					
clean:
	rm -f follow_fof.core follow_fof.mean follow_fof.minpotH follow_fof.potW
	rm -f follow_main.core follow_main.mean follow_main.minpotH follow_main.potW 
	rm -f follow_fof.local follow_sub.local follow_fof.adpt follow_fof.stat checkcenter
	rm -f *.o core.* 
