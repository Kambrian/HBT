include ../../Makefile.runs
BTDIR=../..
GASDIR=../../GasCat

USRINCL=-I /home/kambrain/opt/include
USRLIB=-L /home/kambrain/opt/lib
CFLAGS=-I$(BTDIR) -I$(BTDIR)/$(IODIR) -include $(PARAM) -openmp $(USRINCL) 
LDFLAGS= -L$(BTDIR)/$(IODIR) -l$(RUN_NUM)io $(FTNLIB)  -openmp $(USRLIB)
GASLIB=-L$(GASDIR) -lgasio

VCOM_PATCH =VCoM.patch
OBJS_PATCH =VCoM.patch.o
$(VCOM_PATCH): $(OBJS_PATCH) tree.o $(OBJS_FTN)
	$(CC) $(OBJS_PATCH) tree.o $(OBJS_FTN) $(LDFLAGS)    -o  $(VCOM_PATCH) 

$(OBJS_PATCH): $(INCL)

tree.o: $(BTDIR)/tree.c	Makefile
				$(CC) -c $(BTDIR)/tree.c $(CFLAGS) -g

SAVE_PATCH =save.patch
OBJS_SAVE =savesub.patch.o
$(SAVE_PATCH): $(OBJS_SAVE)  $(OBJS_FTN)
	$(CC) $(OBJS_SAVE)   $(OBJS_FTN) $(LDFLAGS)     -o  $(SAVE_PATCH) 

$(OBJS_SAVE): $(INCL)

SAVE_PATCH2 =save.patch2
OBJS_SAVE2 =savesub.patch2.o
$(SAVE_PATCH2): $(OBJS_SAVE2)  $(OBJS_FTN)
	$(CC) $(OBJS_SAVE2)   $(OBJS_FTN) $(LDFLAGS)     -o  $(SAVE_PATCH2) 

$(OBJS_SAVE2): $(INCL)

sp2pro.patch:sp2pro.patch.c Makefile
	$(CC) -o sp2pro.patch sp2pro.patch.c $(CFLAGS) $(LDFLAGS) 
	
clean: 
	rm -f *.patch *.patch2
	rm -f *.o
	rm -f core.* core
