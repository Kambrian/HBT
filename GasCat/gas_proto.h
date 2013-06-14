//gas_io.c
extern void load_gas_data(HBTInt Nsnap, char *SnapPath);
extern void free_gas_data();
extern void save_gashalocat(HBTInt Nsnap, GASHALOCAT *Cat,char *gaspath);
extern void load_gashalocat(HBTInt Nsnap,GASHALOCAT *Cat,char *GrpPath);
extern void free_gashalocat(GASHALOCAT *A);
extern HBTInt prepare_gasind2halo(GASHALOCAT *A);
extern void fresh_gasID2Index(const void *src, HBTInt src_len);
extern void save_gassrccat(HBTInt Nsnap,GASSRCCAT *Cat,char *outputdir);
extern void load_gassrccat(HBTInt Nsnap,GASSRCCAT *Cat,char *inputdir);
extern void save_gassubcat(HBTInt Nsnap,GASSUBCAT *Cat,char *outputdir);
extern void load_gassubcat(HBTInt Nsnap,GASSUBCAT *Cat,char *outputdir);
extern void create_gassrccat(GASSRCCAT *GSrcCat);
extern void free_gassrccat(GASSRCCAT *GSrcCat);
extern void erase_gassrccat(GASSRCCAT *GSrcCat);
extern void create_gassubcat(GASSUBCAT *GSubCat);
extern void free_gassubcat(GASSUBCAT *GSubCat);
extern void erase_gassubcat(GASSUBCAT *GSubCat);
extern void complete_N_save_gas(GASSRCCAT *GSrcCat,GASSUBCAT *GSubCat,HBTInt SnapshotNum,char *gasdir);

//gas_binding.c
extern HBTInt unbindgas(HBTInt *P2Len,HBTInt **P2PIndex, struct GasProperty *Prop, HBTInt SubLen, HBTInt *SubArr,HBTReal SubCoM[3],HBTReal SubVCoM[3]); /*P2Len=&Len, *P2PIndex=PIndex, 
																*where PIndex[Len] is the array of size Len; 
																* both will be updated after unbind*/
extern void make_gassrccat(GASSRCCAT *GSrcCat,GASSUBCAT *GSubCat,GASHALOCAT *GCat,SUBCATALOGUE *SubCat,HBTInt Nsnap);
extern void dump_gassrccat(GASSRCCAT *GSrcCatTo,GASSRCCAT *GSrcCatFrom);
extern void migrate_ghalosrc(GASSRCCAT *GSrcCat, HBTInt subid, GASHALOCAT *Cat,HBTInt hostid);
extern void migrate_gsubsrc(GASSRCCAT *GSrcCatTo, HBTInt subid, GASSRCCAT *GSrcCatFrom,HBTInt proid);
extern void migrate_gsplsrc(GASSRCCAT *GSrcCatTo, HBTInt subid, GASSRCCAT *GSrcCatFrom,HBTInt proid,HBTInt **PSubArr);
extern void makell_gas();
extern void collect_gas_particles(HBTInt subid,SUBCATALOGUE*SubCat,HBTInt *P2GasSrcLen,HBTInt **P2GasSrc);
