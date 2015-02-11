#define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){ fprintf(logfile,"Error opening file '%s'\n",filename);	fflush(logfile); exit(1);	}
#ifdef PERIODIC_BDR
#define NEAREST(x) (((x)>BOXHALF)?((x)-BOXSIZE):(((x)<-BOXHALF)?((x)+BOXSIZE):(x)))
			/*this macro can well manipulate boundary condition because 
			* usually a halo is much smaller than boxhalf
			* so that any distance within the halo should be smaller than boxhalf */
#endif
#define get_bit(x,k) (((x)&(1<<k))>>k)

//mymath.c
extern HBTInt Dmax_of_vec(double *vec,HBTInt Len);//what if the max is not unique???????
extern HBTInt Fmax_of_vec(float *vec,HBTInt Len);//what if the max is not unique???????
extern HBTInt max_of_vec(HBTReal *vec,HBTInt Len);//what if the max is not unique???????
extern HBTInt Dmin_of_vec(double *vec,HBTInt Len);//what if the max is not unique???????
extern HBTInt Fmin_of_vec(float *vec,HBTInt Len);//what if the max is not unique?,return first min
extern HBTInt min_of_vec(HBTReal *vec,HBTInt Len);//what if the max is not unique?,return first min
extern HBTInt insert_to_array(HBTInt a, HBTInt len, HBTInt *ascd_arr);//insert into an ascending arr
extern HBTInt prepare_ind2halo(CATALOGUE *A);
extern HBTInt * prepare_sub2halo(SUBCATALOGUE *SubCat);
extern HBTInt check_dup(HBTInt l1, HBTInt *p1,HBTInt l2, HBTInt *p2);
extern HBTReal position_modulus(HBTReal x);
extern HBTReal distance(HBTReal x[3],HBTReal y[3]);
extern void center_of_mass(HBTReal CoM[3], HBTInt *PInd, HBTInt np, HBTReal PPos[][3]);
extern HBTReal comoving_virial_radius(HBTInt mass);
extern void *mymalloc(size_t n);
extern void myfree(void *mem);
extern void swap_Nbyte(void *data2swap,size_t nel,size_t mbyte);
extern size_t fread_swap(void *buf,size_t Nsize,size_t Nbuf,FILE *fp, int FlagByteSwap);
extern void spherical_basisD(double dx[3],double er[3],double et[3],double ef[3]);
extern HBTReal vec_prod(HBTReal *a,HBTReal *b,HBTInt dim);
extern void vec_cross(HBTReal a[3],HBTReal b[3],HBTReal c[3]);
extern HBTInt try_readfile(char *filename);
extern HBTInt count_lines(char * filename);
extern HBTReal psort(HBTInt k, HBTInt numel, HBTReal arr[]);
extern HBTReal * GetArrPos(HBTInt pid,void *data);
extern HBTReal * GetSubCoM(HBTInt subid,void *data);
extern HBTInt linklist_round_gridid(HBTInt i,HBTInt ndiv);
extern HBTInt linklist_shift_gridid(HBTInt i,HBTInt ndiv);
extern HBTInt linklist_fix_gridid(HBTInt i, LINKLIST *ll);
extern HBTInt linklist_get_hoc(LINKLIST *ll, HBTInt i,HBTInt j,HBTInt k);
extern HBTInt linklist_get_hoc_safe(LINKLIST *ll, HBTInt i,HBTInt j,HBTInt k);
extern void make_linklist(LINKLIST *ll, HBTInt np,HBTInt ndiv, 
						void *PosData, AccessPosFunc *GetPos, HBTInt UseFullBox);
extern void free_linklist(LINKLIST *ll);
extern HBTInt * linklist_search_sphere(LINKLIST *ll, HBTReal radius, 
							HBTReal searchcenter[3], HBTInt *N_max_and_found);

//user_IO.c
extern int check_grpcat_byteorder(char *filename);
extern void load_particle_header_into(HBTInt Nsnap, char *SnapPath, IO_HEADER *h);
extern void load_particle_header(HBTInt Nsnap,char *SnapPath);
extern void load_particle_data(HBTInt Nsnap,char *SnapPath);
extern void load_particle_data_bypart(HBTInt Nsnap, char *SnapPath, unsigned char loadflags);
extern void free_particle_data();
extern IDatInt *load_PIDs_Sorted();
extern void load_group_catalogue(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath);//the default loader; wrapper for v2 or v3
extern void load_group_catalogue_millimill(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath);//even older format, used by MilliMillennium
extern void load_group_catalogue_v2(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath);//old format
extern void load_group_catalogue_v3(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath);//gadget-3 format
extern void load_group_catalogue_v4(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath);//gadget-4 format
extern void free_catalogue(CATALOGUE *A);
//sub_IO.c
extern void complete_N_save(SUBCATALOGUE *SubCatB,SRCCATALOGUE *SrcCatB,HBTInt SnapshotNum,char *outputdir);
extern void fill_PIDHash2();//general bsearch hash
extern void free_PIDHash2();
extern HBTInt lookup_PIDHash2(HBTInt PID);
extern void fill_PIDHash1();//simple linear hash (PIndex array)
extern void free_PIDHash1();
extern HBTInt lookup_PIDHash1(HBTInt PID);
extern void fill_PIDHash();//wrapper for the above two hash functions
extern void free_PIDHash();
#if defined PID_NEED_HASH && !defined HBTPID_RANKSTYLE
#define lookup_ID2Ind lookup_PIDHash2
#else
#define lookup_ID2Ind lookup_PIDHash1
#endif

#ifdef PID_ORDERED
#define lookup_Ind2ID(i) ((i)<0?-1:((i)+1))
#else
#define lookup_Ind2ID(i) ((i)<0?-1:Pdat.PID[i])
#endif

extern void fresh_ID2Index(const void *src, HBTInt src_len);
extern void create_sub_cat(SUBCATALOGUE *SubCat);
extern void create_src_cat(SRCCATALOGUE *SrcCat);
extern void save_sub_catalogue(HBTInt Nsnap, SUBCATALOGUE *Cat,char *SubCatPath);
extern void load_sub_catalogue(HBTInt Nsnap, SUBCATALOGUE *Cat, char *SubCatPath);
extern void save_src_catalogue(HBTInt Nsnap, SRCCATALOGUE *Cat,char *SrcCatPath);
extern void load_src_catalogue(HBTInt Nsnap, SRCCATALOGUE *Cat,char *SrcCatPath);
extern void load_sp2pro(HBTInt Nsnap_dest,HBTInt *P2Npro,HBTInt *P2Nsplitter_dest,HBTInt **P2sp2pro, char *SubcatPath);
extern void free_sp2pro(HBTInt *sp2pro,HBTInt Npro,HBTInt Nsplitter_dest);
extern void free_src_catalogue(SRCCATALOGUE *SrcCat);
extern void free_sub_catalogue(SUBCATALOGUE *SubCat);
extern void erase_sub_catalogue(SUBCATALOGUE *SubCat);
extern void erase_src_catalogue(SRCCATALOGUE *SrcCat);
extern void fix_spfile(FILE **p2fp,char *buf,HBTInt Snap);
extern void load_pro2dest(HBTInt Nsnap_pro,HBTInt **P2pro2dest,HBTInt *P2Nsubs,char *destdir);
extern void free_pro2dest(HBTInt *pro2dest);
extern void save_pro2dest(HBTInt Nsnap_pro,HBTInt *pro2dest,HBTInt Nsubs,char *destdir);
extern void load_sub_table(HBTInt Nsnap, SUBCATALOGUE *Cat, char *SubCatPath);//only load global properties, without filling particles
extern void free_sub_table(SUBCATALOGUE *SubCat);//corresponding to load_sub_table
extern void save_group_catalogue_HBT(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath);
extern void load_group_catalogue_HBT(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath);

//hierarchy.c
extern void transfer_subcat(SUBCATALOGUE *SubCatTo,SUBCATALOGUE *SubCatFrom);
extern void transfer_srccat(SRCCATALOGUE *SrcCatTo,SRCCATALOGUE *SrcCatFrom);
extern void break_out_sub_family(HBTInt mainsubID,SUBCATALOGUE *SubCat);
extern void markcross_N_kickdeath_recursive(HBTInt mainsubID,SUBCATALOGUE *SubCat);
extern void kickdeath_recursive(HBTInt mainsubID,SUBCATALOGUE *SubCat);
extern void PARAsplit_srccat(CATALOGUE* Cat,SRCCATALOGUE *SrcCat, HBTInt SnapshotNum);
extern void PARAmake_srcsub(SUBCATALOGUE *SubCatA,SUBCATALOGUE *SubCatB,SRCCATALOGUE *SrcCatA,SRCCATALOGUE *SrcCatB);
extern void PARAinit_dessub(SUBCATALOGUE *SubCatA,SUBCATALOGUE *SubCatB, struct LinkInfo *linkinfo);
extern void migrate_sub(SUBCATALOGUE *SubCatA,SUBCATALOGUE *SubCatB,HBTInt proSubID,HBTInt SubRank,HBTInt HostID,HBTInt *pro2dest);
extern void migrate_src(SRCCATALOGUE *SrcCatA, SRCCATALOGUE *SrcCatB,HBTInt proSubID,HBTInt desSubID);
extern void mask_mainsub(CATALOGUE *CatB,SUBCATALOGUE *SubCatB,HBTInt proSubID,HBTInt desID,HBTInt son);
extern void mask_mainsub_check(CATALOGUE *CatB,SUBCATALOGUE *SubCatB,HBTInt proSubID,HBTInt desID,HBTInt son);
extern void mask_mainsrc(CATALOGUE *CatB,SRCCATALOGUE *SrcCatB,HBTInt desID,HBTInt desSubID,HBTInt GrpLen_Sub);
extern void restore_mainsub(SUBCATALOGUE *SubCatA,SUBCATALOGUE *SubCatB,HBTInt proSubID,HBTInt desSubID);
extern void narrow_srccat(SRCCATALOGUE *SrcCat,SUBCATALOGUE *SubCat, HBTInt subid);
extern void PARAinit_mask(CATALOGUE *Cat,HBTInt mtype);
extern void mask_src_recursive(HBTInt subid,SUBCATALOGUE *SubCat, SRCCATALOGUE *SrcCat,char *HaloMaskSrc);
//tree.c
extern void update_internal_nodes(HBTInt no,HBTInt sib,double len,HBTInt * PIndex,HBTReal PPos[][3]);
extern HBTInt maketree(HBTInt halolen,HBTInt * PIndex,HBTReal PPos[][3]);
extern double tree_treeevaluate_potential(HBTReal targetPos[3], HBTInt *PIndex,HBTReal PPos[][3]);
extern size_t tree_tree_allocate(size_t maxnodes, size_t maxpart);	
extern void tree_tree_free(void);
	
//binding_minpot.c
extern HBTInt unbind(HBTInt *P2Len,HBTInt **P2PIndex, struct SubProperty *Prop,HBTInt *P2Len_removed, HBTInt **P2PIndex_removed,HBTReal CoreFrac);
extern void unbind_sub_recursive(HBTInt mainsubID,HBTInt *P2Len_removed,HBTInt **P2PIndex_removed,SUBCATALOGUE *SubCat,SRCCATALOGUE *SrcCat);

//treesearch.c
extern void tree_count_bin(HBTReal cen[3],HBTReal *edges,HBTInt nbin,HBTInt* bin_count,HBTInt *PIndex);
extern HBTReal cutting_sphere_radius_cube(HBTReal cen[3],HBTReal cubecen[3],HBTReal cubelen_half);
extern HBTReal enclosing_sphere_radius_cube(HBTReal cen[3],HBTReal cubecen[3],HBTReal cubelen_half);
extern HBTInt binsert_asc(HBTReal *edges,HBTInt nbin,HBTReal target);
extern void distance_point2cube(HBTReal cen[3],HBTReal cubecen[3],HBTReal cubelen_half,HBTReal dist[3]);
extern HBTReal guess_ngb_range(HBTInt NumNgb);
extern HBTReal guess_ngb_range_halo(HBTInt NumNgb);
extern HBTInt treesearch_nearest(HBTReal cen[3],HBTReal hguess,HBTInt *PIndex,HBTReal PPos[][3]);
extern HBTInt treesearch_sphere(HBTReal searchcenter[3], HBTReal radius,HBTInt *PIndex,HBTReal PPos[][3]);
extern double sph_density(HBTReal cen[3],HBTReal *p2hguess,HBTInt *PIndex,HBTReal PPos[][3]);
extern HBTInt treesearch_linkgrp(HBTReal radius, HBTInt PIndex[], struct GroupData *GrpData);
extern HBTInt treesearch_infect_particles(HBTInt seed, HBTInt grpid, 
											struct ParticleGroup *GrpTags, HBTReal PPos[][3]);







