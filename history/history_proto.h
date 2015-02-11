#define GetHostID(Member) (((Member)==NULL)?-2:((Member)->HostID))
extern HBTReal load_halo_size(HALOSIZE *halosize,HBTInt Ngroups,HBTInt Nsnap);
extern HBTInt load_halo_concentration(HBTReal *halocon,HBTInt Nsnap);
extern HBTInt read_mainsubid(HBTInt Nsnap,HBTInt hostid);
extern HBTInt read_Ngroups(HBTInt Nsnap);
extern HBTInt read_subpos(HBTInt Nsnap,HBTxyz *(*snappos));
extern HBTInt read_subvel(HBTInt Nsnap,HBTxyz *(*snapvel));
extern size_t load_brfcat(HBTInt Nsnap,BRFCAT *Cat);
extern void erase_brfcat(BRFCAT *Cat);

extern void save_sub2hist(HBTInt Nsnap,HBTInt *Sub2Hist,HBTInt Nsubs,char *subcatdir);
extern void load_sub2hist(HBTInt Nsnap,HBTInt **P2sub2hist,HBTInt *P2Nsubs,char *subcatdir);
extern void save_evocat_pre(EVOLUTIONCAT_Pre *EvoCat,char *subcatdir);
extern void load_evocat_pre(EVOLUTIONCAT_Pre *EvoCat,char *subcatdir);
extern void free_evocat_pre(EVOLUTIONCAT_Pre *ECat);
extern void load_history_pre(EVOLUTIONCAT *EvoCat,char *subcatdir);
//same as load_evocat_pre, except for using datatype EVOLUTIONCAT rather than EVOLUTIONCAT_Pre

extern void save_evocat_rev(EVOLUTIONCAT *EvoCat,char *subcatdir);
extern void load_evocat_rev(EVOLUTIONCAT *EvoCat,char *subcatdir);

extern SubNode *GetMember(EVOLUTIONCAT *EvoCat,HBTInt HistID,HBTInt Nsnap);

extern void create_historyshards(struct HistoryShards *HistoryShard,HBTInt NumHist);
extern void save_historyshards(struct HistoryShards *HistoryShard);
extern void load_historyshards(struct HistoryShards *HistoryShard);

extern void fresh_MBDID2Index(MBDCATALOGUE *MbdCat);
extern void save_mbd_catalogue(HBTInt Nsnap, MBDCATALOGUE *MbdCat);
extern void load_mbd_catalogue(HBTInt Nsnap, MBDCATALOGUE *MbdCat);
extern void free_mbd_catalogue(MBDCATALOGUE *MbdCat);
