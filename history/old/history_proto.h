extern void save_history(int NumHist,SubHist *ClusterHistory,char *subcatdir);
extern void load_history(int *P2NumHist,SubHist **P2ClusterHistory,char *subcatdir);
extern void save_sub2hist(int Nsnap,int *Sub2Hist,int Nsubs,char *subcatdir);
extern void load_sub2hist(int Nsnap,int **P2sub2hist,int *P2Nsubs,char *subcatdir);
