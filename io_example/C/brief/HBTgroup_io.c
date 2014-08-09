typedef int HBTInt;
#define HBTIFMT "%d"
typedef float HBTReal;

typedef struct 
{
HBTInt Ngroups; //number of groups
HBTInt Nids; //number of particles
HBTInt *Len; //Length of each group
HBTInt *Offset; //Starting position of each group within PID[]
HBTInt *PID; //particle IDs of all the groups
}CATALOGUE;

void load_group_catalogue_HBT(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{
  FILE *fd;
  HBTInt i,NFiles;
  char buf[1024];
  
#ifdef SNAPLIST
Nsnap=snaplist[Nsnap];
#endif

  sprintf(buf, "%s/group_tab_%03d",GrpPath,(int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&NFiles, sizeof(HBTInt), 1, fd); //useless, just to conform to GADGET's data format
  
  Cat->Len=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  fread(Cat->Len, sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->Offset,sizeof(HBTInt), Cat->Ngroups, fd);
  fclose(fd);

fprintf(logfile,"Snap="HBTIFMT" Ngroups="HBTIFMT" Nids="HBTIFMT"\n", Nsnap, Cat->Ngroups, Cat->Nids);
  sprintf(buf, "%s/group_ids_%03d", GrpPath, (int)Nsnap);
  myfopen(fd,buf,"r");

  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
  fread(&NFiles, sizeof(HBTInt), 1, fd);

  Cat->PID=mymalloc(sizeof(HBTInt)*Cat->Nids);
  fread(Cat->PID,sizeof(HBTInt),Cat->Nids,fd);
    
  fclose(fd);
}
