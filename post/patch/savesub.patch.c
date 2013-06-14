/*to convert v7.5 subcatalogue to v7.6 format*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

void load_sub_catalogue_75(int Nsnap, SUBCATALOGUE *Cat, char *SubCatPath);

int main(int argc,char **argv)
{
char outputdir[1024];
char inputdir[1024];
SUBCATALOGUE SubCat;
int SnapshotNum;
int snapbegin,snapend;

if(argc!=3)
{
	printf("usage: %s [snap_begin] [snap_end]\n",argv[0]);
	return 1;
}
snapbegin=atoi(argv[1]);
snapend=atoi(argv[2]);

sprintf(inputdir,"%s/patched",SUBCAT_DIR);
sprintf(outputdir,"%s/patched_new",SUBCAT_DIR);
for(SnapshotNum=snapbegin;SnapshotNum<=snapend;SnapshotNum++)
{
	load_sub_catalogue_75(SnapshotNum,&SubCat,inputdir);
	save_sub_catalogue(SnapshotNum,&SubCat,outputdir);
	erase_sub_catalogue(&SubCat);
}
return 0;
}

void load_sub_catalogue_75(int Nsnap, SUBCATALOGUE *Cat, char *SubCatPath)
{
FILE *fd;
char buf[1024];
int i;

  sprintf(buf, "%s/subcat_%03d", SubCatPath, Nsnap);
  if(!(fd = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fread(&Cat->Ngroups, sizeof(int), 1, fd);
        Cat->GrpOffset_Sub=mymalloc(sizeof(int)*Cat->Ngroups);   
	Cat->GrpLen_Sub=mymalloc(sizeof(int)*Cat->Ngroups);
  fread(&Cat->Nsubs, sizeof(int), 1, fd);
	Cat->SubLen=mymalloc(sizeof(int)*Cat->Nsubs);
	Cat->SubOffset=mymalloc(sizeof(int)*Cat->Nsubs);
	Cat->SubRank=mymalloc(sizeof(int)*Cat->Nsubs);
	Cat->HaloChains=mymalloc(sizeof(struct Chain_data)*Cat->Nsubs);
	Cat->sub_hierarchy=mymalloc(sizeof(struct Hierarchy)*Cat->Nsubs);
	Cat->Property=mymalloc(sizeof(struct SubProperty)*Cat->Nsubs);
	Cat->PSubArr=mymalloc(sizeof(int *)*Cat->Nsubs);
  fread(&Cat->Nids, sizeof(int), 1, fd);
  fread(Cat->GrpLen_Sub, sizeof(int), Cat->Ngroups, fd);
  fread(Cat->GrpOffset_Sub,sizeof(int), Cat->Ngroups, fd);
  fread(Cat->SubLen, sizeof(int), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(int), Cat->Nsubs, fd);
  fread(Cat->SubRank,sizeof(int), Cat->Nsubs, fd);
  fread(Cat->HaloChains,sizeof(struct Chain_data), Cat->Nsubs, fd);
  for(i=0;i<Cat->Nsubs;i++)
    fread(Cat->Property[i].CoM,sizeof(float), 3, fd);
  for(i=0;i<Cat->Nsubs;i++)
    fread(Cat->Property[i].VCoM,sizeof(float), 3, fd);
  for(i=0;i<Cat->Nsubs;i++)
    fread(&(Cat->Property[i].Pot),sizeof(float), 1, fd);
  for(i=0;i<Cat->Nsubs;i++)
	fread(&(Cat->Property[i].Kin),sizeof(float), 1, fd);
  for(i=0;i<Cat->Nsubs;i++)
	fread(Cat->Property[i].AM,sizeof(float), 3, fd);
  fread(Cat->sub_hierarchy,sizeof(struct Hierarchy),Cat->Nsubs,fd);
for(i=0;i<Cat->Nsubs;i++)
{
	Cat->PSubArr[i]=mymalloc(sizeof(int)*Cat->SubLen[i]);
	fread(Cat->PSubArr[i], sizeof(int), Cat->SubLen[i], fd);
}
fread(&Cat->Nbirth,sizeof(int),1,fd);
fread(&Cat->NQuasi,sizeof(int),1,fd);
fread(&Cat->Ndeath,sizeof(int),1,fd);
fread(&Cat->Nsplitter,sizeof(int),1,fd);

  fclose(fd);
}
