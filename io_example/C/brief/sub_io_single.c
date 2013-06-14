//this file is an assemble of io routines related to subhalo catalogue
//you may want to link against BT library rather than dive into this file

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* for the correct definition of HBTInt and HBTReal,
 * check param/paramRUNNAME.h  to see if HBT_INT8 and HBT_REAL8 are defined.*/
typedef long long HBTInt;
typedef double HBTReal;

FILE * logfile;

struct Chain_data 
{
HBTInt ProSubID;
HBTInt HostID;
};

struct Hierarchy
{
HBTInt nibs;	//its host sub
HBTInt pre;	//its previous sibling
HBTInt next;	//its next sibling
HBTInt sub;	//its sub-in-sub
} ;
struct SubProperty
{
	HBTReal CoM[3];//center of mass position (comoving)
	HBTReal VCoM[3];//CoM velocity (physical)
	HBTReal Pot;//average gavitational potential per unit mass per particle (GM/R_physical)
	HBTReal Kin;//average Kinetic energy per unit mass per particle (0.5*V_physical^2)
	HBTReal AM[3];//average angular momentum per unit mass per particle (R_physical x V_physical)
};
typedef struct
{
	HBTInt Ngroups;
	HBTInt Nsubs;
	HBTInt Nids;
	HBTInt *GrpLen_Sub;//number of subs in a group
	HBTInt *GrpOffset_Sub;//sub index offset of a group
	HBTInt *SubLen;
	HBTInt *SubOffset;

	HBTInt *SubRank;
	struct Chain_data *HaloChains;
	
	struct Hierarchy *sub_hierarchy;
	
	struct SubProperty *Property;   
	
	HBTInt **PSubArr;//PSubArr[Nsubs] as a pointer will point to PIDs block (or PIndices during unbind()) of each sub;i.e.,PSubArr[subhaloid][particle] will give PIDs(or PInd);
	
	HBTInt Nbirth;
	HBTInt NQuasi;
	HBTInt Ndeath;
	HBTInt Nsplitter;//the splitted-out subs from last snap
} SUBCATALOGUE;

#define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){ fprintf(logfile,"Error opening file '%s'\n",filename);	fflush(logfile); exit(1);	}
void *mymalloc(size_t n)
{void * mem;
	if(n)
	{
		if(!(mem=malloc(n)))
		{fprintf(logfile,"failed to allocate memory for %u bytes.\n",(unsigned) n);fflush(logfile);
		exit(1);
		}
	}
	else
	{
		mem=NULL;
	}
	return mem;
}
void myfree(void *mem)
{
	if(mem!=NULL)
	free(mem);
}

void load_sub_catalogue(HBTInt Nsnap, SUBCATALOGUE *Cat, char *SubCatPath)
{
FILE *fd;
char buf[1024];
HBTInt i;

  sprintf(buf, "%s/subcat_%03d", SubCatPath, Nsnap);
  if(!(fd = fopen(buf, "r")))
    {
	fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
	exit(1);
    }

  fread(&Cat->Ngroups, sizeof(HBTInt), 1, fd);
    Cat->GrpOffset_Sub=mymalloc(sizeof(HBTInt)*Cat->Ngroups);   
	Cat->GrpLen_Sub=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
  fread(&Cat->Nsubs, sizeof(HBTInt), 1, fd);
	Cat->SubLen=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->SubOffset=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->SubRank=mymalloc(sizeof(HBTInt)*Cat->Nsubs);
	Cat->HaloChains=mymalloc(sizeof(struct Chain_data)*Cat->Nsubs);
	Cat->sub_hierarchy=mymalloc(sizeof(struct Hierarchy)*Cat->Nsubs);
	Cat->Property=mymalloc(sizeof(struct SubProperty)*Cat->Nsubs);
	Cat->PSubArr=mymalloc(sizeof(HBTInt *)*Cat->Nsubs);
  fread(&Cat->Nids, sizeof(HBTInt), 1, fd);
  fread(Cat->GrpLen_Sub, sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->GrpOffset_Sub,sizeof(HBTInt), Cat->Ngroups, fd);
  fread(Cat->SubLen, sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubOffset,sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->SubRank,sizeof(HBTInt), Cat->Nsubs, fd);
  fread(Cat->HaloChains,sizeof(struct Chain_data), Cat->Nsubs, fd);
  fread(Cat->sub_hierarchy,sizeof(struct Hierarchy),Cat->Nsubs,fd);
  fread(Cat->Property,sizeof(struct SubProperty), Cat->Nsubs, fd);

for(i=0;i<Cat->Nsubs;i++)
{
	Cat->PSubArr[i]=mymalloc(sizeof(HBTInt)*Cat->SubLen[i]);
	fread(Cat->PSubArr[i], sizeof(HBTInt), Cat->SubLen[i], fd);
}
fread(&Cat->Nbirth,sizeof(HBTInt),1,fd);
fread(&Cat->NQuasi,sizeof(HBTInt),1,fd);
fread(&Cat->Ndeath,sizeof(HBTInt),1,fd);
fread(&Cat->Nsplitter,sizeof(HBTInt),1,fd);

  fclose(fd);
}

void load_pro2dest(HBTInt Nsnap_pro,HBTInt **P2pro2dest,HBTInt *P2Nsubs,char *destdir)
{
FILE *fp;
char buf[1024];
HBTInt Nsubs,*pro2dest;

  sprintf(buf, "%s/pro2dest/pro2dest_%03d", destdir, Nsnap_pro);
  myfopen(fp,buf,"r");
  fread(&Nsubs,sizeof(HBTInt),1,fp);//Nsubs=Npro+Nspl
  pro2dest=mymalloc(sizeof(HBTInt)*(Nsubs+1));
  pro2dest++;
  fread(pro2dest,sizeof(HBTInt),Nsubs,fp);
  fclose(fp);
  pro2dest[-1]=-1;//no pro, no dest
  *P2Nsubs=Nsubs;
  *P2pro2dest=pro2dest;
}
void free_pro2dest(HBTInt *pro2dest)
{
	free(pro2dest-1);
}

void free_sub_catalogue(SUBCATALOGUE *SubCat)//patched version
{
	myfree(SubCat->GrpOffset_Sub);
	myfree(SubCat->GrpLen_Sub);
	myfree(SubCat->SubLen);
	myfree(SubCat->SubOffset);
	myfree(SubCat->SubRank);
	myfree(SubCat->HaloChains);
	myfree(SubCat->sub_hierarchy);
	myfree(SubCat->Property);
	myfree(SubCat->PSubArr);//this only frees the pointer array, remember to reuse or free the actual pid array.
}
void erase_sub_catalogue(SUBCATALOGUE *SubCat)
{
	HBTInt i;
	for(i=0;i<SubCat->Nsubs;i++)
	{
		myfree(SubCat->PSubArr[i]);
	}
	free_sub_catalogue(SubCat);
}

void load_sp2pro(HBTInt Nsnap_dest,HBTInt *P2Npro,HBTInt *P2Nsplitter_dest, HBTInt **P2sp2pro,char *SubcatPath)
{/*	Npro=SubCatPro.Nsubs;	Nsplitter_dest=SubCatDest.Nsplitter
	* sp2pro(ProSubID) fixes ProSubID to the valid subid if ProSubID>=Npro;
	* returns NULL if Nsplitter==0
	* */
	
	char buf[1024];
	HBTInt *sp2pro,Npro,Nsplitter;
	FILE *fp;
	
	sprintf(buf, "%s/splitters/sp2pro_%03d",SubcatPath, Nsnap_dest);
	if(!(fp = fopen(buf, "r")))
	{
		fprintf(logfile,"can't open file `%s'\n", buf);fflush(logfile);
		exit(1);
	}
	fread(&Nsplitter,sizeof(HBTInt),1,fp);
	fread(&Npro,sizeof(HBTInt),1,fp);
	*P2Npro=Npro;
	*P2Nsplitter_dest=Nsplitter;
	*P2sp2pro=NULL;
	if(Nsplitter)
	{
	sp2pro=mymalloc(sizeof(HBTInt)*Nsplitter);
	fread(sp2pro,sizeof(HBTInt),Nsplitter,fp);
	*P2sp2pro=sp2pro-Npro;
	}
	fclose(fp);
}
void free_sp2pro(HBTInt *sp2pro,HBTInt Npro,HBTInt Nsplitter_dest)
{
	if(Nsplitter_dest)
	free(sp2pro+Npro);
}



