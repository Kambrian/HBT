#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

typedef struct
{
	int mass;//fof mass
	int Mvir[3];
	float Rvir[3];//[tophat,c200,b200],comoving
	int flag_badvir[3];
	int flag_fakehalo;//set to 1 when halo is not self-bound,to 0 otherwise.
} HALOSIZE;
void load_halo_size(HALOSIZE *halosize,int Ngroups,int Nsnap);
int sum_SubInSub_mass(int subid,SUBCATALOGUE *SubCat);

int main(int argc, char** argv)
{
	CATALOGUE Cat;
	SUBCATALOGUE SubCat;
	HALOSIZE *halosize;
	int Nsnap=0;
	int grpid,subid,pid;
	float *s,*v;

	logfile=stdout;//redirect BT routines' log info to standard output
	
	if(argc!=3)
	{printf("usage: %s [Nsnap],[grpid]\n",argv[0]);fflush(stdout);exit(1);}
	else
	{
		Nsnap=atoi(argv[1]);
		grpid=atoi(argv[2]);
	}
		
	load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	//~ load_particle_data(Nsnap,SNAPSHOT_DIR);
	//~ fresh_ID2Index(&Cat,-1); 	fresh_ID2Index(&SubCat,-2);	
	halosize=mymalloc(sizeof(HALOSIZE)*Cat.Ngroups);
	load_halo_size(halosize,Cat.Ngroups,Nsnap);
	
	subid=SubCat.GrpOffset_Sub[grpid];
	s=SubCat.Property[subid].CoM;
	v=SubCat.Property[subid].VCoM;
	printf("%g,%g,%g,%g,%g,%g\n",s[0],s[1],s[2],v[0],v[1],v[2]);
	printf("%g,%g\n",halosize[grpid].Rvir[0],sqrt(G*halosize[grpid].Mvir[0]*header.mass[1]/halosize[grpid].Rvir[0]/header.time));
	printf("flag:%d,%d,%d,%d\n",halosize[grpid].flag_badvir[0],halosize[grpid].flag_badvir[1],halosize[grpid].flag_badvir[2],halosize[grpid].flag_fakehalo);
	printf("SubLen %d,NextLen %d,FofLen %d, SumLen %d\n",SubCat.SubLen[subid],SubCat.SubLen[subid+1],Cat.Len[grpid],sum_SubInSub_mass(subid,&SubCat));
	printf("SubLen/FofLen %g, SubLen/SumLen %g, NextLen/SumLen %g\n",(float)(SubCat.SubLen[subid])/Cat.Len[grpid],(float)(SubCat.SubLen[subid])/sum_SubInSub_mass(subid,&SubCat),(float)SubCat.SubLen[subid+1]/sum_SubInSub_mass(subid,&SubCat));
	free_catalogue(&Cat);
	erase_sub_catalogue(&SubCat);
	return 0;
}

void load_halo_size(HALOSIZE *halosize,int Ngroups,int Nsnap)
{
	int i;
	char buf[1024];
	FILE *fp;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	load_particle_header(Nsnap,SNAPSHOT_DIR);
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,13*4L,SEEK_CUR);
		fread(&halosize[i].mass,sizeof(int),1,fp);
		fread(halosize[i].Mvir,sizeof(int),3,fp);
		fread(halosize[i].Rvir,sizeof(float),3,fp);
		fread(halosize[i].flag_badvir,sizeof(int),3,fp);
		fread(&halosize[i].flag_fakehalo,sizeof(int),1,fp);
		if(halosize[i].flag_badvir[0])
		{
			halosize[i].Rvir[0]=comoving_virial_radius(halosize[i].mass);
			halosize[i].Mvir[0]=halosize[i].mass;
		}
	}
	fclose(fp);
}
int sum_SubInSub_mass(int subid,SUBCATALOGUE *SubCat)
{
	int id,mass;
	mass=SubCat->SubLen[subid];
	id=SubCat->sub_hierarchy[subid].sub;
	while(id>=0) //add-up sub-in-sub mass
	{
	mass+=sum_SubInSub_mass(id,SubCat);
	id=SubCat->sub_hierarchy[id].next;
	}
	return mass;
}
