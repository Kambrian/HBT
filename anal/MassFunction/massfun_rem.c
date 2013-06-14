//remnant subhalo mass function of the subhalos from SnapLoad
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"
//~ #include "gas_vars.h"
//~ #include "gas_proto.h"

void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);
int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat;
	int *pro2dest;
	FILE *fp;
	char buf[1024];
	int Nsnap,SnapLoad;
	int i,j,pid,grpid,subid, *sublist;
	int Ngroups,cid,Npro,Nsubs;
	double partmass;
	float mass,(*Mvir)[3],(*Rvir)[3];
	char outputdir[1024];
	
	logfile=stdout;
	if(argc!=3)
	{
	printf("usage:%s [Snap] [grpid]\n",argv[0]);
	exit(1);
	}
	SnapLoad=atoi(argv[1]);
	grpid=atoi(argv[2]);

	load_sub_catalogue(SnapLoad,&SubCat,SUBCAT_DIR);
	load_particle_header(SnapLoad,SNAPSHOT_DIR);
	partmass=header.mass[1];
	printf("dm %g,baryon %g, fraction %g,Nsplitter %d, Nquasi %d\n",
		partmass,header.mass[0],header.mass[0]/header.mass[1],SubCat.Nsplitter,SubCat.NQuasi);
	Nsubs=SubCat.GrpLen_Sub[grpid];
	sublist=mymalloc(sizeof(int)*Nsubs);
	for(i=0;i<Nsubs;i++)
		sublist[i]=i+SubCat.GrpOffset_Sub[grpid];
	erase_sub_catalogue(&SubCat);

	sprintf(outputdir,"%s/anal/massfun/rem_%03d",SUBCAT_DIR,SnapLoad);	
	mkdir(outputdir,0755);		
	for(Nsnap=SnapLoad+1;Nsnap<MaxSnap;Nsnap++)
	{
		load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
		if(sublist[0]<0) break; 
		grpid=SubCat.HaloChains[sublist[0]].HostID;
		if(grpid<0) break;
		if(SubCat.SubRank[sublist[0]]!=0) break; //no longer a halo
			
		load_particle_header(Nsnap,SNAPSHOT_DIR);
		partmass=header.mass[1];
		printf("dm %g,baryon %g, fraction %g,Nsplitter %d, Nquasi %d\n",
			partmass,header.mass[0],header.mass[0]/header.mass[1],SubCat.Nsplitter,SubCat.NQuasi);
					
		sprintf(buf,"%s/submass_%03d",outputdir,Nsnap);
		myfopen(fp,buf,"w");
		fwrite(&Nsubs,sizeof(int),1,fp);			
		for(i=0;i<Nsubs;i++)
		{
			subid=sublist[i];
			mass=SubCat.SubLen[subid]*partmass;//+GSubCat.SubLen[i]*header.mass[0];
			fwrite(&mass,sizeof(float),1,fp);
		}
		fwrite(&Nsubs,sizeof(int),1,fp);		
		fclose(fp);
	
		sprintf(buf,"%s/subcom_%03d",outputdir,Nsnap);
		myfopen(fp,buf,"w");
		fwrite(&Nsubs,sizeof(int),1,fp);		
		for(i=0;i<Nsubs;i++)
		{
			subid=sublist[i];
			fwrite(SubCat.Property[subid].CoM,sizeof(float),3,fp);
		}
		fwrite(&Nsubs,sizeof(int),1,fp);		
		fclose(fp);
			
		sprintf(buf,"%s/cid_%03d",outputdir,Nsnap);
		myfopen(fp,buf,"w");
		Ngroups=1;
		fwrite(&Ngroups,sizeof(int),1,fp);
		cid=0;
		fwrite(&cid,sizeof(int),1,fp);
		Ngroups=1;
		fwrite(&Ngroups,sizeof(int),1,fp);
		fclose(fp);
	
		Mvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
		Rvir=mymalloc(sizeof(float)*3*SubCat.Ngroups);
		load_halo_virial_size(Mvir,Rvir,(float)partmass,SubCat.Ngroups,Nsnap);
		sprintf(buf,"%s/grpsizeVIR_%03d",outputdir,Nsnap);
		myfopen(fp,buf,"w");
		fwrite(&Ngroups,sizeof(int),1,fp);
		fwrite(&Mvir[grpid][0],sizeof(float),1,fp);
		fwrite(&Rvir[grpid][0],sizeof(float),1,fp);
		fwrite(&Ngroups,sizeof(int),1,fp);
		fclose(fp);
		sprintf(buf,"%s/grpsizeC200_%03d",outputdir,Nsnap);
		myfopen(fp,buf,"w");
		fwrite(&Ngroups,sizeof(int),1,fp);
		fwrite(&Mvir[grpid][1],sizeof(float),1,fp);
		fwrite(&Rvir[grpid][1],sizeof(float),1,fp);
		fwrite(&Ngroups,sizeof(int),1,fp);
		fclose(fp);
		sprintf(buf,"%s/grpsizeB200_%03d",outputdir,Nsnap);
		myfopen(fp,buf,"w");
		fwrite(&Ngroups,sizeof(int),1,fp);
		fwrite(&Mvir[grpid][2],sizeof(float),1,fp);
		fwrite(&Rvir[grpid][2],sizeof(float),1,fp);
		fwrite(&Ngroups,sizeof(int),1,fp);
		fclose(fp);
		myfree(Mvir);
		myfree(Rvir);	
		erase_sub_catalogue(&SubCat);
		
		if(Nsnap<MaxSnap-1)
		{
			int Npro,Nalive;
			load_pro2dest(Nsnap,&pro2dest,&Npro,SUBCAT_DIR);
			Nalive=0;
			for(i=0;i<Nsubs;i++)
			{
				subid=pro2dest[sublist[i]];
				if(subid>=0)
				{
				sublist[Nalive]=subid;
				Nalive++;
				}			
			}
			Nsubs=Nalive;
			free_pro2dest(pro2dest);
		}	
	}							
	return 0;
}

void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
{
	char buf[1024];
	FILE *fp;
	int Nvir[3],i,j;
	sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
	myfopen(fp,buf,"r");
	for(i=0;i<Ngroups;i++)
	{
		fseek(fp,14*4L,SEEK_CUR);
		fread(Nvir,sizeof(int),3,fp);
		for(j=0;j<3;j++)
		Mvir[i][j]=Nvir[j]*partmass;
		fread(Rvir+i,sizeof(float),3,fp);
		fseek(fp,4*4L,SEEK_CUR);
	}
	fclose(fp);
}
