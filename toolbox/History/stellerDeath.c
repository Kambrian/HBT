//progenitor mass function for those subhalos which are dead before a give snapshot
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

typedef struct 
{
	int *mass;
	int len;
} MassList;
int trace_death_progenitors(MassList *Mlist,int Nsnap,int grpid)
{
	/* output: Mlist 
	 * return value:  haloid at Nsnap-1
	 * input: Nsnap, snapshot just before death
	 *        grpid, haloid at Nsnap		  
	 * */
	SUBCATALOGUE SubCat;
	//~ CATALOGUE Cat;
	int i,j,subid,Nsub,flag_continue;
	int *sp2pro,*pro2dest;
	int Npro,Nsplitter,hgrpid,hsubid;
	int *sublist,*ranklist;

	load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
	//~ load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
	load_pro2dest(Nsnap,&pro2dest,&Nsub,SUBCAT_DIR);
	load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);
	hsubid=SubCat.GrpOffset_Sub[grpid];
	Nsub=SubCat.GrpLen_Sub[grpid];
	sublist=mymalloc(sizeof(int)*Nsub);
	ranklist=mymalloc(sizeof(int)*Nsub);
	Mlist->mass=mymalloc(sizeof(int)*Nsub);
	j=0;
	for(i=0;i<Nsub;i++)
	{
		subid=hsubid+i;
		if(pro2dest[subid]<0)
		{
			sublist[j]=SubCat.HaloChains[subid].ProSubID;
			if(sublist[j]>=Npro)
				sublist[j]=sp2pro[sublist[j]];
			ranklist[j]=SubCat.SubRank[subid];
			if(!(SubCat.SubRank[subid]))// a central 
			//~ Mlist->mass[j]=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]):(Cat.Len[SubCat.HaloChains[subid].HostID]));
				Mlist->mass[j]=SubCat.SubLen[subid];
			j++;
		}
	}
	Nsub=j;
	free_sp2pro(sp2pro,Npro,Nsplitter);
	free_pro2dest(pro2dest);
	hsubid=SubCat.HaloChains[hsubid].ProSubID;
	erase_sub_catalogue(&SubCat);
	//~ free_catalogue(&Cat);
	
	if(hsubid<0)
	hgrpid=-1;
	else
	{
	load_sub_catalogue(Nsnap-1,&SubCat,SUBCAT_DIR);
	if(SubCat.SubRank[hsubid])
	hgrpid=-1;
	else
	hgrpid=SubCat.HaloChains[hsubid].HostID;//host grpid at Nsnap-1;
	erase_sub_catalogue(&SubCat);
	}
	
	Mlist->mass=realloc(Mlist->mass,sizeof(int)*Nsub);
	Mlist->len=Nsub;
	
	flag_continue=Nsub;
	while(flag_continue)
	{
		Nsnap--;
		if(Nsnap<0) break;
		flag_continue=0;
		load_sub_catalogue(Nsnap,&SubCat,SUBCAT_DIR);
		//~ load_group_catalogue(Nsnap,&Cat,GRPCAT_DIR);
		load_sp2pro(Nsnap,&Npro,&Nsplitter,&sp2pro, SUBCAT_DIR);
		for(i=0;i<Nsub;i++)
		{				
			if(sublist[i]>=0)//have a progenitor
			{
			subid=sublist[i];
			sublist[i]=SubCat.HaloChains[subid].ProSubID;
			if(sublist[i]>=Npro)
				sublist[i]=sp2pro[sublist[i]];
			if(sublist[i]>=0) flag_continue=1;//have something to go on
			if(ranklist[i])//a satellite in last snapshot
			{
				if(!(ranklist[i]=SubCat.SubRank[subid])) //ok, becomes a central at this snap
					//~ Mlist->mass[i]=((SubCat.HaloChains[subid].HostID<0)?(SubCat.SubLen[subid]):(Cat.Len[SubCat.HaloChains[subid].HostID]));
					Mlist->mass[i]=SubCat.SubLen[subid];
			}
			}
		}
		free_sp2pro(sp2pro,Npro,Nsplitter);
		//~ free_catalogue(&Cat);
		erase_sub_catalogue(&SubCat);
	}
	myfree(sublist);
	myfree(ranklist);
	
	return hgrpid;	
} 
int main(int argc,char **argv)
{
	char outputdir[1024];
	FILE *fp;
	char buf[1024];
	
	SUBCATALOGUE SubCat;
	MassList Mlist;
	int maxlen;
	int i,grpid,hsubid,SnapDeath,Nsnap;
	
	sprintf(outputdir,"%s/anal/steller",SUBCAT_DIR);	
	mkdir(outputdir,0755);		
	logfile=stdout;
	
	if(argc!=3)
	{
		printf(" %s [SnapDeath] [grpid]\n",argv[0]);
		exit(1);
	}
	SnapDeath=atoi(argv[1]);
	grpid=atoi(argv[2]);
	if(SnapDeath<1) return 0;
	
	sprintf(buf, "%s/DeathInfall_first_%03d_%d",outputdir,SnapDeath,grpid);
	myfopen(fp,buf,"w");
	
	load_sub_catalogue(SnapDeath,&SubCat,SUBCAT_DIR);
	hsubid=SubCat.GrpOffset_Sub[grpid];
	hsubid=SubCat.HaloChains[hsubid].ProSubID;
	erase_sub_catalogue(&SubCat);
	load_sub_catalogue(SnapDeath-1,&SubCat,SUBCAT_DIR);
	grpid=SubCat.HaloChains[hsubid].HostID;
	maxlen=SubCat.GrpLen_Sub[grpid]*2;
	erase_sub_catalogue(&SubCat);
	
	Mlist.len=0;
	Mlist.mass=mymalloc(sizeof(int)*maxlen);
	for(Nsnap=SnapDeath-1;Nsnap>0;Nsnap--)
	{
		MassList Slist;
		grpid=trace_death_progenitors(&Slist,Nsnap,grpid);
		if(Slist.len+Mlist.len>maxlen)
		{
			maxlen*=2;
			Mlist.mass=realloc(Mlist.mass,sizeof(int)*maxlen);
		}
		memcpy(Mlist.mass+Mlist.len,Slist.mass,sizeof(int)*Slist.len);
		Mlist.len+=Slist.len;
		myfree(Slist.mass);
		if(grpid<0) break;
	}
	
	
	fwrite(&Mlist.len,sizeof(int),1,fp);
	
	load_particle_header(0,SNAPSHOT_DIR);
	float m;
	for(i=0;i<Mlist.len;i++)
	{
	m=Mlist.mass[i]*header.mass[1];
	fwrite(&m,sizeof(float),1,fp);
	}
	
	fwrite(&Mlist.len,sizeof(int),1,fp);
	
	return 0;
}
