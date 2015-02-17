/*sample io*/
#define HBTIFMT "%d"

int load_CenVCen(HBTReal Cen[][3],HBTReal VCen[][3], HBTInt Nsnap)
{
	char buf[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;
	
	if(!Cen||!VCen)
	{
		printf("error: allocate cen, vcen first \n");
		exit(1);
	}
	#if CEN_TYPE==CEN_COM
	sprintf(buf,"%s/profile/CenVCen_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
	#elif CEN_TYPE==CEN_MBD
	sprintf(buf,"%s/profile/CenVCen_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
	#elif CEN_TYPE==CEN_MPT
	sprintf(buf,"%s/profile/CenVCen_"HBTIFMT".MPT",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	for(i=0;i<Nsubs;i++)
	fread(Cen[i],sizeof(HBTReal),3,fp);
	for(i=0;i<Nsubs;i++)
	fread(VCen[i],sizeof(HBTReal),3,fp);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
	return Nsubs;
}

int load_RmaxVmax(HBTReal *rmax,HBTReal *vmax, HBTReal *rhalf, HBTInt Nsnap)
{
	char buf[1024];
	FILE *fp;
	HBTInt Nsubs,dummy;
	//~ HBTReal *rsig, *r2sig, *r3sig, *rpoisson; /* declare this as input var if you want them!!! */
	
	if(!rmax||!vmax||!rhalf)
	{
		printf("error: allocate rmax , vmax and rhalf first \n");
		exit(1);
	}
	#if CEN_TYPE==CEN_COM
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".COM",SUBCAT_DIR,Nsnap);
	#elif CEN_TYPE==CEN_MBD
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".MBD",SUBCAT_DIR,Nsnap);
	#elif CEN_TYPE==CEN_MPT
	sprintf(buf,"%s/profile/RmaxVmax_"HBTIFMT".MPT",SUBCAT_DIR,Nsnap);
	#endif
	myfopen(fp,buf,"r");	
	fread(&Nsubs,sizeof(HBTInt),1,fp);
	fread(rmax,sizeof(HBTReal),Nsubs,fp);
	fread(vmax,sizeof(HBTReal),Nsubs,fp);
	fread(rhalf,sizeof(HBTReal),Nsubs,fp);
	fseek(fp,sizeof(HBTReal)*Nsubs*4,SEEK_CUR);
	//~ fread(rsig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r2sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(r3sig,sizeof(HBTReal),Nsubs,fp);
	//~ fread(rpoisson,sizeof(HBTReal),Nsubs,fp);
	fread(&dummy,sizeof(HBTInt),1,fp);
	if(dummy!=Nsubs) 
	{
		printf("error loading %s:\n size not consistent "HBTIFMT"!="HBTIFMT"\n",buf,Nsubs,dummy);
		exit(1);
	}
	fclose(fp);
	return Nsubs;
}
