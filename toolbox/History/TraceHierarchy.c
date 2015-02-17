//print hierarchy history of given subhalos
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"

#define LEN(x) ((x)<0?0:SubCat.SubLen[(x)])

struct sub_info
{
	int subid;
	int sublen;
	int subrank;
	struct Hierarchy hierarchy;
	struct Chain_data chain;
};

void read_subinfo(HBTInt Nsnap,struct sub_info * s)
{
	FILE *fd;
	char buf[1024];
	HBTInt Ngroups,Nsubs,Nids;

	  sprintf(buf, "%s/subcat_%03d", SUBCAT_DIR, (int)Nsnap);
	    myfopen(fd,buf,"r");
	    fread(&Ngroups,sizeof(HBTInt),1,fd);
	    fread(&Nsubs,sizeof(HBTInt),1,fd);
	    fread(&Nids,sizeof(HBTInt),1,fd);

	    fseek(fd,sizeof(HBTInt)*(Ngroups*2),SEEK_CUR);
	    fseek(fd,sizeof(HBTInt)*s->subid,SEEK_CUR);
	    fread(&s->sublen,sizeof(HBTInt),1,fd);
	    fseek(fd,sizeof(HBTInt)*(2*Nsubs-1),SEEK_CUR);
	    fread(&s->subrank,sizeof(HBTInt),1,fd);
	    fseek(fd,sizeof(HBTInt)*(Nsubs-s->subid-1),SEEK_CUR);

	    fseek(fd,sizeof(struct Chain_data)*s->subid,SEEK_CUR);
	    fread(&s->chain,sizeof(struct Chain_data),1,fd);
	    fseek(fd,sizeof(struct Chain_data)*(Nsubs-s->subid-1),SEEK_CUR);
	    
	    fseek(fd,sizeof(struct Hierarchy)*s->subid,SEEK_CUR);
	    fread(&s->hierarchy,sizeof(struct Hierarchy),1,fd);

	    fclose(fd);
}

int main(int argc,char **argv)
{
	SUBCATALOGUE SubCat;
	CATALOGUE Cat;
	HBTInt Nsnap,i;
	struct sub_info s[3];
	
	logfile=stdout;
	
	s[0].subid=334453;
	s[1].subid=334457;
	s[2].subid=334459;
	
	FILE *fp;
	char buf[1024];
	sprintf(buf,"%s/anal/check_hierarchy.dat",SUBCAT_DIR);
	myfopen(fp,buf,"w");
	
	for(Nsnap=99;Nsnap>0;Nsnap--)
	{
	fprintf(fp,"Nsnap=%d\n",Nsnap);
	//load_sub_table(Nsnap,&SubCat,SUBCAT_DIR);
	for(i=0;i<3;i++)
	{
		if(s[i].subid>0)
		{
			read_subinfo(Nsnap,&s[i]);
			fprintf(fp,"%d: %d,%d,%d\t(%d,%d,%d,%d)\t(%d,%d)\n",
					i,s[i].subid,s[i].sublen,s[i].subrank,s[i].hierarchy.nibs,s[i].hierarchy.pre,s[i].hierarchy.next,s[i].hierarchy.sub,s[i].chain.ProSubID,s[i].chain.HostID);
			s[i].subid=s[i].chain.ProSubID;
		}
	}
	fprintf(fp,"\n");
	if(s[0].subid<0&&s[1].subid<0&&s[2].subid<0) break;
	}
	
	fclose(fp);
	
	//free_sub_table(&SubCat);
	return 0;
}
