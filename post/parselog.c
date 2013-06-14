#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

int get_new_word(char *s,FILE *fp,int strlen_max)
{
	int n;
	n=0;
	while(!(isalpha(*s=fgetc(fp))))
	{
		if(feof(fp))
		return -1;
	}
	n++;
	s++;
	while(isalpha(*s=fgetc(fp)))
	{
		s++;
		n++;
		if(n==strlen_max-1)
		{
			*s='\0';
			return n;
		}
	}
	*s='\0';
	return n;
}
int seek_new_line(FILE *fp)
{
	char tmp;
	tmp=fgetc(fp);
	while(tmp!='\n')
	{
		if(feof(fp))
			return -1;
		tmp=fgetc(fp);
	}
	return 1;
}
int main(int argc, char** argv)
{
	FILE * fp, *logfile;
		char inputdir[512]="/SANdisk5/kambrain/Sim6702/SubCat6";
		char outputdir[1024]="/SANdisk5/kambrain/Sim6702/SubCat6/anal";
		char buf[1024];
		char fofstr[20]="Snap";
		char substr[20]="Nsublast";
		char failstr[20]="Nff";
		char wordtmp[20];
		int Nsnap,Ngroups,Nsubs,Nbirth,Ndeath,Nquasi,Nsplitter,Nff,Nfm;
		int junk;
		
		sprintf(buf,"%s/substat",outputdir);
	if(!(fp=fopen(buf,"w")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
			sprintf(buf,"%s/logfile",inputdir);
	if(!(logfile=fopen(buf,"r")))
	{
		printf("Error opening file '%s'\n",buf);
		exit(1);
	}
	while((get_new_word(wordtmp,logfile,20))>0)
	{
		if(!strcmp(wordtmp,fofstr))
			{
				fscanf(logfile,"%d Ngroups=%d",&Nsnap,&Ngroups);
				seek_new_line(logfile);
			}
		else if(!strcmp(wordtmp,substr))
			fscanf(logfile,"%d\tNsub=%d\tBirth=%d\tDeath=%d\tQuasi=%d\tSplitter=%d\n",&junk,&Nsubs,&Nbirth,&Ndeath,&Nquasi,&Nsplitter);
		else if(!strcmp(wordtmp,failstr))
		{
			fscanf(logfile,"%d\tNfm=%d\n",&Nff,&Nfm);
			fprintf(fp,"%d,%d,%d,%d,%d,%d,%d,%d,%d\n",Nsnap,Ngroups,Nsubs,Nbirth,Ndeath,Nquasi,Nsplitter,Nff,Nfm);
		}
		else
			seek_new_line(logfile);
	}
	fclose(fp);
	fclose(logfile);
	
	return 0;
}


	
