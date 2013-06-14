#include <stdlib.h>
#include <stdio.h>

int count_lines(char * filename);
int main()
{
int i;
FILE *fp;
char ch='a';
fp=fopen("testread.dat","r");//an empty file
putchar(ch);
printf("%d\n",ch);
if(feof(fp)) printf("EoF\n");
ch=fgetc(fp);//get -1 when eof is encountered
putchar(ch);
printf("%d\n",ch);
fclose(fp);

fp=fopen("testlines.dat","w");
fprintf(fp,"1\n2\n3");
fclose(fp);
i=count_lines("testlines.dat");
printf("%d lines\n",i);
return 0;
}

int count_lines(char * filename)
{// return number of lines
	FILE *fp;
	if(!(fp=fopen(filename,"r")))
	{
		 return -1;
	 }
	 int i=0;
	 char ch;
	 while(1)
	 {
		 ch=fgetc(fp);
		 if(feof(fp))
		 {
			fseek(fp,-sizeof(char),SEEK_END);
			ch=fgetc(fp);
			if(ch!='\n') i++; //in case the last line is not terminated with '\n'
			break;
		 }
		 if(ch=='\n')
			i++;
	 }
	 fclose(fp);
	return i;
}
