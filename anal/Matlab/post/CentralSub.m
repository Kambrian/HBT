file='/work/Projects/HBT/code/data/CentralSub_059';
fp=fopen(file);
n=fread(fp,1,'int32');
DCen=fread(fp,n,'float32');
DMin=fread(fp,n,'float32');
Mvir=fread(fp,n,'float32');
Rvir=fread(fp,n,'float32');
RankMin=fread(fp,n,'int32');
n2=fread(fp,1,'int32');
if n~=n2
    error([num2str(n),',',num2str(n2)]);
end

f=DCen>60000;
DCen(f)=[];
DMin(f)=[];
Mvir(f)=[];
Rvir(f)=[];
RankMin(f)=[];

f=~Mvir;
DCen(f)=[];
DMin(f)=[];
Mvir(f)=[];
Rvir(f)=[];
RankMin(f)=[];
%%
f=Mvir<1000*min(Mvir);
DCen(f)=[];
DMin(f)=[];
Mvir(f)=[];
Rvir(f)=[];
RankMin(f)=[];

