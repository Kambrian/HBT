%% for the data given to Yang

clear;
runnum='6113';

file=['/mnt/A4700/data/',runnum,'/subcat/anal/infalldata'];
fid=fopen(file,'r');
Mgrp=fread(fid,400,'float32');
GrpLen=fread(fid,400,'int32');
n=sum(GrpLen);
Mhost=fread(fid,n,'float32');
Msub=fread(fid,n,'float32');
tinf=fread(fid,n,'float32');
Mmax=fread(fid,n,'float32');
tmax=fread(fid,n,'float32');
MRvir=fread(fid,n,'float32');
tRvir=fread(fid,n,'float32');
DirectInfall=fread(fid,n,'int32');
fclose(fid);
%%
n1=sum(GrpLen(1:100));
n=n1+sum(GrpLen(101:200));
M=Msub(n1+1:n);
% t=tmax(1:n);
xmin=min(M(M>0));xmax=max(M);
nbin=10;
x=logspace(log10(xmin),log10(xmax),nbin);
% subinR=sub_mass_099(find(Rsub>rmin*Rvir&Rsub<rmax*Rvir),1)/M0;
y=histc(M,x);
% figure;
loglog(x,y,'g*-');
hold on;
% set(gca,'xlim',[1e-2,1e4]);

%%
file=['/mnt/A4700/data/',runnum,'/subcat/anal/mainsample'];
fid=fopen(file,'r');
Mgrp=fread(fid,400,'float32');
GrpLen=fread(fid,400,'int32');
n=sum(GrpLen);
Mhost=fread(fid,n,'float32');
Msub=fread(fid,n,'float32');
tinf=fread(fid,n,'float32');
fclose(fid);
