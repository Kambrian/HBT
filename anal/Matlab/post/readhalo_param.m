function par=readhalo_param(basedir,Nsnap)

sizefile=fullfile(basedir,['halo_param_',num2str(Nsnap,'%03d')])

fid=fopen(sizefile,'r');
ngroups=fread(fid,1,'int32');

par.stat=zeros(ngroups,1);
par.rhos=zeros(ngroups,2);
par.rs=zeros(ngroups,2);
par.c=zeros(ngroups,2);
par.Merr=zeros(ngroups,1);

for i=1:ngroups
par.stat(i)=fread(fid,1,'int32');
par.rhos(i,:)=fread(fid,[1,2],'float32');
par.rs(i,:)=fread(fid,[1,2],'float32');
par.c(i,:)=fread(fid,[1,2],'float32');
par.Merr(i)=fread(fid,1,'float32');
end

ngroups2=fread(fid,1,'int32');
fclose(fid);
if(ngroups~=ngroups2)
    error(['reading error: ngroups=',ngroups,',or ',ngroups2]);
end
