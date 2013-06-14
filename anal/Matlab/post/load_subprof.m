function prof=load_subprof(proffile,nbin)
nsub=numel(nbin);
prof.ns=cell(nsub,1);
prof.no=cell(nsub,1);
prof.nb=cell(nsub,1);
prof.rs=cell(nsub,1);
prof.ro=cell(nsub,1);
prof.rb=cell(nsub,1);

fid=fopen(proffile,'r');
for i=1:nsub
prof.ns{i}=fread(fid,nbin(i),'int32');
prof.no{i}=fread(fid,nbin(i),'int32');
prof.nb{i}=fread(fid,nbin(i),'int32');
prof.rs{i}=fread(fid,nbin(i),'float32');
prof.ro{i}=fread(fid,nbin(i),'float32');
prof.rb{i}=fread(fid,nbin(i),'float32');
end
fclose(fid);