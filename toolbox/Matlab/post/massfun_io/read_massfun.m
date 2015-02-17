function [data,redshift]=read_massfun(file)
fid=fopen(file);
redshift=fread(fid,1,'float32');
nfun=fread(fid,1,'int32');
virtype=fread(fid,1,'int32');virtype
rmin=fread(fid,1,'float32');rmin
rmax=fread(fid,1,'float32');rmax
data=struct([]);
for i=1:nfun
nbin=fread(fid,1,'int32');
data(i).nbin=nbin;
data(i).Mbin=fread(fid,2,'float32');
data(i).Nhost=fread(fid,1,'int32');
data(i).Mhost=fread(fid,1,'float64');
data(i).xmass=fread(fid,[2,nbin],'float32')';
data(i).mfunspec=fread(fid,[2,nbin],'float32')';
data(i).mfunspecln=fread(fid,[2,nbin],'float32')';
data(i).mfuncum=fread(fid,[2,nbin],'float32')';
end
nfun2=fread(fid,1,'int32');
if nfun~=nfun2
    error(['error reading ',file,': check failed']);
end
fclose(fid);