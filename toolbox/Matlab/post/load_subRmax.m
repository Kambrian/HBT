function s=load_subRmax(sizefile)
fid=fopen(sizefile,'r');
% int Nsub
% 	float rmax,vmax,rhalf,r1sig,r2sig,r3sig,rpoisson
% int Nsub
nsub=fread(fid,1,'int32');
s.nsub=nsub;
s.rmax=fread(fid,nsub,'float32');
s.vmax=fread(fid,nsub,'float32');
s.rhalf=fread(fid,nsub,'float32');
s.r1sig=fread(fid,nsub,'float32');
s.r2sig=fread(fid,nsub,'float32');
s.r3sig=fread(fid,nsub,'float32');
s.rpoisson=fread(fid,nsub,'float32');
n=fread(fid,1,'int32');
if n~=nsub
    error(['error reading: ',sizefile, num2str(n), '!=', numstr(nsub)]);
end

