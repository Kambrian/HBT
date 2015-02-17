function mass=read_biggestsub(file)
fid=fopen(file);
disp(['reading ',file]);
ngrps=fread(fid,1,'int32');
mass=fread(fid,[3,ngrps],'float32');
mass=mass';
ngrps2=fread(fid,1,'int32');
if ngrps~=ngrps2
    error(['error reading ','file :',num2str(ngrps),',',num2str(ngrps2)]);
end
fclose(fid);