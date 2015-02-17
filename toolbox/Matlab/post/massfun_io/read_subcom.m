function subcom=read_subcom(file)
fid=fopen(file);
disp(['reading ',file]);
nsubs=fread(fid,1,'int32');
subcom=fread(fid,[3,nsubs],'float32');
subcom=subcom';
nsubs2=fread(fid,1,'int32');
if nsubs~=nsubs2
    error(['error reading ','file :',num2str(nsubs),',',num2str(nsubs2)]);
end
fclose(fid);