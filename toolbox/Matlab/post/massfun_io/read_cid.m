function cenID=read_cid(file)
fid=fopen(file);
disp(['reading ',file]);
nsubs=fread(fid,1,'int32');
cenID=fread(fid,nsubs,'int32');
nsubs2=fread(fid,1,'int32');
if nsubs~=nsubs2
    error(['error reading ','file :',num2str(nsubs),',',num2str(nsubs2)]);
end
fclose(fid);