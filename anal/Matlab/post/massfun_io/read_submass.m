function submass=read_submass(file)
fid=fopen(file);
disp(['reading ',file]);
nsubs=fread(fid,1,'int32');
submass=fread(fid,nsubs,'float32');
nsubs2=fread(fid,1,'int32');
if nsubs~=nsubs2
    error(['error reading ','file :',num2str(nsubs),',',num2str(nsubs2)]);
end
fclose(fid);