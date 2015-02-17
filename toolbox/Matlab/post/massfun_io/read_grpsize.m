function grpsize=read_grpsize(file)
fid=fopen(file);
disp(['reading ',file]);
nsubs=fread(fid,1,'int32');
grpsize=fread(fid,[2,nsubs],'float32');
grpsize=grpsize';
nsubs2=fread(fid,1,'int32');
if nsubs~=nsubs2
    error(['error reading ','file :',num2str(nsubs),',',num2str(nsubs2)]);
end
fclose(fid);