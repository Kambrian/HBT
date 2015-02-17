function M=load_deathpro(basedir,SnapDeath,grpid)
file=[basedir,'steller/DeathInfall_first_',num2str(SnapDeath,'%03d'),'_',num2str(grpid)];
disp(['reading ',file]);
fid=fopen(file,'r');
N=fread(fid,1,'int32');
M=fread(fid,[N,1],'float32');
N2=fread(fid,1,'int32');
fclose(fid);
if N~=N2
    error('error reading file,check fail')
end
