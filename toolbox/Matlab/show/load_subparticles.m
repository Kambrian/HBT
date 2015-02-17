function [x,n]=load_subparticles(file);
% file='/work/Projects/HBT/code/data/subparticles/6702DM/subparticle_099_1.bin';
disp(file);
fp=fopen(file);
n=fread(fp,1,'int32');
x=fread(fp,[3,n],'float32');
tmp=fread(fp,1,'int32');
if n~=tmp
    error('load error');
    disp([n,tmp]);
end
fclose(fp);
