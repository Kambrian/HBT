function [Snappar,Mratepar,Mhostpar,Rhostpar,Kpar,j2par,Kerrpar,j2errpar,Chostpar,Csatpar]=load_hist_par(file)
fid=fopen(file);
tmp=fread(fid,1,'int32');
Nhist=tmp;
Snappar=zeros(Nhist,3);
Mratepar=zeros(Nhist,2);
Mhostpar=zeros(Nhist,2);
Rhostpar=zeros(Nhist,2);
Kpar=zeros(Nhist,2);
j2par=zeros(Nhist,2);
Kerrpar=zeros(Nhist,2);
j2errpar=zeros(Nhist,2);
Chostpar=zeros(Nhist,2);
Csatpar=zeros(Nhist,2);
for i=1:Nhist
Snappar(i,:)=fread(fid,[1,3],'int32');
Mratepar(i,:)=fread(fid,[1,2],'float32');
Mhostpar(i,:)=fread(fid,[1,2],'float32');
Rhostpar(i,:)=fread(fid,[1,2],'float32');
Kpar(i,:)=fread(fid,[1,2],'float32');
j2par(i,:)=fread(fid,[1,2],'float32');
Kerrpar(i,:)=fread(fid,[1,2],'float32');
j2errpar(i,:)=fread(fid,[1,2],'float32');
Chostpar(i,:)=fread(fid,[1,2],'float32');
Csatpar(i,:)=fread(fid,[1,2],'float32');
end
tmp2=fread(fid,1,'int32');
if tmp2~=tmp
    error('Nhsit');
end
fclose(fid);
