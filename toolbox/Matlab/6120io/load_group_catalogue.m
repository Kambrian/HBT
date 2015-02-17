function cat=load_group_catalogue(snapnum)
% cat=load_group_catalogue(snapnum)

global fofdir header RUN_NUM snaplist
% typedef struct 
% {
% int Ngroups;
% int Nids;
% int *Len;
% int *Offset;
% float *HaloCen[3];//HaloCen[3][Ngroups], center of halos
% int *PIDorIndex; //stores PID when loading FOF then changes to PIndex after loading particles
% short *HaloMask; //HaloMask[NP_DM],HaloMask==1 means the particle does not belong to any sub,i.e, it's free.
% short *HaloMaskSrc;
% int *ID2Halo;//for index2halo, this is ID2halo[NP_DM];
% }CATALOGUE;

cat=struct('Ngroups',0,'Nids',0,'Len',[],'Offset',[],'HaloCen',[],'PIDorIndex',[],'property',struct('type','cat','particles','matlab-ind','snap',-1));

grpfile=fullfile(fofdir,['fof.b20.',num2str(RUN_NUM),'.',num2str(snaplist(snapnum+1),'%04d')]);
fid=fopen(grpfile,'r','ieee-be');

dummy=fread(fid,1,'int32');
b=fread(fid,1,'float32') %#ok<NOPRT,NASGU>
cat.Ngroups=fread(fid,1,'int32');
dummy2=fread(fid,1,'int32');
if dummy~=dummy2
    error(['error reading group header:',grpfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
end

cat.HaloCen=zeros(cat.Ngroups,3);
for i=1:3
  dummy=fread(fid,1,'int32');
  cat.HaloCen(:,i)=fread(fid,cat.Ngroups,'float32');
  dummy2=fread(fid,1,'int32');
  if dummy~=dummy2
      error(['error reading group center:',grpfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
  end
end
cat.HaloCen=cat.HaloCen*header.rLbox*1000;

cat.Nids=0;
% cat.PIDorIndex=zeros(header.Np,1);
for i=1:cat.Ngroups
    dummy=fread(fid,1,'int32');
    cat.Len(i)=fread(fid,1,'int32');
    dummy2=fread(fid,1,'int32');
    if dummy~=dummy2
        error(['error reading group Len:',grpfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
    end
    cat.Offset(i)=cat.Nids;
    cat.Nids=cat.Nids+cat.Len(i);
    
    dummy=fread(fid,1,'int32');
%     cat.PIDorIndex(cat.Offset(i)+(1:cat.Len(i)))=fread(fid,cat.Len(i),'int32');
fseek(fid,cat.Len(i)*4,'cof');
    dummy2=fread(fid,1,'int32');
    if dummy~=dummy2
        error(['error reading group ids:',grpfile,';\nbrackets mismatch:',num2str(dummy),',',num2str(dummy2)]);
    end
end
% cat.PIDorIndex(cat.Nids+1:end)=[];

cat.property.snap=snapnum;