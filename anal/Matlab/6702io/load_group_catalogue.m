function cat=load_group_catalogue(snapnum)
% cat=load_group_catalogue(snapnum)

global fofdir
% typedef struct 
% {
% int Ngroups;
% int Nids;
% int *Len;
% int *Offset;
% //~ float *HaloCen[3];//HaloCen[3][Ngroups], center of halos
% int *PIDorIndex; //stores PID when loading FOF then changes to PIndex after loading particles
% short *HaloMask; //HaloMask[NP_DM],HaloMask==1 means the particle does not belong to any sub,i.e, it's free.
% short *HaloMaskSrc;
% int *ID2Halo;//for index2halo, this is ID2halo[NP_DM];
% }CATALOGUE;

cat=struct('Ngroups',0,'Nids',0,'Len',[],'Offset',[],'PIDorIndex',[],'property',struct('type','cat','particles','id','snap',-1));
grptab=fullfile(fofdir,['group_tab_',num2str(snapnum,'%03d')]);
fid=fopen(grptab,'r');
cat.Ngroups=fread(fid,1,'int32');
cat.Nids=fread(fid,1,'int32');
TotNgroups=fread(fid,1,'int32'); %#ok<NASGU>
Nfiles=fread(fid,1,'int32'); %#ok<NASGU>
cat.Len=fread(fid,cat.Ngroups,'int32');
cat.Offset=fread(fid,cat.Ngroups,'int32');
fclose(fid);

grpids=fullfile(fofdir,['group_ids_',num2str(snapnum,'%03d')]);
fid=fopen(grpids,'r');
fseek(fid,4*4,'bof');
cat.PIDorIndex=fread(fid,cat.Nids,'int32');
fclose(fid);

cat.property.snap=snapnum;