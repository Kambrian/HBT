function srccat=load_src_catalogue(snapnum)
% srccat=load_src_catalogue(snapnum)
global subcatdir
% typedef struct
% {
% 	int Nsubs;
% 	int Nids;
% 	int *SubLen;
% 	int *SubLen2;
% 	float *CoreFrac;
% 	int *SubOffset;//? necessary?
% 	int **PSubArr;//PSubArr[Nsubs] as a pointer will point to PIDs block (or PIndices during unbind()) of each sub;i.e.,PSubArr[subhaloid][particle] will give PIDs(or PInd);
% 	int **PSubArr2;
% 	struct Chain_data *HaloChains;//this is only create during splitt_srccat() and then passed away to subcat during make_srcsub(); at any other times it remains NULL;
% 																			//add it here only to make it act like a global var; 
% 	int NDeathSp;//splitted to death
% }SRCCATALOGUE;

srccat=struct('Nsubs',0,'Nids',0,'SubLen',[],'SubLen2',[],'CoreFrac',[],'SubOffset',[],...
               'PSubArr',[],'PSubArr2',[],'NDeathSp',0,'property',struct('type','srccat','particles','id','snap',-1));

subfile=fullfile(subcatdir,['srccat_',num2str(snapnum,'%03d')]);
fid=fopen(subfile,'r');        
srccat.Nsubs=fread(fid,1,'int32');
srccat.Nids=fread(fid,1,'int32');
srccat.SubLen=fread(fid,srccat.Nsubs,'int32');
srccat.SubLen2=fread(fid,srccat.Nsubs,'int32');
srccat.CoreFrac=fread(fid,srccat.Nsubs,'float32');
srccat.SubOffset=fread(fid,srccat.Nsubs,'int32');
srccat.PSubArr=fread(fid,srccat.Nids,'int32');
srccat.PSubArr=mat2cell(srccat.PSubArr,srccat.SubLen,1);
sublen2=srccat.SubLen2;
sublen2(sublen2<0)=0;
nids2=sum(sublen2,1);
if nids2>0
    subarr2=fread(fid,nids2,'int32');
    srccat.PSubArr2=mat2cell(subarr2,sublen2,1);
else
    srccat.PSubArr2=cell(srccat.Nsubs,1);
end
srccat.NDeathSp=fread(fid,1,'int32');
fclose(fid);

srccat.property.snap=snapnum;