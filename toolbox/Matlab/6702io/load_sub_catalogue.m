function subcat=load_sub_catalogue(snapnum)
% function subcat=load_sub_catalogue(snapnum)

global subcatdir
% struct Chain_data 
% {
% int ProSubID;
% int HostID;
% };
% 
% struct Hierarchy
% {
% int nibs;	//its host sub
% int pre;	//its previous sibling
% int next;	//its next sibling
% int sub;	//its sub-in-sub
% } ;
% struct SubProperty
% {
% 	float CoM[3];//center of mass position (comoving)
% 	float VCoM[3];//CoM velocity (physical)
% 	float Pot;//average gavitational potential per unit mass per particle (GM/R_physical)
% 	float Kin;//average Kinetic energy per unit mass per particle (0.5*V_physical^2)
% 	float AM[3];//average angular momentum per unit mass per particle (R_physical x V_physical)
% };
% typedef struct
% {
% 	int Ngroups;
% 	int Nsubs;
% 	int Nids;
% 	int *GrpLen_Sub;//number of subs in a group
% 	int *GrpOffset_Sub;//sub index offset of a group
% 	int *SubLen;
% 	int *SubOffset;
% 
% 	int *SubRank;
% 	struct Chain_data *HaloChains;
% 	
% 	struct Hierarchy *sub_hierarchy;
% 	
% 	struct SubProperty *Property;   
% 	
% 	int **PSubArr;//PSubArr[Nsubs] as a pointer will point to PIDs block (or PIndices during unbind()) of each sub;i.e.,PSubArr[subhaloid][particle] will give PIDs(or PInd);
% 	
% 	int Nbirth;
% 	int NQuasi;
% 	int Ndeath;
% 	int Nsplitter;//the splitted-out subs from last snap
% } SUBCATALOGUE;

subcat=struct('Ngroups',0,'Nsubs',0,'Nids',0,'GrpLen_Sub',[],'GrpOffset_Sub',[],'SubLen',[],'SubOffset',[],'SubRank',[],...
                'HaloChains',struct('ProSubID',[],'HostID',[]),'sub_hierarchy',struct('nibs',[],'pre',[],'next',[],'sub',[]),...
                'SubProp',struct('CoM',[],'VCoM',[],'Pot',[],'Kin',[],'AM',[]),'PSubArr',[],'Nbirth',0,'NQuasi',0,'Ndeath',0,'Nsplitter',0,...
                'property',struct('type','subcat','particles','id','snap',-1));

subfile=fullfile(subcatdir,['subcat_',num2str(snapnum,'%03d')]);
fid=fopen(subfile,'r');        
subcat.Ngroups=fread(fid,1,'int32');
subcat.Nsubs=fread(fid,1,'int32');
subcat.Nids=fread(fid,1,'int32');
subcat.GrpLen_Sub=fread(fid,subcat.Ngroups,'int32');
subcat.GrpOffset_Sub=fread(fid,subcat.Ngroups,'int32');
subcat.SubLen=fread(fid,subcat.Nsubs,'int32');
subcat.SubOffset=fread(fid,subcat.Nsubs,'int32');
subcat.SubRank=fread(fid,subcat.Nsubs,'int32');
halochains=fread(fid,[2,subcat.Nsubs],'int32');
subcat.HaloChains.ProSubID=halochains(1,:)';
subcat.HaloChains.HostID=halochains(2,:)';
% clear halochains;
sub_hier=fread(fid,[4,subcat.Nsubs],'int32');
subcat.sub_hierarchy.nibs=sub_hier(1,:)';
subcat.sub_hierarchy.pre=sub_hier(2,:)';
subcat.sub_hierarchy.next=sub_hier(3,:)';
subcat.sub_hierarchy.sub=sub_hier(4,:)';
% clear sub_hier;
subprop=fread(fid,[3+3+1+1+3,subcat.Nsubs],'float32');
subcat.SubProp.CoM=subprop(1:3,:)';
subcat.SubProp.VCoM=subprop(4:6,:)';
subcat.SubProp.Pot=subprop(7,:)';
subcat.SubProp.Kin=subprop(8,:)';
subcat.SubProp.AM=subprop(9:11,:)';
% clear subprop;
subcat.PSubArr=fread(fid,subcat.Nids,'int32');
subcat.PSubArr=mat2cell(subcat.PSubArr,subcat.SubLen,1);
subcat.Nbirth=fread(fid,1,'int32');
subcat.NQuasi=fread(fid,1,'int32');
subcat.Ndeath=fread(fid,1,'int32');
subcat.Nsplitter=fread(fid,1,'int32');
fclose(fid);

subcat.property.snap=snapnum;
