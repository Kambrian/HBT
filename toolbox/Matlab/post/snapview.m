RunNum='CoCo';SnapLoad=0;grpid=0;
SnapPlot=SnapLoad;
% datadir=['/mnt/uv/HBT/data/',RunNum,'/subcat/anal/image/'];
% datadir=['/mnt/A4700/data/',RunNum,'/subcat/anal/image/'];
datadir=['/mnt/charon/HBT/data/',RunNum,'/test/anal/image/'];
% datadir=['/mnt/A4700/data/',RunNum,'/subcatS/anal/image/'];
% fofmap=load([datadir,'fofmapxy_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
% subpos=load([datadir,'subcen_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);  %position of most bound particle better marks subhalo than CoM
% subcom=load([datadir,'subcom_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
% fofsize=load([datadir,'fofsize_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
fofmap=load([datadir,'fofmap2xy_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(SnapPlot,'%03d')]);
subpos=load([datadir,'subcen_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(SnapPlot,'%03d')]);  %position of most bound particle better marks subhalo than CoM
subcom=load([datadir,'subcom_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
% subpos=subcom(:,2:4);
fofsize=load([datadir,'fofsize_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(SnapPlot,'%03d')]);
addpath('../post');

% halo=readhalo_size(['/mnt/A4700/data/',RunNum,'/subcatS/profile/logbin'],SnapLoad,'halo');
% halo=readhalo_size(['/mnt/A4700/data/',RunNum,'/subcat/profile/logbin'],SnapLoad,'halo');
% halo=readhalo_size(['/mnt/charon/HBT/data/',RunNum,'/subcat/profile/logbin'],SnapLoad,'halo');
% halo=readhalo_size([datadir,'../../profile/logbin'],SnapLoad,'halo');
% rvir=halo.Rvir(grpid+1,1);
% rvir=245.5625;
% rvir=233;
% rvir=0.0443;
% rvir=0.0599;
% rvir=1609.3;
% rvir=1517.0;
% rivr=1452;
% rvir=nan;

dim=[1,2];
massmin=10;
%%
subid=0;
% file=fullfile(datadir,['SubInSub_S',num2str(SnapLoad),'B',num2str(subid),'.COM.2'])
% SpHier=importdata(file);
file=fullfile(datadir,['SubInSub_S',num2str(SnapLoad),'B',num2str(subid),'.DYN.2'])
DyHier=importdata(file);
R2sig=DyHier(:,end);
%% load subhalo tidal size
datadir=['/mnt/A4700/data/',runnum,'/subcat/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(SnapLoad,'%03d')]);
btsize=load_subsize(sizefile);
sizefile=fullfile(datadir,basedir,['main_size_',num2str(SnapLoad,'%03d')]);
btmsize=load_subsize(sizefile);
% rtidal=[btmsize.rtidal(grpid+1);btsize.rtidal(btsize.id>btmsize.id(grpid+1)&btsize.id<btmsize.id(grpid+2))];

datadir=['/mnt/A4700/data/',runnum,'/subcatS/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(SnapLoad,'%03d')]);
sfsize=load_subsize(sizefile);
sizefile=fullfile(datadir,basedir,['main_size_',num2str(SnapLoad,'%03d')]);
sfmsize=load_subsize(sizefile);
% rtidal=[sfmsize.rtidal(grpid+1);sfsize.rtidal(sfsize.id>sfmsize.id(grpid+1)&sfsize.id<sfmsize.id(grpid+2))];
%%
% outputdir='/home/kam/Projects/HBT/presentation/baposter/images';
outputdir='/work/Projects/HBT/code/data/show/images';
myfigure;%('visible','off');
cmap=contrast(log(fofmap));imagesc(fofsize(dim(1),:),fofsize(dim(2),:),(log(fofmap)'));colormap(cmap);colormap('bone');
% imagesc(fofsize(1,:),fofsize(2,:),(log(fofmap+1)'));colormap('gray');
hold on;
tc=0:0.1:2*pi+0.1;
plot(rvir*sin(tc)+subpos(1,dim(1)),rvir*cos(tc)+subpos(1,dim(2)),'r--');

% for i=[2,3,8,9]  % grpid=11  overlap
% for i=[7,13]  % grpid=0 switch
% for i=[3,14]  % grpid=3
% for i=[5,24]  % grpid=26
for i=1:sum(subcom(:,1)>massmin)
plot(subpos(i,dim(1)),subpos(i,dim(2)),'or','MarkerSize',40./i.^0.6,'displayname',num2str(i-1));
% plot(rtidal(i)*sin(tc)+subpos(i,1),rtidal(i)*cos(tc)+subpos(i,2),'r-');
end
% plot(subpos(4,1),subpos(4,3),'xr','MarkerSize',8,'displayname',num2str(i-1));

axis equal

%datadir2=['/mnt/A4700/data/',RunNum,'/subcat/anal/image/'];
%fofsize2=load([datadir2,'fofsize_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
% axis([fofsize(1,:),fofsize(2,:)]);

% set(gca(),'plotboxaspectratio',diff(fofsize1,1,2));
% xlim([subpos(1,dim(1))-rvir,subpos(1,dim(1))+rvir]);
% ylim([subpos(1,dim(2))-rvir,subpos(1,dim(2))+rvir]);
% xlim(fofsize(dim(1),:));
% ylim([subpos(1,dim(2))-diff(fofsize(dim(1),:))/2,subpos(1,dim(2))+diff(fofsize(dim(1),:))/2])
axis tight
% set(gca,'xtick',[],'ytick',[]);title(RunNum)
title(RunNum);
% set(gca,'ydir','reverse');
% f=[f;getframe];
% print('-dpng',[outputdir,'/halo_image',RunNum,'S',num2str(SnapLoad),'G',num2str(grpid),'.','dim',num2str(dim(1)),num2str(dim(2)),'.png']);
% print('-depsc',[outputdir,'/halo_image',RunNum,'S',num2str(SnapLoad),'G',num2str(grpid),'.sf.eps']);
%hgsave([outputdir,'/halo_image_',RunNum,'.fig']);
