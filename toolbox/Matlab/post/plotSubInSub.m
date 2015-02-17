RunNum='AqA5';SnapLoad=127;grpid=0;
SnapPlot=SnapLoad;
% datadir=['/mnt/A4700/data/',RunNum,'/subcat/anal/image/'];
datadir=['/mnt/charon/HBT/data/',RunNum,'/subcat/anal/image/'];
% datadir=['/mnt/A4700/data/',RunNum,'/subcatS/anal/image/'];
% fofmap=load([datadir,'fofmapxy_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
% subpos=load([datadir,'subcen_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);  %position of most bound particle better marks subhalo than CoM
% fofsize=load([datadir,'fofsize_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
fofmap=load([datadir,'fofmap3xy_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(SnapPlot,'%03d')]);
subpos=load([datadir,'subcen_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(SnapPlot,'%03d')]);  %position of most bound particle better marks subhalo than CoM
subcom=load([datadir,'subcom_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
% subpos=subcom(:,2:4);
fofsize=load([datadir,'fofsize_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(SnapPlot,'%03d')]);
addpath('../post');

% halo=readhalo_size(['/mnt/A4700/data/',RunNum,'/subcatS/profile/logbin'],SnapLoad,'halo');
% halo=readhalo_size(['/mnt/A4700/data/',RunNum,'/subcat/profile/logbin'],SnapLoad,'halo');
halo=readhalo_size(['/mnt/charon/HBT/data/',RunNum,'/subcat/profile/logbin'],SnapLoad,'halo');
rvir=halo.Rvir(grpid+1,1);
% rvir=245.5625;
% rvir=233;
% rvir=0.0443;
% rvir=0.0599;

dim=[1,2];
massmin=20;
%% load subhalo tidal size
subid=2;
file=fullfile(datadir,['SubInSub_S',num2str(SnapLoad),'B',num2str(subid),'.COM.2'])
SpHier=importdata(file);
file=fullfile(datadir,['SubInSub_S',num2str(SnapLoad),'B',num2str(subid),'.DYN.2'])
DyHier=importdata(file);
%%
outputdir='/home/kam/Projects/HBT/presentation/baposter/images';
% outputdir='/home/kam/Projects/HBT/code/data/show/massfun/Aqua';
n1=size(SpHier,1);
nc=max(max(SpHier(:,2)),max(DyHier(:,2)))+1;
h=figure;
colors=colormap(lines(nc));
% n2=size(DyHier,1);
% colors2=colormap(lines(nc));
close(h);

figure;
cmap=contrast(log(fofmap));imagesc(fofsize(dim(1),:),fofsize(dim(2),:),(log(fofmap)'));colormap(cmap);colormap('bone');
hold on;
tc=0:0.1:2*pi+0.1;
plot(rvir*sin(tc)+subpos(1,dim(1)),rvir*cos(tc)+subpos(1,dim(2)),'r--');

for i=1:n1%size(subpos,1)
% plot(subpos(i,dim(1)),subpos(i,dim(2)),'or','MarkerSize',40./i.^0.6,'displayname',num2str(i-1));
l=plot_circle(SpHier(i,dim+2),SpHier(i,6));
set(l,'color',colors1(SpHier(i,2)+1,:),'displayname',num2str(SpHier(i,1)));
end


% for i=1:n2
% % plot(subpos(i,dim(1)),subpos(i,dim(2)),'or','MarkerSize',40./i.^0.6,'displayname',num2str(i-1));
% l=plot_circle(DyHier(i,dim+2),DyHier(i,6));
% set(l,'color',colors2(DyHier(i,2)+1,:),'linestyle','none','marker','.');
% end

axis equal
colormap('bone');
%datadir2=['/mnt/A4700/data/',RunNum,'/subcat/anal/image/'];
%fofsize2=load([datadir2,'fofsize_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
% axis([fofsize(1,:),fofsize(2,:)]);

% set(gca(),'plotboxaspectratio',diff(fofsize1,1,2));
xlim([subpos(1,dim(1))-rvir,subpos(1,dim(1))+rvir]);
ylim([subpos(1,dim(2))-rvir,subpos(1,dim(2))+rvir]);
% xlim(fofsize(dim(1),:));
% ylim([subpos(1,dim(2))-diff(fofsize(dim(1),:))/2,subpos(1,dim(2))+diff(fofsize(dim(1),:))/2])
% axis tight
set(gca,'xtick',[],'ytick',[]);%title('subfind')
% title(RunNum);
% set(gca,'ydir','reverse');
% f=[f;getframe];
% print('-depsc',[outputdir,'/halo_image',RunNum,'S',num2str(SnapLoad),'G',num2str(grpid),'.eps']);
% print('-depsc',[outputdir,'/halo_image',RunNum,'S',num2str(SnapLoad),'G',num2str(grpid),'.sf.eps']);
%hgsave([outputdir,'/halo_image_',RunNum,'.fig']);
%%
outputdir='/home/kam/Projects/HBT/presentation/baposter/images';
% outputdir='/home/kam/Projects/HBT/code/data/show/massfun/Aqua';
n=size(DyHier,1);
nc=max(max(SpHier(:,2)),max(DyHier(:,2)))+1;
h=figure;
colors=colormap(lines(nc));
close(h);

figure;
cmap=contrast(log(fofmap));imagesc(fofsize(dim(1),:),fofsize(dim(2),:),(log(fofmap)'));colormap(cmap);colormap('bone');
hold on;
tc=0:0.1:2*pi+0.1;
plot(rvir*sin(tc)+subpos(1,dim(1)),rvir*cos(tc)+subpos(1,dim(2)),'r--');

for i=1:n%size(subpos,1)
% plot(subpos(i,dim(1)),subpos(i,dim(2)),'or','MarkerSize',40./i.^0.6,'displayname',num2str(i-1));
l=plot_circle(DyHier(i,dim+2),DyHier(i,6));
set(l,'color',colors(DyHier(i,2)+1,:),'linestyle','-');
end

axis equal
colormap('bone');
%datadir2=['/mnt/A4700/data/',RunNum,'/subcat/anal/image/'];
%fofsize2=load([datadir2,'fofsize_',num2str(SnapLoad,'%03d'),'_',num2str(grpid,'%d')]);
% axis([fofsize(1,:),fofsize(2,:)]);

% set(gca(),'plotboxaspectratio',diff(fofsize1,1,2));
xlim([subpos(1,dim(1))-rvir,subpos(1,dim(1))+rvir]);
ylim([subpos(1,dim(2))-rvir,subpos(1,dim(2))+rvir]);
% xlim(fofsize(dim(1),:));
% ylim([subpos(1,dim(2))-diff(fofsize(dim(1),:))/2,subpos(1,dim(2))+diff(fofsize(dim(1),:))/2])
% axis tight
set(gca,'xtick',[],'ytick',[]);%title('subfind')
% title(RunNum);
% set(gca,'ydir','reverse');
% f=[f;getframe];
% print('-depsc',[outputdir,'/halo_image',RunNum,'S',num2str(SnapLoad),'G',num2str(grpid),'.eps']);
% print('-depsc',[outputdir,'/halo_image',RunNum,'S',num2str(SnapLoad),'G',num2str(grpid),'.sf.eps']);
%hgsave([outputdir,'/halo_image_',RunNum,'.fig']);
