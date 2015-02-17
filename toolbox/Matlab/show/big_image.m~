clear;
Nsnap=99;grpid=0;

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';

datadir='/mnt/A4700/data/6600/subcat/anal/image/';
bigmap=load([datadir,'bigmapxy_',num2str(Nsnap,'%03d')]);
bigsize=load([datadir,'bigsize_',num2str(Nsnap,'%03d')]);
%%
fofmap=load([datadir,'fofmapxy_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
fofbmap=load([datadir,'fofbmapxy_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
% subpos=load([datadir,'subcen_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
subcom=load([datadir,'subcom_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
subpos=subcom(:,2:4);
fofsize=load([datadir,'fofsize_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
fofbsize=load([datadir,'fofbsize_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
addpath('../post');
halo=readhalo_size('/mnt/A4700/data/6600/subcat/profile/logbin/',Nsnap,'halo');
rvir=halo.Rvir(grpid+1,1);
%%
figure;
h3=axes('position',[.3,.05,.6,.9]);
imagesc(bigsize(1,:),bigsize(2,:),log(bigmap)');colormap('gray');
hold on;
tc=0:0.1:2*pi+0.1;
plot(rvir*sin(tc)+subpos(1,1),rvir*cos(tc)+subpos(1,2),'r');
axis equal
axis tight
% set(gca(),'plotboxaspectratio',diff(fofsize1,1,2));
set(gca,'xtick',[],'ytick',[]);
% grid;

h1=axes('position',[0.1,0.5,0.3,0.45]);
imagesc(fofsize(1,:),fofsize(2,:),log(fofmap)');colormap('gray');hold on;

for i=1:100
plot(subpos(i,1),subpos(i,2),'og','MarkerSize',30./i^0.5,'displayname',num2str(i-1));
end
hold on;
plot(rvir*sin(tc)+subpos(1,1),rvir*cos(tc)+subpos(1,2),'r');
axis equal
set(gca,'xtick',[],'ytick',[]);
% title('HBT without sat accretion');

h2=axes('position',[.1,.05,.3,.45]);
imagesc(fofbsize(1,:),fofbsize(2,:),(fofbmap)');colormap('gray');hold on;

% for i=1:100
% plot(subpos(i,1),subpos(i,2),'og','MarkerSize',30./i^0.5,'displayname',num2str(i-1));
% end
plot(rvir*sin(tc)+subpos(1,1),rvir*cos(tc)+subpos(1,2),'r');
axis equal
set(gca,'xtick',[],'ytick',[]);

% annotation('rectangle','Color','r',...
%     'Position',[8/15*0.6+0.3, 8/15*0.9+0.05,1/5*0.6,1/5*0.9]);
% 
% print('-depsc',[outputdir,'/lostsub_image.eps']);
% hgsave([outputdir,'/lostsub_image.fig']);