clear;
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
% subid=8 at snap40,grpid=0
% subid=4484 at snap36,grpid=33,grplen_sub=51, with sat accr
% subid=4414 at snap36,grpid=33,grplen_sub=51, without sat accr
% subid=5211 at snap36,grpid=33,grplen_sub=45, with sat accr,6702dm
Nsnap=40;grpid=0;
datadir1='/mnt/A4700/data/6702DMNSM/subcat/anal/image/';
fofmap1=load([datadir1,'fofmapxy_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
% subpos1=load([datadir2,'subcen_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
subcom1=load([datadir1,'subcom_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
subpos1=subcom1(:,2:4);
fofsize1=load([datadir1,'fofsize_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
%%
Nsnap=40;grpid=0;
datadir2='/mnt/A4700/data/6702DM/subcat/anal/image/';
fofmap2=load([datadir2,'fofmapxy_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
% subpos2=load([datadir2,'subcen_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
subcom2=load([datadir2,'subcom_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
subpos2=subcom2(:,2:4);
fofsize2=load([datadir2,'fofsize_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
addpath('../post');
halo=readhalo_size('/mnt/A4700/data/6702DM/subcat/profile/logbin/',Nsnap,'halo');
rvir=halo.Rvir(grpid+1,1);
%%
figure;
h3=axes('position',[.3,.05,.6,.9]);
imagesc(fofsize1(1,:),fofsize1(2,:),log(fofmap1)');colormap('gray');
hold on;
tc=0:0.1:2*pi+0.1;
plot(rvir*sin(tc)+subpos2(1,1),rvir*cos(tc)+subpos2(1,2),'r');
axis equal
axis tight
% set(gca(),'plotboxaspectratio',diff(fofsize1,1,2));
set(gca,'xtick',[],'ytick',[]);
% grid;

h1=axes('position',[0.1,0.5,0.3,0.45]);
imagesc(fofsize1(1,:),fofsize1(2,:),log(fofmap1)');colormap('gray');hold on;

for i=1:500
plot(subpos1(i,1),subpos1(i,2),'or','MarkerSize',30./i^0.5,'displayname',num2str(i-1));
end
axis equal
axis([1.49,1.50,1.71,1.72]*1e5);
set(gca,'xtick',[],'ytick',[]);
% title('HBT without sat accretion');


h2=axes('position',[.1,.05,.3,.45]);
imagesc(fofsize2(1,:),fofsize2(2,:),log(fofmap2)');colormap('gray');hold on;

for i=1:500
plot(subpos2(i,1),subpos2(i,2),'or','MarkerSize',30./i^0.5,'displayname',num2str(i-1));
end
plot(subpos2(9,1),subpos2(9,2),'ob','MarkerSize',30./9^0.5,'displayname','lost');
axis equal
axis([1.49,1.50,1.71,1.72]*1e5);
set(gca,'xtick',[],'ytick',[]);
% title('HBT with sat accretion','color','red');

annotation('rectangle','Color','r',...
    'Position',[8/15*0.6+0.3, 8/15*0.9+0.05,1/5*0.6,1/5*0.9]);

print('-depsc',[outputdir,'/lostsub_image.eps']);
hgsave([outputdir,'/lostsub_image.fig']);
%%
Nsnap=36;grpid=33;
datadir1='/mnt/Altix/sd8/6702/subcat/NoSatMerg/anal/follow/';
lostsub1=load([datadir1,'lostsub_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
datadir2='/mnt/A4700/data/6702/subcat/anal/follow/';
lostsub2=load([datadir2,'lostsub_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
figure;
plot(lostsub2(1:31,1)-lostsub1(:,1),'s');hold on;plot(sum(lostsub2(1,2:4))-sum(lostsub2(:,2:4),2),'*');
plot(lostsub1(:,1),'.')
figure;
plot(lostsub2(:,1),'s');hold on;plot(sum(lostsub2(1,2:4))-sum(lostsub2(1:31,2:4),2)+lostsub1(:,1),'*');
plot(sum(lostsub2(:,2:5),2),'.');