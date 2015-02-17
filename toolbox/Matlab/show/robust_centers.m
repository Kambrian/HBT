RunNum='6702DM';Nsnap=51;grpid=86;
datadir=['/mnt/A4700/data/',RunNum,'/subcat/anal/image/'];
% datadir=['/mnt/A4700/data/',RunNum,'/subcatS/anal/image/'];
fofmap=load([datadir,'fofmapxy_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
%subpos=load([datadir,'subcen_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);  %position of most bound particle better marks subhalo than CoM
subcom=load([datadir,'subcom_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
subpos=subcom(:,2:4);
fofsize=load([datadir,'fofsize_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
addpath('../post');

% halo=readhalo_size(['/mnt/A4700/data/',RunNum,'/subcatS/profile/logbin'],Nsnap,'halo');
halo=readhalo_size(['/mnt/A4700/data/',RunNum,'/subcat/profile/logbin'],Nsnap,'halo');
rvir=halo.Rvir(grpid+1,1);
%%
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
figure;%('visible','off');
imagesc(fofsize(1,:),fofsize(2,:),log(fofmap)');colormap('gray');
hold on;
tc=0:0.1:2*pi+0.1;
plot(rvir*sin(tc)+subpos(1,1),rvir*cos(tc)+subpos(1,2),'r--');

for i=1:size(subpos,1)
plot(subpos(i,1),subpos(i,2),'or','MarkerSize',30./i^0.5,'displayname',num2str(i-1));
end

axis equal

%datadir2=['/mnt/A4700/data/',RunNum,'/subcat/anal/image/'];
%fofsize2=load([datadir2,'fofsize_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
%axis([fofsize2(1,:),fofsize2(2,:)]);
% axis tight
% set(gca(),'plotboxaspectratio',diff(fofsize1,1,2));
set(gca,'xtick',[],'ytick',[]);
title(RunNum);

% print('-depsc',[outputdir,'/halo_image_100_',RunNum,'.eps']);
%hgsave([outputdir,'/halo_image_',RunNum,'.fig']);
