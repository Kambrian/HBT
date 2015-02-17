RunNum='8213';Nsnap=59;grpid=3;
datadir=['/mnt/A4700/data/',RunNum,'/subcat/anal/image/'];
% fofmap=load([datadir,'fofmap3yz_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(Nsnap,'%03d')]);
fofmap=load([datadir,'fofmap2xz_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(Nsnap,'%03d')]);
% fofmap=load([datadir,'fofmapxz_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
btpos=load([datadir,'subcen_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);  %position of most bound particle better marks subhalo than CoM
% subcom=load([datadir,'subcom_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
% subpos=subcom(:,2:4);
% fofsize=load([datadir,'fofsize_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);
fofsize=load([datadir,'fofsize_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d'),'.',num2str(Nsnap,'%03d')]);

datadir=['/mnt/A4700/data/',RunNum,'/subcatS/anal/image/'];
sfpos=load([datadir,'subcen_',num2str(Nsnap,'%03d'),'_',num2str(grpid,'%d')]);

% halo=readhalo_size(['/mnt/A4700/data/',RunNum,'/subcatS/profile/logbin'],Nsnap,'halo');
% halo=readhalo_size(['/mnt/A4700/data/',RunNum,'/subcat/profile/logbin'],Nsnap,'halo');
% rvir=halo.Rvir(grpid+1,1);
dims=[1,3];
%%
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
outputdir='/home/kam/Projects/HBT/code/data/show';

figure;
I=imfilter((fofmap),fspecial('disk',2));
I=log(I+1);
% I=histeq(I);
% cmap=contrast(I);
imagesc(fofsize(dims(1),:),fofsize(dims(2),:),I');
% colormap(cmap);
colormap('gray');
hold on;
%%
I=bin_image(fofmap,4);
I=log(I);
% I=histeq(I);
imagesc(fofsize(dims(1),:),fofsize(dims(2),:),I');
cmap=contrast(log(fofmap+1));
imagesc(fofsize(dims(1),:),fofsize(dims(2),:),(log(fofmap+1)'));
colormap(cmap);
colormap('gray');
hold on;
% tc=0:0.1:2*pi+0.1;
% plot(rvir*sin(tc)+subpos(1,1),rvir*cos(tc)+subpos(1,2),'r--');
%%
plot(btpos(2,dims(1)),btpos(2,dims(2)),'ro','markersize',10);%subid 2571
plot(btpos(7,dims(1)),btpos(7,dims(2)),'rx');%subid 2576
plot(sfpos(2,dims(1)),sfpos(2,dims(2)),'go','markersize',10); %subid 2196
plot(sfpos(33,dims(1)),sfpos(33,dims(2)),'gx');%subid 2221
axis equal
axis([fofsize(dims(1),:),fofsize(dims(2),:)]);
set(gca,'xtick',[],'ytick',[]);
% title(RunNum);

% print('-depsc',[outputdir,'/match_switch.eps']);
% print('-depsc',[outputdir,'/halo_image',RunNum,'S',num2str(Nsnap),'G',num2str(grpid),'.sf.eps']);
%hgsave([outputdir,'/halo_image_',RunNum,'.fig']);
