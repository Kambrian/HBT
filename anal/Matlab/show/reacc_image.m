clear;
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
outputdir='/home/kam/Projects/HBT/code/data/show';
Nsnap=58;subid=7;
datadir='/mnt/A4700/data/6702DM/subcat/anal/image/';
inmap=load([datadir,'inmapxy_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d')]);
outmap=load([datadir,'outmapxy_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d')]);
outsize=load([datadir,'outsize_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d')]);
%%
figure;
h1=axes('position',[.05,.05,.3,.9]);
imagesc(outsize(1,:),outsize(2,:),log(inmap+outmap)');%colormap('gray');
axis equal
axis tight
set(gca,'xtick',[],'ytick',[]);

h2=axes('position',[0.35,0.05,0.3,0.9]);
imagesc(outsize(1,:),outsize(2,:),log(inmap)');%colormap('gray');

axis equal
axis tight
set(gca,'xtick',[],'ytick',[]);


h3=axes('position',[.65,.05,.3,.9]);
imagesc(outsize(1,:),outsize(2,:),log(outmap)');colormap('gray');
axis equal
axis tight
set(gca,'xtick',[],'ytick',[]);


% print('-depsc',[outputdir,'/lostsub_image.eps']);
%%
clear;
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
Nsnap=58;subid=7;
datadir='/mnt/A4700/data/6702DM/subcat/anal/image/';
indata=load([datadir,'inpart_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d')]);
outdata=load([datadir,'outpart_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d')]);
InPos=indata(:,1:3);
InKin=indata(:,4);
InPot=indata(:,5);
OutPos=outdata(:,1:3);
OutKin=outdata(:,4);
OutPot=outdata(:,5);
Pm=-mean(InPot);
InKin=InKin/Pm;
InPot=-InPot/Pm;
OutKin=OutKin/Pm;
OutPot=-OutPot/Pm;
%%
figure;
plot3(OutPos(:,1),OutPos(:,2),OutPos(:,3),'b.','markersize',2);
hold on;
plot3(InPos(:,1),InPos(:,2),InPos(:,3),'r.','markersize',2);
hold on;

figure;plot(InKin,InPot,'.');
figure;plot(OutKin,OutPot,'.');

x=0:0.03:3;
figure;hist(InKin,x);
figure;hist(OutKin,x);

figure;hist(InPot,100);
figure;hist(OutPot,100);

figure;hist(InKin./InPot,x);
figure;hist(log10(OutKin./OutPot),100);
figure;hist(InKin-InPot,50);
figure;hist(OutKin-OutPot,50);


