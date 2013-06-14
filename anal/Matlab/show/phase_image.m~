clear;
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
outputdir='/home/kam/Projects/HBT/code/data/show';
Nsnap=99;subid=1;
datadir='/home/kam/Downloads/HBTdata/';
% datadir='/mnt/A4700/data/6702DM/subcat/anal/image/';
E=load([datadir,'energymap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.',num2str(Nsnap)]);
x=load([datadir,'relposmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.',num2str(Nsnap)]);
v=load([datadir,'relvelmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.',num2str(Nsnap)]);
% datadir='/mnt/A4700/data/6702DM/subcatS/anal/image/';
datadir='/home/kam/Downloads/HBTdata/phaseS/';
Es=load([datadir,'energymap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.',num2str(Nsnap)]);
xs=load([datadir,'relposmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.',num2str(Nsnap)]);
vs=load([datadir,'relvelmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.',num2str(Nsnap)]);
%%
r=sum(x.^2,2).^(1/2);
s=sum(v.^2,2).^(1/2);
rs=sum(xs.^2,2).^(1/2);
ss=sum(vs.^2,2).^(1/2);

% figure;plot(r,s,'.');
myfigure;
set(gcf,'DefaultLineMarkerSize',2);
subplot(1,2,1);
[xx,yy,n]=densitygrid(r(:),s(:),[50,50],[0,1500],[0,2500]);
p=log10(n+1);
l=0:0.5:3.5;
pcolor(xx,yy,p);
hold on;
contour(xx,yy,p,l);
subplot(1,2,2);
[xx,yy,n]=densitygrid(rs(:),ss(:),[50,50],[0,1500],[0,2500]);
p=log10(n+1);
pcolor(xx,yy,p);
hold on;
contour(xx,yy,p,l);
myfigure;
set(gcf,'DefaultLineMarkerSize',2);
h1=subplot(1,2,1);
plot(r,s,'.');
h2=subplot(1,2,2);
plot(rs,ss,'.');
%%
myfigure;
set(gcf,'DefaultLineMarkerSize',2);
subplot(1,2,1);
plot(x(:),v(:),'.');
axis([-1500,1500,-2500,2500]);
subplot(1,2,2);
plot(xs(:),vs(:),'.');
axis([-1500,1500,-2500,2500]);

myfigure;
set(gcf,'DefaultLineMarkerSize',2);
subplot(1,2,1);
[xx,yy,n]=densitygrid(x(:),v(:),[30,30],[-1500,1500],[-2500,2500]);
contourf(xx,yy,n);
subplot(1,2,2);
[xx,yy,n]=densitygrid(xs(:),vs(:),[30,30],[-1500,1500],[-2500,2500]);
contourf(xx,yy,n);
%%
myfigure;
set(gcf,'DefaultLineMarkerSize',2);
h1=axes('position',[.05,.5,.3,.4]);
plot(x(:,1),v(:,1),'.');
axis([-1500,1500,-2500,2500]);

h2=axes('position',[0.35,0.5,0.3,0.4]);
plot(x(:,2),v(:,2),'.');
axis([-1500,1500,-2500,2500]);

h3=axes('position',[.65,.5,.3,.4]);
plot(x(:,3),v(:,3),'.');
axis([-1500,1500,-2500,2500]);

h4=axes('position',[.05,.1,.3,.4]);
plot(xs(:,1),vs(:,1),'.');
axis([-1500,1500,-2500,2500]);

h5=axes('position',[0.35,0.1,0.3,0.4]);
plot(xs(:,2),vs(:,2),'.');
axis([-1500,1500,-2500,2500]);

h6=axes('position',[.65,.1,.3,.4]);
plot(xs(:,3),vs(:,3),'.');
axis([-1500,1500,-2500,2500]);


% print('-depsc',[outputdir,'/lostsub_image.eps']);
%%
myfigure;
set(gcf,'DefaultLineMarkerSize',2);
h1=subplot(1,2,1);
plot(E(:,1),r,'.');
axis([0,4e6,0,1600]);
h2=subplot(1,2,2);
plot(Es(:,1),rs,'.');
axis([0,4e6,0,1600]);
%%
myfigure;
set(gcf,'DefaultLineMarkerSize',2);
plot(x(:,1),v(:,2),'.');
hold on;
plot(xs(:,1),vs(:,2),'r.');
%% read phase map for SF,BT-only and Background particles.
outputdir='/home/kam/Projects/HBT/code/data/show';
Nsnap=99;subid=1;
datadir='/home/kam/Downloads/HBTdata/';
% datadir='/mnt/A4700/data/6702DM/subcat/anal/image/';
fid=fopen([datadir,'posmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.SF']);
x=fread(fid,inf,'float32');
xSF=reshape(x,3,numel(x)/3)';
fclose(fid);
fid=fopen([datadir,'velmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.SF']);
v=fread(fid,inf,'float32');
vSF=reshape(v,3,numel(v)/3)';
fclose(fid);
fid=fopen([datadir,'posmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.BT']);
x=fread(fid,inf,'float32');
xBT=reshape(x,3,numel(x)/3)';
fclose(fid);
fid=fopen([datadir,'velmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.BT']);
v=fread(fid,inf,'float32');
vBT=reshape(v,3,numel(v)/3)';
fclose(fid);
fid=fopen([datadir,'posmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.BK']);
x=fread(fid,inf,'float32');
xBK=reshape(x,3,numel(x)/3)';
fclose(fid);
fid=fopen([datadir,'velmap_',num2str(Nsnap,'%03d'),'_',num2str(subid,'%d'),'.BK']);
v=fread(fid,inf,'float32');
vBK=reshape(v,3,numel(v)/3)';
fclose(fid);
xcSF=xSF(1,:);vcSF=vSF(1,:);
xcBT=xBT(1,:);vcBT=vBT(1,:);
xSF(1,:)=[];
vSF(1,:)=[];
xBT(1,:)=[];
vBT(1,:)=[];
xBK(1,:)=[];
vBK(1,:)=[];
rBK=sqrt(sum((xBK-repmat(xcBT,size(xBK,1),1)).^2,2));
sBK=sqrt(sum((vBK-repmat(vcBT,size(vBK,1),1)).^2,2));
rBT=sqrt(sum((xBT-repmat(xcBT,size(xBT,1),1)).^2,2));
sBT=sqrt(sum((vBT-repmat(vcBT,size(vBT,1),1)).^2,2));
rSF=sqrt(sum((xSF-repmat(xcBT,size(xSF,1),1)).^2,2));
sSF=sqrt(sum((vSF-repmat(vcBT,size(vSF,1),1)).^2,2));
%%
figure;
set(gcf,'defaultlinemarkersize',2);
h1=subplot(1,3,1);
plot(xBK(:,1),vBK(:,1),'b.');hold on;
h2=subplot(1,3,2);
plot(xBT(:,1),vBT(:,1),'y.');
h3=subplot(1,3,3);
plot(xSF(:,1),vSF(:,1),'r.');
set([h2,h3],'xlim',get(h1,'xlim'),'ylim',get(h1,'ylim'));
%%
myfigure;
set(gcf,'defaultlinemarkersize',2);
plot(rBK(1:100:end),sBK(1:100:end),'b.');hold on;
[xx,yy,n]=densitygrid(rBK,sBK,[50,50],[0,1600],[0,7000]);
[c,h1]=contour(xx,yy,n,'b');hold on;

% plot(rBT(1:100:end),sBT(1:100:end),'k.');
% [xx,yy,n]=densitygrid(rBT,sBT,[50,50],[0,1600],[0,7000]);
[xx,yy,n]=densitygrid([rBT;rSF],[sBT;sSF],[50,50],[0,1600],[0,7000]);
[c,h2]=contour(xx,yy,n,'k');

% plot(rSF(1:100:end),sSF(1:100:end),'r.');
[xx,yy,n]=densitygrid(rSF,sSF,[50,50],[0,1600],[0,7000]);
[c,h3]=contour(xx,yy,n,'r');

xlabel('r/(kpc/h)');
ylabel('v/(km/s)');
legend([h2,h3,h1],'HBT','SUBFIND','Other particles');
print('-depsc',[outputdir,'/phase_map_6702DM_SN99SUB1.eps']);
%%
figure;
subplot(1,3,1);
[xx,yy,n]=densitygrid(rBK,sBK,[50,50],[0,1600],[0,7000]);
p=log10(n+1);
l=0:0.5:3.5;
pcolor(xx,yy,p);
hold on;
contour(xx,yy,p,l);
subplot(1,3,2);
[xx,yy,n]=densitygrid(rBT,sBT,[50,50],[0,1600],[0,7000]);
p=log10(n+1);
l=0:0.5:4;
% pcolor(xx,yy,p);
% hold on;
contour(xx,yy,p,l);
subplot(1,3,3);
[xx,yy,n]=densitygrid(rSF,sSF,[50,50],[0,1600],[0,7000]);
p=log10(n+1);
l=0:0.5:3.5;
pcolor(xx,yy,p);
hold on;
contour(xx,yy,p,l);