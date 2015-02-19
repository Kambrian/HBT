clear;
addpath('../post');

runnum='8213';Nsnap=59;

datadir=['/mnt/A4700/data/',runnum,'/subcat/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(Nsnap,'%03d')]);
btsize=load_subsize(sizefile);

datadir=['/mnt/A4700/data/',runnum,'/subcatS/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(Nsnap,'%03d')]);
sfsize=load_subsize(sizefile);

f=@(x) x.*(3*(log(1+x)-x./(1+x))-(x./(1+x)).^2).^(-1.0/3);%ritdal_sub/rvir_sub as a function of rcen_sub/rs_host
%%
nplotBT=1:btsize.nsubs;
nplotBT=nplotBT(~btsize.flag_badbins(nplotBT)&logical(btsize.rcen(nplotBT)));%all the satellites
nplotSF=1:sfsize.nsubs;
nplotSF=nplotSF(~btsize.flag_badbins(nplotSF)&logical(btsize.rcen(nplotSF)));%all the satellites
%==calculate median line==%
x=logspace(1,3,10);
[n,bin]=histc(btsize.rtidal(nplotBT),x);
bty1=zeros(9,1);bty02=bty1;btxm=bty1;
for i=1:9
bty1(i)=median(btsize.req_bk_1(nplotBT(bin==i)));
bty02(i)=median(btsize.req_bk_02(nplotBT(bin==i)));
btxm(i)=median(btsize.rtidal(nplotBT(bin==i)));
end
x=logspace(1,3,10);[n,bin]=histc(sfsize.rtidal(nplotSF),x);
sfy1=zeros(9,1);sfy02=sfy1;sfxm=sfy1;
for i=1:9
sfy1(i)=median(sfsize.req_bk_1(nplotSF(bin==i)));
sfy02(i)=median(sfsize.req_bk_02(nplotSF(bin==i)));
sfxm(i)=median(sfsize.rtidal(nplotSF(bin==i)));
end
%===plots===%
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
set(gcf,'DefaultLineMarkerSize',5);
axes('position',[.1,.5,.4,.4]);loglog(btsize.rtidal(nplotBT),btsize.req_bk_1(nplotBT),'.');hold on;loglog([10,1e3],[10,1e3]);
ylabel('$Req/(kpc/h)$');set(gca,'xaxislocation','top','xtick',[10,100],'ytick',[10,100]);axis([10,1e3,10,1e3]);text(400,20,'HBT');
loglog(btxm,bty1,'r');
axes('position',[.5,.5,.4,.4]);loglog(sfsize.rtidal(nplotSF),sfsize.req_bk_1(nplotSF),'.');hold on;loglog([10,1e3],[10,1e3]);
ylabel('$Req/(kpc/h)$');set(gca,'xaxislocation','top','xtick',[10,100,1000],'yaxislocation','right');axis([10,1e3,10,1e3]);text(400,20,'HBT');
loglog(sfxm,sfy1,'r');
axes('position',[.1,.1,.4,.4]);loglog(btsize.rtidal(nplotBT),btsize.req_bk_02(nplotBT),'.');hold on;loglog([10,1e3],[10,1e3]);
ylabel('$Req_{0.02}/(kpc/h)$');axis([10,1e3,10,1e3]);set(gca,'ytick',[10,100],'xtick',[10,100]);xlabel('$Rtidal/(kpc/h)$');text(300,20,'SUBFIND');
loglog(btxm,bty02,'r');
loglog(sfxm,sfy02,'c--');
axes('position',[.5,.1,.4,.4]);loglog(sfsize.rtidal(nplotSF),sfsize.req_bk_02(nplotSF),'.');hold on;loglog([10,1e3],[10,1e3]);
ylabel('$Req_{0.02}/(kpc/h)$');set(gca,'yaxislocation','right','ytick',[10,100],'xtick',[10,100,1000]);axis([10,1e3,10,1e3]);text(300,20,'SUBFIND');
xlabel('$Rtidal/(kpc/h)$');
loglog(sfxm,sfy02,'r');

set(gcf,'paperposition',[0.6,6,20,18]);
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/subsize_bk_aqua_',runnum,'.eps']);
%%
clear;
addpath('../post');

runnum='6702DM';Nsnap=99;

datadir=['/mnt/A4700/data/',runnum,'/subcat/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(Nsnap,'%03d')]);
btsize=load_subsize(sizefile);

datadir=['/mnt/A4700/data/',runnum,'/subcatS/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(Nsnap,'%03d')]);
sfsize=load_subsize(sizefile);
%%
% req_02 for SUBFIND is 20~30 percent smaller than HBT
nplotBT=1:btsize.nsubs;
nplotBT=nplotBT(~btsize.flag_badbins(nplotBT)&logical(btsize.rcen(nplotBT)));%all the satellites
nplotSF=1:sfsize.nsubs;
nplotSF=nplotSF(~btsize.flag_badbins(nplotSF)&logical(btsize.rcen(nplotSF)));%all the satellites
%==calculate median line==%
x=logspace(1,3,10);
[n,bin]=histc(btsize.rtidal(nplotBT),x);
bty1=zeros(9,1);bty02=bty1;btxm=bty1;
for i=1:9
bty1(i)=median(btsize.req_bk_1(nplotBT(bin==i)));
bty02(i)=median(btsize.req_bk_02(nplotBT(bin==i)));
btxm(i)=median(btsize.rtidal(nplotBT(bin==i)));
end
x=logspace(1,3,10);[n,bin]=histc(sfsize.rtidal(nplotSF),x);
sfy1=zeros(9,1);sfy02=sfy1;sfxm=sfy1;
for i=1:9
sfy1(i)=median(sfsize.req_bk_1(nplotSF(bin==i)));
sfy02(i)=median(sfsize.req_bk_02(nplotSF(bin==i)));
sfxm(i)=median(sfsize.rtidal(nplotSF(bin==i)));
end
%===plots===%
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
set(gcf,'DefaultLineMarkerSize',5);
axes('position',[.1,.5,.4,.4]);loglog(btsize.rtidal(nplotBT),btsize.req_bk_1(nplotBT),'.');hold on;loglog([10,1e3],[10,1e3]);
ylabel('$Req/(kpc/h)$');set(gca,'xaxislocation','top','xtick',[10,100],'ytick',[10,100]);axis([10,1e3,10,1e3]);text(400,20,'HBT');
loglog(btxm,bty1,'r');
axes('position',[.5,.5,.4,.4]);loglog(sfsize.rtidal(nplotSF),sfsize.req_bk_1(nplotSF),'.');hold on;loglog([10,1e3],[10,1e3]);
ylabel('$Req/(kpc/h)$');set(gca,'xaxislocation','top','xtick',[10,100,1000],'yaxislocation','right');axis([10,1e3,10,1e3]);text(400,20,'HBT');
loglog(sfxm,sfy1,'r');
axes('position',[.1,.1,.4,.4]);loglog(btsize.rtidal(nplotBT),btsize.req_bk_02(nplotBT),'.');hold on;loglog([10,1e3],[10,1e3]);
ylabel('$Req_{0.02}/(kpc/h)$');axis([10,1e3,10,1e3]);set(gca,'ytick',[10,100],'xtick',[10,100]);xlabel('$Rtidal/(kpc/h)$');text(300,20,'SUBFIND');
loglog(btxm,bty02,'r');
loglog(sfxm,sfy02,'c--');
axes('position',[.5,.1,.4,.4]);loglog(sfsize.rtidal(nplotSF),sfsize.req_bk_02(nplotSF),'.');hold on;loglog([10,1e3],[10,1e3]);
ylabel('$Req_{0.02}/(kpc/h)$');set(gca,'yaxislocation','right','ytick',[10,100],'xtick',[10,100,1000]);axis([10,1e3,10,1e3]);text(300,20,'SUBFIND');
xlabel('$Rtidal/(kpc/h)$');
loglog(sfxm,sfy02,'r');

set(gcf,'paperposition',[0.6,6,20,18]);
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/subsize_bk_aqua_',runnum,'.eps']);