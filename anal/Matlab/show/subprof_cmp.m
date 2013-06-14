clear;
addpath('../post');

TIDALMAX=2;
runnum='6702DM';Nsnap=99;SofteningHalo=1.5;
% runnum='8213';Nsnap=59;SofteningHalo=3.3;

datadir=['/mnt/A4700/data/',runnum,'/subcat/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(Nsnap,'%03d'),'.',num2str(TIDALMAX)]);
proffile=fullfile(datadir,basedir,['sat_prof_',num2str(Nsnap,'%03d'),'.',num2str(TIDALMAX)]);
btsize=load_subsize(sizefile);
btprof=load_subprof(proffile,btsize.nbin);
sizefile=fullfile(datadir,basedir,['main_size_',num2str(Nsnap,'%03d'),'.',num2str(TIDALMAX)]);
proffile=fullfile(datadir,basedir,['main_prof_',num2str(Nsnap,'%03d'),'.',num2str(TIDALMAX)]);
btmsize=load_subsize(sizefile);
btmprof=load_subprof(proffile,btmsize.nbin);

datadir=['/mnt/A4700/data/',runnum,'/subcatS/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['sat_size_',num2str(Nsnap,'%03d'),'.',num2str(TIDALMAX)]);
proffile=fullfile(datadir,basedir,['sat_prof_',num2str(Nsnap,'%03d'),'.',num2str(TIDALMAX)]);
sfsize=load_subsize(sizefile);
sfprof=load_subprof(proffile,sfsize.nbin);
sizefile=fullfile(datadir,basedir,['main_size_',num2str(Nsnap,'%03d'),'.',num2str(TIDALMAX)]);
proffile=fullfile(datadir,basedir,['main_prof_',num2str(Nsnap,'%03d'),'.',num2str(TIDALMAX)]);
sfmsize=load_subsize(sizefile);
sfmprof=load_subprof(proffile,sfmsize.nbin);

%% sat profile
myfigure;
btsubind=1;sfsubind=1;
% btsubind=727;sfsubind=636;
% btsubind=254;sfsubind=223;
vs=4*pi/3*diff(logspace(log10(max(SofteningHalo,btsize.rmax(btsubind)*1e-2)),log10(btsize.rmax(btsubind)),btsize.nbin(btsubind)+1).^3)';
svs=cumsum(vs);
na=btprof.no{btsubind}+btprof.nb{btsubind};
sna=cumsum(na);
sns=cumsum(btprof.ns{btsubind});
% rscale=btsize.rtidal(btsubind);
rscale=1;
% rscale=btsize.req_all_1(btsubind);
% rscale=sfsize.rtidal(sfsubind);

loglog(btprof.rs{btsubind}/rscale,btprof.ns{btsubind}./vs,'r-','displayname','HBT sub');
hold on;
loglog(btprof.ro{btsubind}/rscale,(btprof.no{btsubind}+btprof.nb{btsubind})./vs,'g-','displayname','HBT back');
% loglog(btprof.ro{btsubind}/rscale,btprof.no{btsubind}./vs,'g-','displayname','HBT other');
% loglog(btprof.rb{btsubind}/rscale,btprof.nb{btsubind}./vs,'b-','displayname','HBT back');
% loglog(btprof.rs{btsubind}/rscale,na./vs,'k-. .','displayname','HBT other+back');

vs=4*pi/3*diff(logspace(log10(max(SofteningHalo,sfsize.rmax(sfsubind)*1e-2)),log10(sfsize.rmax(sfsubind)),sfsize.nbin(sfsubind)+1).^3)';
svs=cumsum(vs);
na=sfprof.no{sfsubind}+sfprof.nb{sfsubind};
sna=cumsum(na);
sns=cumsum(sfprof.ns{sfsubind});
loglog(sfprof.rs{sfsubind}/rscale,sfprof.ns{sfsubind}./vs,'r--','displayname','SF sub');
loglog(sfprof.ro{sfsubind}/rscale,(sfprof.no{sfsubind}+sfprof.nb{sfsubind})./vs,'g--','displayname','SF back');
% loglog(sfprof.ro{sfsubind}/rscale,sfprof.no{sfsubind}./vs,'g--','displayname','SF other');
% loglog(sfprof.rb{sfsubind}/rscale,sfprof.nb{sfsubind}./vs,'b--','displayname','SF back');
% plot([sfsize.rtidal(sfsubind)/rscale,sfsize.rtidal(sfsubind)/rscale],[1e-6,1e-1],':','displayname','SF tidal radius');
% loglog(sfprof.rs{sfsubind}/rscale,na./vs,'k-. o','displayname','SF other+back');
% hold off;
l=legend('show');
loglog([btsize.rtidal(btsubind)/rscale,btsize.rtidal(btsubind)/rscale],[1e-6,10],'k-');
loglog([sfsize.rtidal(sfsubind)/rscale,sfsize.rtidal(sfsubind)/rscale],[1e-6,10],'k--');
loglog([SofteningHalo/rscale,SofteningHalo/rscale],[1e-6,10],'k:');
set(l,'fontsize',15);
ylabel('$\rho (10^{10}M_{\odot}kpc^{\frac{\ }{\ }3}h^2)$');
% xlabel('$R/R_{tidal,HBT}$');
xlabel('$R /(\rm{kpc/h})$');
% set(gca,'xscale','linear');
% 
set(gcf,'paperposition',[0.6,6,20,17]);
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
% xlim([1e-3,1e1]);
% set(gca,'xtick',logspace(-3,1,5));
xlim([1e0,1e4]);
set(gca,'xtick',logspace(0,4,5));
set(gca,'xminortick','on');
ylim([1e-8,1e2]);
set(gca,'ytick',logspace(-8,2,6));
set(gca,'yminortick','on');
outputdir='/home/kam/Projects/HBT/code/data/show';
print('-depsc',[outputdir,'/satprof_',runnum,'_',num2str(btsubind),'_',num2str(sfsubind),'.new.eps']);
%% central profile
subind=7;
vs=4*pi/3*diff(logspace(log10(max(SofteningHalo,btmsize.rmax(subind)*1e-2)),log10(btmsize.rmax(subind)),btmsize.nbin(subind)+1).^3)';
svs=cumsum(vs);
na=btmprof.ns{subind}+btmprof.no{subind};%+btmprof.nb{subind};
sna=cumsum(na);
sns=cumsum(btmprof.ns{subind});
rscale=btmsize.rvir(subind);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
loglog(btmprof.rs{subind}/rscale,btmprof.ns{subind}./vs,'r-','displayname','HBT sub');
hold on;
loglog(btmprof.ro{subind}/rscale,btmprof.no{subind}./vs,'g-','displayname','HBT other');
loglog(btmprof.rb{subind}/rscale,btmprof.nb{subind}./vs,'b-','displayname','HBT back');
% loglog(btmprof.rs{subind}/rscale,na./vs,'k-. .','displayname','HBT back+sub');

vs=4*pi/3*diff(logspace(log10(max(SofteningHalo,sfmsize.rmax(subind)*1e-2)),log10(sfmsize.rmax(subind)),sfmsize.nbin(subind)+1).^3)';
svs=cumsum(vs);
na=sfmprof.ns{subind}+sfmprof.no{subind};%+sfmprof.nb{subind};
sna=cumsum(na);
sns=cumsum(sfmprof.ns{subind});
loglog(sfmprof.rs{subind}/rscale,sfmprof.ns{subind}./vs,'r--','displayname','SF sub');
loglog(sfmprof.ro{subind}/rscale,sfmprof.no{subind}./vs,'g--','displayname','SF other');
loglog(sfmprof.rb{subind}/rscale,sfmprof.nb{subind}./vs,'b--','displayname','SF back');
% loglog(sfmprof.rs{subind}/rscale,na./vs,'k-. o','displayname','SF other+sub');
hold off;
l=legend('show');
set(l,'fontsize',15);
ylabel('$\rho (10^{10}M_{\odot}kpc^{\frac{\ }{\ }3}h^2)$');
xlabel('$R/R_{vir,HBT}$');

set(gcf,'paperposition',[0.6,6,20,17]);
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/mainprof_',runnum,'_',num2str(subind-1),'.eps']);
%% sat profile of an overlapped subhalo
btsubind=621;

figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');

btsubind=627;
rscale=btsize.rtidal(btsubind);
vs=4*pi/3*diff(logspace(log10(max(SofteningHalo,btsize.rmax(btsubind)*1e-2)),log10(btsize.rmax(btsubind)),btsize.nbin(btsubind)+1).^3)';
svs=cumsum(vs);
na=btprof.no{btsubind}+btprof.ns{btsubind};
loglog(btprof.rs{btsubind}/rscale,btprof.ns{btsubind}./vs,'r:','displayname','small');
hold on;
loglog(btprof.ro{btsubind}/rscale,btprof.no{btsubind}./vs,'g-','displayname','big');
loglog(btprof.rs{btsubind}/rscale,na./vs,'k--','displayname','total');

sfsubind=554;
hold off;
l=legend('show');
set(l,'fontsize',15);
ylabel('$\rho (10^{10}M_{\odot}kpc^{\frac{\ }{\ }3}h^2)$');
xlabel('$R/R_{tidal,small}$');

set(gcf,'paperposition',[0.6,6,20,17]);
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/overlap.eps']);
%% an immersed subhalo
btsubind=727;  %id=7272
vs=4*pi/3*diff(logspace(log10(max(SofteningHalo,btsize.rmax(btsubind)*1e-2)),log10(btsize.rmax(btsubind)),btsize.nbin(btsubind)+1).^3)';
svs=cumsum(vs);
na=btprof.no{btsubind}+btprof.nb{btsubind};
rscale=btsize.rtidal(btsubind);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
loglog(btprof.rs{btsubind}/rscale,btprof.ns{btsubind}./vs,'r:','displayname','immersed sub');
hold on;
loglog(btprof.rs{btsubind}/rscale,na./vs,'g-','displayname','background');
l=legend('show');
set(l,'fontsize',15);
ylabel('$\rho (10^{10}M_{\odot}kpc^{\frac{\ }{\ }3}h^2)$');
xlabel('$R/R_{tidal}$');

set(gcf,'paperposition',[0.6,6,20,17]);
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/immerse.eps']);

%% sat profile of an overlapped subhalo
btsubind=1009;

figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
% rscale=btsize.rtidal(btsubind);
rscale=SofteningHalo;
vs=4*pi/3*diff(logspace(log10(max(SofteningHalo,btsize.rmax(btsubind)*1e-2)),log10(btsize.rmax(btsubind)),btsize.nbin(btsubind)+1).^3)';
svs=cumsum(vs);
% na=btprof.no{btsubind}+btprof.ns{btsubind};
loglog(btprof.rs{btsubind}/rscale,btprof.ns{btsubind}./vs,'k-','displayname','HBT sub');
hold on;
loglog(btprof.ro{btsubind}/rscale,btprof.no{btsubind}./vs,'b-','displayname','HBT other');
% loglog(btprof.rb{btsubind}/rscale,btprof.nb{btsubind}./vs,'b-','displayname','HBT back');

btsubind=1028;
% rscale=btsize.rtidal(btsubind);
vs=4*pi/3*diff(logspace(log10(max(SofteningHalo,btsize.rmax(btsubind)*1e-2)),log10(btsize.rmax(btsubind)),btsize.nbin(btsubind)+1).^3)';
svs=cumsum(vs);
na=btprof.no{btsubind}+btprof.ns{btsubind};
loglog(btprof.rs{btsubind}/rscale,btprof.ns{btsubind}./vs,'r--','displayname','HBT sub');
hold on;
loglog(btprof.ro{btsubind}/rscale,btprof.no{btsubind}./vs,'g--','displayname','HBT other');
% loglog(btprof.rb{btsubind}/rscale,btprof.nb{btsubind}./vs,'b-','displayname','HBT back');


sfsubind=2028;
hold off;
l=legend('show');
set(l,'fontsize',15);
ylabel('$\rho (10^{10}M_{\odot}kpc^{\frac{\ }{\ }3}h^2)$');
xlabel('$R/R_{tidal,small}$');

% set(gcf,'paperposition',[0.6,6,20,17]);
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
% print('-depsc',[outputdir,'/overlap.eps']);
