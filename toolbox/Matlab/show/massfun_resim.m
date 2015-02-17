%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Projects/HBT/code/data/show/massfun/resim';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

virtype=0;
rmin=0;rmax=1;
nbin=8;
Nsnap=162;grpid=0;
RunName='AqA2';
% scaleF_file=['/mnt/charon/HBT/data/',RunName,'/subcat/Redshift.dat'];
% a=load_scaleF(scaleF_file);
% z=1/a(Nsnap+1)-1;
z=0;

% datadir=['/mnt/A4700/data/',RunName,'/subcat/anal/massfun/'];
datadir=['/mnt/charon/HBT/data/',RunName,'/subcatmore/anal/massfun/'];
[Mlist,M0,R0]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','');M0,R0
[xmass,mfunspec,mfunspecln,mfuncum]=mass_count(Mlist,nbin);
% datadir=['/mnt/A4700/data/',RunName,'/subcatS/anal/massfun/'];
datadir=['/mnt/charon/HBT/data/',RunName,'/subcatmore/anal/subfind/massfun/'];
[Mlist,M0,R0]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','');M0,R0
[xmass2,mfunspec2,mfunspecln2,mfuncum2]=mass_count(Mlist,nbin);
%% plot specific mass function
figure;%('visible','off');
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfunspec/M0;
    errorbar(log10(xmass(:,2))+10,mfun(:,1),mfun(:,2),'c^-',...
        'displayname','HBT');
    hold on;
    
    mfun=mfunspec2/M0;
    errorbar(log10(xmass2(:,2))+10,mfun(:,1),mfun(:,2),'mo-',...
        'displayname','SUBFIND');
    hold on;
    
xref=logspace(10,13,5);
yref=10^-3.03*(xref).^-1.9*10^20;
plot(log10(xref),yref,'-k','displayname','G10');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/dM_{sub}/M_{host}\times(10^{10}M_{\odot}/h)^2$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfun_',RunName,'.eps']);
% hgsave([outputdir,'/msfun_',RunName,'.fig']);
%% plot logspaced specific mass function
figure;%('visible','off');
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfunspecln/M0;
    errorbar(log10(xmass(:,2))+10,mfun(:,1),mfun(:,2),'c^-',...
        'displayname','HBT');
    hold on;

    mfun=mfunspecln2/M0;
    errorbar(log10(xmass2(:,2))+10,mfun(:,1),mfun(:,2),'mo-',...
        'displayname','SubFind');
    
xref=logspace(5,10,5);
yref=10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','G10');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/d\ln M_{sub}/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfunln_',RunName,'.eps']);
% hgsave([outputdir,'/msfunln_',RunName,'.fig']);
%% plot cumulative mass function
figure;%('visible','off');
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfuncum;
    errorbar(log10(xmass(:,1))+10,mfun(:,1),mfun(:,2),'c^-',...
        'displayname','HBT');
    hold on;
     
    mfun=mfuncum2;
    errorbar(log10(xmass2(:,1))+10,mfun(:,1),mfun(:,2),'mo-',...
        'displayname','SUBFIND');
    hold on;
    
xref=logspace(6,10,5);
yref=1/0.9*10^-3.03*(xref).^-0.9*10^10;
% plot(log10(xref),yref,'-k','displayname','G10');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$N(>M_{sub})$','interpreter','latex');
% ylabel('$N(>M_{sub})/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfuncum_',RunName,'.eps']);
hgsave([outputdir,'/msfuncum_',RunName,'.fig']);
