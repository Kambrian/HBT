%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

virtype=0;
rmin=0;rmax=1;
nbin=9;
Nsnap=99;grpid=0;
RunName='6702DM';
a=load_scaleF();
z=1/a(Nsnap+1)-1;

xref=logspace(-6,-1,5);
datadir=['/mnt/A4700/data/',RunName,'/subcat/anal/massfun/'];
[Mlist,M0,R0]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','norm');
xmin=min(Mlist);xmax=max(Mlist)*1.001;
[xmass,mfunspec,mfunspecln,mfuncum]=mass_count(Mlist,nbin,xmin,xmax);
datadir=['/mnt/A4700/data/',RunName,'/subcatS/anal/massfun/'];
[Mlist,M0,R0]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','norm');
[xmass2,mfunspec2,mfunspecln2,mfuncum2]=mass_count(Mlist,nbin,xmin,xmax);
%% plot specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfunspec;
    errorbar(log10(xmass(:,2)),mfun(:,1),mfun(:,2),'c^-',...
        'displayname','HBT');
    hold on;

    mfun=mfunspec2;
    errorbar(log10(xmass2(:,2)),mfun(:,1),mfun(:,2),'mo-',...
        'displayname','SubFind');
    
yref=10^-3.02*(M0*1e10)^0.1*(xref).^-1.9;
% yref=10^-1*(xref).^-1.9;
plot(log10(xref),yref,'-k','displayname','G10');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d(M_{sub}/M_{vir})$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfunN_',RunName,'.eps']);
% hgsave([outputdir,'/msfunN_',RunName,'.fig']);
%% plot logspaced specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

axes('position',[0.15,0.45,0.8,0.5]);

    mfun=mfunspecln;
    errorbar(log10(xmass(:,2)),mfun(:,1),mfun(:,2),'c^-',...
        'displayname','HBT');
    hold on;

    mfun=mfunspecln2;
    errorbar(log10(xmass2(:,2)),mfun(:,1),mfun(:,2),'mo-',...
        'displayname','SubFind');
    
yref=10^-3.02*(M0*1e10)^0.1*(xref).^-0.9;
plot(log10(xref),yref,'-k','displayname','G10');hold off;
set(gca,'yscale','log','yminortick','on');
set(gca,'xlim',[-6,-0.5],'xticklabel',[]);
% xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d\ln(M_{sub}/M_{vir})$','interpreter','latex');
hl=legend('show','location','northeast');set(hl,'interpreter','latex');
% title(['z=',num2str(z,'%2.1f')]);

axes('position',[0.15,0.15,0.8,0.3]);
mfun=mfunspecln./mfunspecln2;
plot(log10(xmass2(:,2)),mfun(:,1),'-','displayname','data');
hold on;
% plot([-6,-0.5],[1,1],':');
ymean=mean(mfun(1:end-3,1));
plot([-6,-0.5],[ymean,ymean],':','displayname',['y=',num2str(ymean,'%3.2f')]);
set(gca,'xlim',[-6,-0.5],'ytick',[0,1,2]);
ystr={'$dN_{HBT}/dN_{SF}$','\ '};
ylabel(ystr,'interpreter','latex');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
h2=legend('show','location','northeast');set(h2,'interpreter','latex');


print('-depsc',[outputdir,'/msfunNln_',RunName,'.eps']);
% hgsave([outputdir,'/msfunNln_',RunName,'.fig']);
%% plot cumulative mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfuncum;
    errorbar(log10(xmass(:,1)),mfun(:,1),mfun(:,2),'c^-',...
        'displayname','HBT');
    hold on;

    mfun=mfuncum2;
    errorbar(log10(xmass2(:,1)),mfun(:,1),mfun(:,2),'mo-',...
        'displayname','SubFind');
    
yref=1/0.9*10^-3.02*(M0*1e10)^0.1*(xref).^-0.9;
plot(log10(xref),yref,'-k','displayname','G10');hold off;
set(gca,'yscale','log','yminortick','on');
ylabel('$N(>M_{sub}/M_{vir})$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);

xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');

print('-depsc',[outputdir,'/msfunNcum_',RunName,'.eps']);
% hgsave([outputdir,'/msfunNcum_',RunName,'.fig']);