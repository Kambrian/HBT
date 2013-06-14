%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

markers=['o--';'s--';'d--';'^--'];
colors=['r';'g';'b';'c'];

virtype=0;
rmin=0;rmax=1;
nbin=10;
Nsnap=99;grpid=0;
runs=[6402,6404,6409;6500,6501,6506;6600,6601,6602;6700,6701,6702];
a=load_scaleF();
z=1/a(Nsnap+1)-1;

M0=zeros(size(runs));
xmass=cell(1,4);mfunspec=cell(1,4);mfunspecln=cell(1,4);mfuncum=cell(1,4);
for i=1:4
    Mlist=[];
    for j=1:3
        datadir=['/mnt/A4700/data/',num2str(runs(i,j)),'/subcat/anal/massfun/'];
        [Msubs,M0(i,j)]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','norm');
        Mlist=[Mlist;Msubs];
    end
    [xmass{i},mfunspec{i},mfunspecln{i},mfuncum{i}]=mass_count(Mlist,nbin);
end
% %% plot specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:4
    mfun=mfunspec{i}/3;
    errorbar(log10(xmass{i}(:,2)),mfun(:,1),mfun(:,2),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(mean(M0(i,:))*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
end
    
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d(M_{sub}/M_{vir})/N_{host}$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
% print('-depsc',[outputdir,'/msfunN_resim_bin.eps']);
% hgsave([outputdir,'/msfunN_resim_bin.fig']);
%% plot logspaced specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:4
    mfun=mfunspecln{i}/3;
    errorbar(log10(xmass{i}(:,2)),mfun(:,1),mfun(:,2),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(mean(M0(i,:))*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
end
    
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d\ln(M_{sub}/M_{vir})/N_{host}$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfunNln_resim_bin.eps']);
hgsave([outputdir,'/msfunNln_resim_bin.fig']);
%% plot cumulative mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:4
    mfun=mfuncum{i}/3;
    errorbar(log10(xmass{i}(:,1)),mfun(:,1),mfun(:,2),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(mean(M0(i,:))*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
end
    
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$N(>M_{sub}/M_{vir})/N_{host}$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfunNcum_resim_bin.eps']);
hgsave([outputdir,'/msfunNcum_resim_bin.fig']);
