%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim/';
addpath(genpath('../post'));

markers=['o--';'s--';'d--';'^--';'<--'];
colors=['r';'g';'b';'c';'m'];

virtype=0;
rmin=0;rmax=1;
nbin=10;
Nsnap=99;grpid=0;
runs=[6402,6404,6409;6500,6501,6506;6600,6601,6602;6700,6701,6702];
% runs=[6402;6500;6600;6700;];
% runs=[6404;6501;6601;6701;];
% runs=[6409;6506;6602;6702;];
%  only 6500 and 67** seems acceptable, others bend down severely in low
%  mass of their subhalo massfunction
a=load_scaleF();
z=1/a(Nsnap+1)-1;

M0=zeros(size(runs));
xmass=cell(1,size(runs,1));mfunspec=cell(1,size(runs,1));mfunspecln=cell(1,size(runs,1));mfuncum=cell(1,size(runs,1));
for i=1:size(runs,1);
    Mlist=[];
    for j=1:size(runs,2)
        datadir=['/mnt/A4700/data/',num2str(runs(i,j)),'/subcat/anal/massfun/'];
        [Msubs,M0(i,j)]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','');
        Mlist=[Mlist;Msubs];
    end
    [xmass{i},mfunspec{i},mfunspecln{i},mfuncum{i}]=mass_count(Mlist,nbin);
end
%% plot specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:size(runs,1)
    mfun=mfunspec{i}/sum(M0(i,:));
    errorbar(log10(xmass{i}(:,2))+10,mfun(:,1),mfun(:,2),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(mean(M0(i,:))*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
end
    
xref=logspace(8,14,5);
yref=10^-3.03*(xref).^-1.9*10^20;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/dM_{sub}/M_{host}\times(10^{10}M_{\odot}/h)^2$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);

fname=['msfun_resim_binS',num2str(Nsnap)];
print('-depsc',[outputdir,fname,'.eps']);
hgsave([outputdir,fname,'.fig']);
%% plot logspaced specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:size(runs,1)
    mfun=mfunspecln{i}/sum(M0(i,:));
    errorbar(log10(xmass{i}(:,2))+10,mfun(:,1),mfun(:,2),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(mean(M0(i,:))*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
end
    
xref=logspace(8,14,5);
yref=10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/d\ln M_{sub}/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);

fname=['msfunln_resim_binS',num2str(Nsnap)];
print('-depsc',[outputdir,fname,'.eps']);
hgsave([outputdir,fname,'.fig']);
%% plot cumulative mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:size(runs,1)
    mfun=mfuncum{i}/sum(M0(i,:));
    errorbar(log10(xmass{i}(:,1))+10,mfun(:,1),mfun(:,2),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(mean(M0(i,:))*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
end
    
xref=logspace(8,14,5);
yref=1/0.9*10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$N(>M_{sub})/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);

fname=['msfuncum_resim_binS',num2str(Nsnap)];
print('-depsc',[outputdir,fname,'.eps']);
hgsave([outputdir,fname,'.fig']);
