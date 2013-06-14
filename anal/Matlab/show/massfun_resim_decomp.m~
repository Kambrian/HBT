%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

markers=['s--';'o--';'d--';'^--';'<--';'>--';'p--';'x--';'+--';'.--'];
colors=['r';'g';'b';'c';'m';'k';'r';'g';'b';'c'];

virtype=0;
rmin=0;rmax=1;
xmin=10^-0.6;xmax=10^3;
nbin=15;
SnapStart=39;SnapEnd=40;
grpid=0;
% runs={'6702','6702DM'};
runs='6702DM';
RunName=[runs,'DecompS',num2str(SnapStart),'_S',num2str(SnapEnd),'V',num2str(virtype)];
a=load_scaleF();
z1=1./a(SnapStart+1)-1;
z2=1./a(SnapEnd+1)-1;
% z=1./a(Nsnap+1)-1;

xref=logspace(10,13,5);
M0=zeros(2,1);R0=M0;
xmass=cell(2,1);mfunspec=cell(2,1);mfunspecln=cell(2,1);mfuncum=cell(2,1);
datadir=['/mnt/A4700/data/',runs,'/subcat/anal/massfun/'];
[Mlist,M0(1),R0(1)]=Msublist_in_radii(datadir,virtype,SnapEnd,grpid,rmin,rmax,'bindata','');
[xmass{1},mfunspec{1},mfunspecln{1},mfuncum{1}]=mass_count(Mlist,nbin,xmin,xmax);
datadir=['/mnt/A4700/data/',runs,'/subcat/anal/massfun/rem_',num2str(SnapStart,'%03d'),'/'];
[Mlist,M0(2),R0(2)]=Msublist_in_radii(datadir,virtype,SnapEnd,grpid,rmin,rmax,'bindata','');
[xmass{2},mfunspec{2},mfunspecln{2},mfuncum{2}]=mass_count(Mlist,nbin,xmin,xmax);
%% plot specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    
    mfun=mfunspec{1}/M0(1)-mfunspec{2}/M0(2);
    plot(log10(xmass{1}(:,2))+10,mfun(:,1),markers(1,:),...
        'color',colors(1),...
        'displayname','new accretion');
    hold on;
    
    mfun=mfunspec{2}/M0(2);
    plot(log10(xmass{2}(:,2))+10,mfun(:,1),markers(2,:),...
        'color',colors(2),...
        'displayname','remnant');
    hold on;

  
yref=10^-3.03*(xref).^-1.9*10^20;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/dM_{sub}/M_{host}\times(10^{10}M_{\odot}/h)^2$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title([runs,' z=',num2str(z1,'%2.1f'),'~',num2str(z2,'%2.1f')]);
 
% print('-depsc',[outputdir,'/msfun_',RunName,'.eps']);
% hgsave([outputdir,'/msfun_',RunName,'.fig']);
%% plot logspaced specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfunspecln{1}/M0(1)-mfunspecln{2}/M0(2);
    plot(log10(xmass{1}(:,2))+10,mfun(:,1),markers(1,:),...
        'color',colors(1),...
        'displayname','new accretion');
    hold on;
    
    mfun=mfunspecln{2}/M0(2);
    plot(log10(xmass{2}(:,2))+10,mfun(:,1),markers(2,:),...
        'color',colors(2),...
        'displayname','remnant');
    hold on;


    
yref=10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/d\ln M_{sub}/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title([runs,' z=',num2str(z1,'%2.1f'),'~',num2str(z2,'%2.1f')]);

% print('-depsc',[outputdir,'/msfunln_',RunName,'.eps']);
% hgsave([outputdir,'/msfunln_',RunName,'.fig']);
%% plot cumulative mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfuncum{1}/M0(1)-mfuncum{2}/M0(2);
    plot(log10(xmass{1}(:,2))+10,mfun(:,1),markers(1,:),...
        'color',colors(1),...
        'displayname','new accretion');
    hold on;
    
    mfun=mfuncum{2}/M0(2);
    plot(log10(xmass{2}(:,2))+10,mfun(:,1),markers(2,:),...
        'color',colors(2),...
        'displayname','remnant');
    hold on;
    
yref=1/0.9*10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$N(>M_{sub})/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title([runs,' z=',num2str(z1,'%2.1f'),'~',num2str(z2,'%2.1f')]);
% 
% print('-depsc',[outputdir,'/msfuncum_',RunName,'.eps']);
% hgsave([outputdir,'/msfuncum_',RunName,'.fig']);
