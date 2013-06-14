%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

markers=['s--';'o--';'d--';'^--'];
colors=['r';'g';'r';'g'];

virtype=0;
rmin=0;rmax=1;
nbin=15;
Nsnap=99;grpid=0;
runs={'6402','6402Low','6402','6402Low'};
dname={'6402 Unevo','6402Low Unevo','6402','6402Low'};
DMpmass=[0.008848,0.0104089];
RunName='6402ProCmp';
a=load_scaleF();
z=1/a(Nsnap+1)-1;

xref=logspace(10,13,5);
N=numel(runs);
M0=zeros(1,N);R0=M0;
xmass=cell(1,N);mfunspec=cell(1,N);mfunspecln=cell(1,N);mfuncum=cell(1,N);
for i=1:2
    basedir=['/mnt/A4700/data/',runs{i},'/subcat/anal/'];
    [Mlist,M0(i),R0(i)]=Mprolist_in_radii(basedir,virtype,Nsnap,grpid,rmin,rmax,'bindata','');
    Mlist=Mlist*DMpmass(i);
%     Mlist=[];
    MlistD=load_deathpro(basedir,Nsnap,grpid);
    Mlist=[Mlist;MlistD(MlistD>0)];
    [xmass{i},mfunspec{i},mfunspecln{i},mfuncum{i}]=mass_count(Mlist,nbin);
end
for i=3:4
    datadir=['/mnt/A4700/data/',runs{i},'/subcat/anal/massfun/'];
    [Mlist,M0(i),R0(i)]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','');
    [xmass{i},mfunspec{i},mfunspecln{i},mfuncum{i}]=mass_count(Mlist,nbin);
end
%% plot specific mass function
figure('visible','off');
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:N
    mfun=mfunspec{i}/M0(i);
    plot(log10(xmass{i}(:,2))+10,mfun(:,1),markers(i,:),...
        'color',colors(i),...
        'displayname',dname{i});
    hold on;
end
   
yref=10^-3.03*(xref).^-1.9*10^20;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/dM_{sub}/M_{host}\times(10^{10}M_{\odot}/h)^2$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);

print('-depsc',[outputdir,'/msfun_',RunName,'.eps']);
hgsave([outputdir,'/msfun_',RunName,'.fig']);
%% plot logspaced specific mass function
figure('visible','off');
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:N
    mfun=mfunspecln{i}/M0(i);
    plot(log10(xmass{i}(:,2))+10,mfun(:,1),markers(i,:),...
        'color',colors(i),...
        'displayname',dname{i});
    hold on;
end
   
yref=10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/d\ln M_{sub}/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);

print('-depsc',[outputdir,'/msfunln_',RunName,'.eps']);
hgsave([outputdir,'/msfunln_',RunName,'.fig']);
%% plot cumulative mass function
figure('visible','off');
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:N
    mfun=mfuncum{i}/M0(i);
    plot(log10(xmass{i}(:,1))+10,mfun(:,1),markers(i,:),...
        'color',colors(i),...
        'displayname',dname{i});
    hold on;
end

    
yref=1/0.9*10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','Giocoli09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$N(>M_{sub})/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);

print('-depsc',[outputdir,'/msfuncum_',RunName,'.eps']);
hgsave([outputdir,'/msfuncum_',RunName,'.fig']);
