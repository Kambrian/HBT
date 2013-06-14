%% data preparation
% to add: 6702low, cosmo

% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
outputdir='/work/Projects/HBT/code/data/show/massfun';
addpath(genpath('../post'));

rmin=0;rmax=1;
nbin=10;
Nsnap=59;
RunNum=8213;

datadir=['/mnt/A4700/data/',num2str(RunNum),'/subcat/anal/massfun/'];
submass=read_submass([datadir,'submass_',num2str(Nsnap,'%03d')]);
subcom=read_subcom([datadir,'subcom_',num2str(Nsnap,'%03d')]);
massdata=[submass,subcom];
clear submass subcom
cenid=read_cid([datadir,'cid_',num2str(Nsnap,'%03d')]);
grpsize=read_grpsize([datadir,'grpsizeVIR_',num2str(Nsnap,'%03d')]);
Mhost=grpsize(:,1);
Rvir=grpsize(:,2);
clear grpsize
%%
gid=100;
M0=Mhost(gid);
HMlist=Msublists_cosmo(massdata,cenid(gid)+1,rmin*Rvir(gid),rmax*Rvir(gid));
Hlen=numel(HMlist);
% HMlist=HMlist/M0;
[xmass,mfunspec,mfunspecln,mfuncum]=mass_count(HMlist,nbin);
%%% plot specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',4);

%     mfun=mfunspec;
%     errorbar(log10(xmass(:,2)),mfun(:,1),mfun(:,2),'c^--',...
%         'displayname',[num2str(M0*1e10,'%2.1e'),'$M_{\odot}/h$']);
%     hold on;
% xref=logspace(-5,-1,5);
% yref=10^-3.02*(M0*1e10)^0.1*(xref).^-1.9;

    mfun=mfunspec/M0;
    errorbar(log10(xmass(:,2))+10,mfun(:,1),mfun(:,2),'c^--',...
        'displayname',[num2str(M0*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
xref=logspace(8,14,5);
yref=10^-3.03*(xref).^-1.9*10^20;
plot(log10(xref),yref,'-k','displayname','G09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/(M_{\odot}/h)$)','interpreter','latex');
ylabel('$dN/dM_{sub}/M_{host}\times(10^{10}M_{\odot}/h)^2$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');

% print('-depsc',[outputdir,'/msfun_cosm_bin.eps']);
% hgsave([outputdir,'/msfun_cosm_bin.fig']);
%% plot logspaced specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',4);

    mfun=mfunspecln/M0;
    errorbar(log10(xmass(:,2))+10,mfun(:,1),mfun(:,2),'c^--',...
        'displayname',[num2str(M0*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;

    
xref=logspace(8,14,5);
yref=10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','G09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h)$)','interpreter','latex');
ylabel('$dN/d\ln M_{sub}/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');

% print('-depsc',[outputdir,'/msfunln_cosm_bin.eps']);
% hgsave([outputdir,'/msfunln_cosm_bin.fig']);
%% plot cumulative mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',4);

    mfun=mfuncum/M0;
    errorbar(log10(xmass(:,1))+10,mfun(:,1),mfun(:,2),'c^--',...
        'displayname',[num2str(M0*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;

    
xref=logspace(8,14,5);
yref=1/0.9*10^-3.03*(xref).^-0.9*10^10;
plot(log10(xref),yref,'-k','displayname','G09');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h)$)','interpreter','latex');
ylabel('$N(>M_{sub})/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');

% print('-depsc',[outputdir,'/msfuncum_cosm_bin.eps']);
% hgsave([outputdir,'/msfuncum_cosm_bin.fig']);
