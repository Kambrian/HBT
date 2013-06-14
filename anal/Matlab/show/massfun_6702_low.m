%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

markers=['s--';'o--';'d--';'^--'];
colors=['r';'g';'b';'c'];

virtype=0;
rmin=0;rmax=1;
nbin=15;
nbinl=10;
Nsnap=99;grpid=0;
RunName='6702Low';

datadir=['/mnt/A4700/data/6702/subcat/anal/massfun/'];
[Mlist,M0,R0]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','');
[xmass,mfunspec,mfunspecln,mfuncum]=mass_count(Mlist,nbin);
%% subfind data for 6702
subfinddir='~/BT/data/6702LowMSFun';
load([subfinddir,'/SubFind_mass_099']);
load([subfinddir,'/SubFindgroup_offset_099']);
partmass=0.008848;
Rsub_s=sqrt(sum((SubFind_mass_099(:,3:5)-repmat(SubFind_mass_099(1,3:5),size(SubFind_mass_099,1),1)).^2,2));
id=1:size(SubFind_mass_099,1);id=id';
Mlist_S=SubFind_mass_099(find(Rsub_s>rmin*R0&Rsub_s<rmax*R0&id~=1),1)*partmass;
[xmassS,mfunspecS,mfunspeclnS,mfuncumS]=mass_count(Mlist_S,nbin);
%%
load([subfinddir,'/SubFind_mass_099_l']);
load([subfinddir,'/SubFindgroup_offset_099_l']);
partmass=0.897207;
Rsub_sL=sqrt(sum((SubFind_mass_099_l(:,3:5)-repmat(SubFind_mass_099_l(1,3:5),size(SubFind_mass_099_l,1),1)).^2,2));
id=1:size(SubFind_mass_099_l,1);id=id';
Mlist_SL=SubFind_mass_099_l(find(Rsub_sL>rmin*R0&Rsub_sL<rmax*R0&id~=1),1)*partmass;
[xmassSL,mfunspecSL,mfunspeclnSL,mfuncumSL]=mass_count(Mlist_SL,nbinl);
%%
load([subfinddir,'/sub_mass_099_l']);
load([subfinddir,'/group_offset_099_l']);
G=43007.1;
Hz=0.1;Omega=0.3;
partmass=0.897207;
Rsub_l=sqrt(sum((sub_mass_099_l(:,2:4)-repmat(sub_mass_099_l(1,2:4),size(sub_mass_099_l,1),1)).^2,2));
id=1:size(sub_mass_099_l,1);id=id';
Mlist_l=sub_mass_099_l(find(Rsub_l>rmin*R0&Rsub_l<rmax*R0&id~=1),1);
[xmassL,mfunspecL,mfunspeclnL,mfuncumL]=mass_count(Mlist_l,nbinl);
%% plot specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfunspec/M0;
    plot(log10(xmass(:,2))+10,mfun(:,1),'k-',...
        'displayname','6702');
    hold on;
    
    mfun=mfunspecL/M0;
    plot(log10(xmassL(:,2))+10,mfun(:,1),'k^',...
        'displayname','6702Low');
    hold on;


    mfun=mfunspecS/M0;
    plot(log10(xmassS(:,2))+10,mfun(:,1),'m--',...
        'displayname','SubFind');
    
    mfun=mfunspecSL/M0;
    plot(log10(xmassSL(:,2))+10,mfun(:,1),'mo',...
        'displayname','SubFindLow');
    
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/dM_{sub}/M_{host}\times(10^{10}M_{\odot}/h)^2$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');

print('-depsc',[outputdir,'/msfun_',RunName,'.eps']);
hgsave([outputdir,'/msfun_',RunName,'.fig']);
%% plot logspaced specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfunspecln/M0;
    plot(log10(xmass(:,2))+10,mfun(:,1),'k-',...
        'displayname','6702');
    hold on;
    
    mfun=mfunspeclnL/M0;
    plot(log10(xmassL(:,2))+10,mfun(:,1),'k^',...
        'displayname','6702Low');
    hold on;


    mfun=mfunspeclnS/M0;
    plot(log10(xmassS(:,2))+10,mfun(:,1),'m-',...
        'displayname','SubFind');
    
    mfun=mfunspeclnSL/M0;
    plot(log10(xmassSL(:,2))+10,mfun(:,1),'mo',...
        'displayname','SubFindLow');
    
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dN/d\ln M_{sub}/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');

print('-depsc',[outputdir,'/msfunln_',RunName,'.eps']);
hgsave([outputdir,'/msfunln_',RunName,'.fig']);
%% plot cumulative mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

    mfun=mfuncum/M0;
    plot(log10(xmass(:,2))+10,mfun(:,1),'k-',...
        'displayname','6702');
    hold on;
    
    mfun=mfuncumL/M0;
    plot(log10(xmassL(:,2))+10,mfun(:,1),'k^',...
        'displayname','6702Low');
    hold on;


    mfun=mfuncumS/M0;
    plot(log10(xmassS(:,2))+10,mfun(:,1),'m-',...
        'displayname','SubFind');
    
    mfun=mfuncumSL/M0;
    plot(log10(xmassSL(:,2))+10,mfun(:,1),'mo',...
        'displayname','SubFindLow');
    
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$N(>M_{sub})/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');

print('-depsc',[outputdir,'/msfuncum_',RunName,'.eps']);
hgsave([outputdir,'/msfuncum_',RunName,'.fig']);
