%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

markers=['s--';'o--';'d--';'^--'];
colors=['r';'g';'b';'c'];

virtype=0;
rmin=0;rmax=1;
nbin=20;
Nsnap=99;grpid=0;
runs={'6702','6702NSM'};
RunName='6702NSM';
a=load_scaleF();
z=1/a(Nsnap+1)-1;

xref=logspace(-4,-2,5);
N=numel(runs);
M0=zeros(1,N);R0=M0;
xmass=cell(1,N);mfunspec=cell(1,N);mfunspecln=cell(1,N);mfuncum=cell(1,N);
for i=1:N
    datadir=['/mnt/A4700/data/',runs{i},'/subcat/anal/massfun/'];
    [Mlist,M0(i),R0(i)]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','norm');
    [xmass{i},mfunspec{i},mfunspecln{i},mfuncum{i}]=mass_count(Mlist,nbin);
end
%% subfind data for 6702
% subfinddir='/home/kam/Documents/research/Galaxy/code/BoundTracing/v7.1/data/massfun/subfind';
% load([subfinddir,'/SubFind_mass_099']);
% load([subfinddir,'/SubFindgroup_offset_099']);
% G=43007.1;
% Hz=0.1;Omega=0.3;
% partmass=0.008848;
% % MS=SubFindgroup_offset_099(1,3);
% % Rvir=SubFindgroup_offset_099(1,4);%Critical_200
% Rsub_s=sqrt(sum((SubFind_mass_099(:,3:5)-repmat(SubFind_mass_099(1,3:5),size(SubFind_mass_099,1),1)).^2,2));
% id=1:size(SubFind_mass_099,1);id=id';
% Mlist_S=SubFind_mass_099(find(Rsub_s>rmin*R0(1)&Rsub_s<rmax*R0(1)&id~=1),1)*partmass;
% [xmassS,mfunspecS,mfunspeclnS,mfuncumS]=mass_count(Mlist_S,nbin);
%% plot specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:N
    mfun=mfunspec{i};
    plot(log10(xmass{i}(:,2)),mfun(:,1),markers(i,:),...
        'color',colors(i),...
        'displayname',runs{i});
    hold on;
end

%     mfun=mfunspecS/M0(1);
%     plot(log10(xmassS(:,2))+10,mfun(:,1),'mo-',...
%         'displayname','SubFind');
 
yref=10^-1*(xref).^-1.9;
plot(log10(xref),yref,'-k','displayname','slope=-0.9');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d(M_{sub}/M_{vir})$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfunN_',RunName,'.eps']);
hgsave([outputdir,'/msfunN_',RunName,'.fig']);
%% plot logspaced specific mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:N
    mfun=mfunspecln{i};
    plot(log10(xmass{i}(:,2)),mfun(:,1),markers(i,:),...
        'color',colors(i),...
        'displayname',runs{i});
    hold on;
end

%     mfun=mfunspeclnS/M0(1);
%     plot(log10(xmassS(:,2))+10,mfun(:,1),'mo-',...
%         'displayname','SubFind');
    
yref=10^-1*(xref).^-0.9;
plot(log10(xref),yref,'-k','displayname','slope=-0.9');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d\ln(M_{sub}/M_{vir})$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfunNln_',RunName,'.eps']);
hgsave([outputdir,'/msfunNln_',RunName,'.fig']);
%% plot cumulative mass function
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

for i=1:N
    mfun=mfuncum{i};
    plot(log10(xmass{i}(:,1)),mfun(:,1),markers(i,:),...
        'color',colors(i),...
        'displayname',runs{i});
    hold on;
end

%     mfun=mfuncumS/M0(1);
%     plot(log10(xmassS(:,1))+10,mfun(:,1)'mo-',...
%         'displayname','SubFind');
    
yref=1/0.9*10^-1*(xref).^-0.9;
plot(log10(xref),yref,'-k','displayname','slope=-0.9');hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$N(>M_{sub}/M_{vir})$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
print('-depsc',[outputdir,'/msfunNcum_',RunName,'.eps']);
hgsave([outputdir,'/msfunNcum_',RunName,'.fig']);
