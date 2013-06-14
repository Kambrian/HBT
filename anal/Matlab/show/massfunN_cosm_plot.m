%% load data

addpath(genpath('../post'));

markers=['o--';'s--';'d--';'^--';'<--';'x--'];
% markers=['.-';'.-';'.-';'.-';'.-';'.-'];
colors=['r';'g';'b';'c';'m';'k'];

Nsnap=49;RunNum='BJLGR';name='msfunN';
virtype=0;

outputdir=['/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/',num2str(RunNum),'/'];
datadir=['/mnt/A4700/data/',num2str(RunNum),'/subcat/anal/massfun/'];
[data,redshift]=read_massfun([datadir,'massfunN_',num2str(Nsnap,'%03d'),'.',num2str(virtype,'%d')]);
nfun=numel(data);

Mhost=zeros(nfun,1);
for i=1:nfun
    Mhost(i)=data(i).Nhost;
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

for i=1:nfun
    mfun=data(i).mfunspec/Mhost(i);
    xmass=data(i).xmass(:,2);
    plot(log10(xmass),mfun(:,1),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(log10(data(i).Mbin(1))+10,'%2.1f'),'$\sim$',num2str(log10(data(i).Mbin(2))+10,'%2.1f')]);
    hold on;
end

hold off;
set(gca,'yscale','log','xminortick','on','yminortick','on');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d(M_{sub}/M_{vir})/N_{host}$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(redshift,'%2.1f')]);
% 
fname=[name,'_',num2str(RunNum),'S',num2str(Nsnap),'V',num2str(virtype)];
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


for i=1:nfun
    mfun=data(i).mfunspecln/Mhost(i);
%     mfun=data(i).mfunspec.*repmat(data(i).xmass(:,2),1,2)/Mhost(i);
    xmass=data(i).xmass(:,2);
    plot(log10(xmass),mfun(:,1),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(log10(data(i).Mbin(1))+10,'%2.1f'),'$\sim$',num2str(log10(data(i).Mbin(2))+10,'%2.1f')]);
%     errorbar(log10(xmass)+10,mfun(:,1),mfun(:,2),markers(i,:),...
%         'color',colors(i),'markerfacecolor',colors(i),...
%         'displayname',[num2str(data(i).Mhost/data(i).Nhost*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
end  

hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d\ln(M_{sub}/M_{vir})/N_{host}$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(redshift,'%2.1f')]);

fname=[name,'ln_',num2str(RunNum),'S',num2str(Nsnap),'V',num2str(virtype)];
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


for i=1:nfun
    mfun=data(i).mfuncum/Mhost(i);
    xmass=data(i).xmass(:,1);
    plot(log10(xmass),mfun(:,1),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(log10(data(i).Mbin(1))+10,'%2.1f'),'$\sim$',num2str(log10(data(i).Mbin(2))+10,'%2.1f')]);
%     errorbar(log10(xmass)+10,mfun(:,1),mfun(:,2),markers(i,:),...
%         'color',colors(i),'markerfacecolor',colors(i),...
%         'displayname',[num2str(data(i).Mhost/data(i).Nhost*1e10,'%2.1e'),'$M_{\odot}/h$']);
    hold on;
end
   
hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$N(>M_{sub}/M_{vir})/N_{host}$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(redshift,'%2.1f')]);

fname=[name,'cum_',num2str(RunNum),'S',num2str(Nsnap),'V',num2str(virtype)];
print('-depsc',[outputdir,fname,'.eps']);
hgsave([outputdir,fname,'.fig']);
