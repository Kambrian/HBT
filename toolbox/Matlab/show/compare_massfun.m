%% load data
addpath(genpath('/work/Projects/Lensing/code/v8.0/Matlab'));
addpath(genpath('../post'));
outputdir=['/work/Projects/HBT/code/data/show/massfun/'];
colors=['r';'g';'b';'c';'m';'k'];
alpha=1; %multiply by power alpha
xref=logspace(9,14.,5);
%%
myfigure;
% yref=10^-3.2*(xref).^-0.9*10^10;plot(log10(xref),yref.*(xref/1e10).^alpha,'-k','linewidth',2,'displayname','Gao04');hold on;
% yref2=10^-3.03*(xref).^-0.9*10^10;plot(log10(xref),yref2.*(xref/1e10).^alpha,'-k','linewidth',2,'displayname','G10,z=0');hold on;
yref3=10^-4.15*(xref).^-0.6*10^10/100/1.;plot(log10(xref),yref3.*(xref/1e10).^alpha,'--k','linewidth',2,'displayname','$\alpha$=-0.6');hold on;
yref3=10^-3.*(xref).^-0.9*10^10/1.;plot(log10(xref),yref3.*(xref/1e10).^alpha,'-k','linewidth',2,'displayname','$\alpha$=-0.9');hold on;

virtype=0;

Nsnap=99; RunNum='6116';
datadir=['/mnt/A4700/data/',RunNum,'/Analysis/'];
[data,redshift]=read_massfun([datadir,'massfun_',num2str(Nsnap,'%03d'),'.',num2str(virtype,'%d')]);%,'.R1.0-2.0']);
nfun=numel(data);
Mhost=zeros(nfun,1);
for i=1:nfun
    Mhost(i)=data(i).Mhost;
end
for i=1:2:nfun
    mfun=data(i).mfunspecln/Mhost(i);
    xmass=data(i).xmass(:,2);
    plot(log10(xmass)+10,mfun(:,1).*(xmass).^alpha,'-',...
        'color',colors(i), 'linewidth', 5,...
        'displayname',[num2str(log10(data(i).Mbin(1))+10,'%2.1f'),'$\sim$',num2str(log10(data(i).Mbin(2))+10,'%2.1f')]);
end
hl=legend('show','location','southwest');set(hl,'interpreter','latex');

% Nsnap=99; RunNum='8511';
% datadir=['/mnt/uv/',RunNum,'/subcat/anal/massfun/'];
% [data,redshift]=read_massfun([datadir,'massfun_',num2str(Nsnap,'%03d'),'.',num2str(virtype,'%d')]);%,'.R1.0-2.0']);
% nfun=numel(data);
% Mhost=zeros(nfun,1);
% for i=1:2:nfun
%     Mhost(i)=data(i).Mhost;
% end
% for i=1:nfun
%     mfun=data(i).mfunspecln/Mhost(i);
%     xmass=data(i).xmass(:,2);
%     plot(log10(xmass)+10,mfun(:,1).*(xmass).^alpha,'-',...
%         'color',colors(i), 'linewidth', 2,...
%         'displayname',[num2str(log10(data(i).Mbin(1))+10,'%2.1f'),'$\sim$',num2str(log10(data(i).Mbin(2))+10,'%2.1f')]);
% end
% 
% Nsnap=59; RunNum='8213';
% datadir=['/mnt/A4700/data/',RunNum,'/Analysis/'];
% [data,redshift]=read_massfun([datadir,'massfun_',num2str(Nsnap,'%03d'),'.',num2str(virtype,'%d')]);%,'.R1.0-2.0']);
% nfun=numel(data);
% Mhost=zeros(nfun,1);
% for i=1:nfun
%     Mhost(i)=data(i).Mhost;
% end
% for i=1:2:nfun
%     mfun=data(i).mfunspecln/Mhost(i);
%     xmass=data(i).xmass(:,2);
%     plot(log10(xmass)+10,mfun(:,1).*(xmass).^alpha,'--',...
%         'color',colors(i), 'linewidth', 3)
% end

Nsnap=99; RunNum='6114'; name='msfun'; skip='';
datadir=['/mnt/A4700/data/',RunNum,'/Analysis/'];
% datadir=['/mnt/A4700/data/',RunNum,'/subcat',skip,'/anal/massfun/'];
% datadir=['/mnt/uv/',RunNum,'/subcat',skip,'/anal/massfun/'];
% datadir=['/mnt/charon/HBT/data/',RunNum,'/subcat',skip,'/anal/massfun/1Rvir/'];
[data,redshift]=read_massfun([datadir,'massfun_',num2str(Nsnap,'%03d'),'.',num2str(virtype,'%d')]);%,'.R1.0-2.0']);
nfun=numel(data);
Mhost=zeros(nfun,1);
for i=1:2:nfun
    Mhost(i)=data(i).Mhost;
end
for i=1:nfun
    mfun=data(i).mfunspecln/Mhost(i);
    xmass=data(i).xmass(:,2);
    plot(log10(xmass)+10,mfun(:,1).*(xmass).^alpha,'--',...
        'color',colors(i),'linewidth',3,...
        'displayname',[num2str(log10(data(i).Mbin(1))+10,'%2.1f'),'$\sim$',num2str(log10(data(i).Mbin(2))+10,'%2.1f')]);
    hold on;
end  

hold off;
set(gca,'yscale','log','yminortick','on');
xlabel('$\log(M_{sub}/(M_{\odot}/h)$)','interpreter','latex');
% ylabel('$dN/d\ln M_{sub}/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
ylabel('$M_{sub}dN/d\ln M_{sub}/M_{host}$','interpreter','latex');

% title([RunNum, ',z=',num2str(redshift,'%2.1f')]);
ylim([1e-3,5e-2]);
xlim([9.,14]);

% print('-depsc',[outputdir,'msfunln_6114vs6116S99V0.R1-2.eps']);
% fname=[name,'ln_',num2str(RunNum),'S',num2str(Nsnap),'V',num2str(virtype)];
% print('-depsc',[outputdir,fname,'.eps']);
% hgsave([outputdir,fname,'.fig']);