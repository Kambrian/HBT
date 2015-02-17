%% data preparation
% to add: 6702low, cosmo
outputdir='/home/kam/Projects/HBT/code/data/show/massfun/resim';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

markers=['s--';'o--';'d--';'^--';'<--';'>--';'p--';'x--';'+--';'.--'];
colors=['r';'g';'b';'c';'m';'k';'r';'g';'b';'c'];

virtype=0;
rmin=0;rmax=1;
nbin=8;
Nsnap=[37,39:20:99];grpid=0;
runs={'6702DM','6702DMNSM'};
RunName='6702DMNSMrat';
scaleF_file=['/mnt/A4700/data/',runs{1},'/subcat/Redshift.dat'];
a=load_scaleF(scaleF_file);
z=1./a(Nsnap+1)-1;
 
xmin=0.2;xmax=1.2e4;type='';
% xmin=1e-5;xmax=1;type='norm';
xref=logspace(10,13,5);
N=numel(Nsnap);
M0=zeros(1,N);R0=M0;
xmass=cell(1,N);mfunspec=cell(1,N);mfunspecln=cell(1,N);mfuncum=cell(1,N);
datadir=['/mnt/A4700/data/',runs{1},'/subcat/anal/massfun/'];
for i=1:N
    [Mlist,M0(i),R0(i)]=Msublist_in_radii(datadir,virtype,Nsnap(i),grpid,rmin,rmax,'bindata',type);
    [xmass{i},mfunspec{i},mfunspecln{i},mfuncum{i}]=mass_sumcount(Mlist,nbin,xmin,xmax);
end

M0=zeros(1,N);R0=M0;
xmass2=cell(1,N);mfunspec2=cell(1,N);mfunspecln2=cell(1,N);mfuncum2=cell(1,N);
datadir=['/mnt/A4700/data/',runs{2},'/subcat/anal/massfun/'];
for i=1:N
    [Mlist,M0(i),R0(i)]=Msublist_in_radii(datadir,virtype,Nsnap(i),grpid,rmin,rmax,'bindata',type);
    [xmass2{i},mfunspec2{i},mfunspecln2{i},mfuncum2{i}]=mass_sumcount(Mlist,nbin,xmin,xmax);
end
%% plot specific mass function
figure;%('visible','off');
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

ax1=subplot(2,1,1);
for i=2:N
    mfun=mfunspec{i}./mfunspec2{i};
    plot(log10(xmass{i}(:,2))+10,mfun(:,1),markers(i,:),...
        'color',colors(i),...
        'displayname',['z=',num2str(z(i),'%2.1f')]);
    hold on;
end

% xlabel('$log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$dM/dM_{NSM}$','interpreter','latex');
hl=legend('show','location','northwest');set(hl,'interpreter','latex');

% print('-depsc',[outputdir,'/msfun_',RunName,'.eps']);
%% plot cumulative mass function
% figure;%('visible','off');
% set(gcf,...
%     'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
%     'DefaultAxesFontName','Helvetica',...
%     'DefaultAxesFontSize',20,...
%     'DefaultAxesTickLength',[0.02,0.02],... 
%     'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
% set(gcf,'DefaultLineMarkerSize',6);

ax2=subplot(2,1,2);
for i=2:N
    mfun=mfuncum{i}./mfuncum2{i};
    plot(log10(xmass{i}(:,1))+10,mfun(:,1),markers(i,:),...
        'color',colors(i),...
        'displayname',['z=',num2str(z(i),'%2.1f')]);
    hold on;
end
    
xlabel('$\log(M_{sub}/(M_{\odot}/h))$','interpreter','latex');
ylabel('$M(>M_{sub})/M_{NSM}(>M_{sub})$','interpreter','latex');
% hl=legend('show','location','northwest');set(hl,'interpreter','latex');
TightenSubplots([ax1;ax2]);
set([ax1,ax2],'xlim',[9,14.5]);

% print('-depsc',[outputdir,'/msfuncum_',RunName,'.eps']);
print('-depsc',[outputdir,'/msfuncombine_',RunName,'.eps']);