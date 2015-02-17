clear;
runnum='8213';
datadir=['/mnt/A4700/data/',runnum,'/subcat/anal/'];
data=importdata([datadir,'populations.dat'],'\t',1);
pop=data.data;
scaleF_file=['/mnt/A4700/data/',runnum,'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);    
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',10,...
    'DefaultLineMarkerSize',10,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
% semilogy(1./a-1,pop(:,3),'g-','displayname','subhalo');
loglog(a,pop(:,3),'g-','displayname','subhalo');
hold on;
loglog(a,pop(:,2),'r-','displayname','halo');
loglog(a,pop(:,4),'k-','displayname','birth');
loglog(a,pop(:,5),'b-','displayname','death');
loglog(a,pop(:,6),'c-','displayname','quasi-halos');
loglog(a,pop(:,7),'m-','displayname','splinters');
xlabel('Scale Factor');
ylabel('Number');
set(gca,'xlim',[0.05,1]);
legend('show','location','southeast');
outputdir='/home/kam/BT/data/show';
print('-depsc',[outputdir,'/population',runnum,'.eps']);