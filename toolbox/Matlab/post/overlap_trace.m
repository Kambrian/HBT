runnum='6702DM';Nsnap=40;
% runnum='8213';Nsnap=59;
scaleF_file=['/mnt/A4700/data/',runnum,'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);   
% a=a(end:-1:1);

file=['/mnt/A4700/data/',runnum,'/subcat/anal/trace_overlap_',num2str(Nsnap,'%03d')];
tmp=load(file);
d=tmp(:,1:end/3);
% m1=tmp(:,end/3+1:end*2/3);
% m2=tmp(:,end*2/3+1:end);
n=size(d,1);
% m0=m2(1,:);

% dd=d(:,m0>1000);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',10,...
    'DefaultLineMarkerSize',10,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
plot(d);
set(gca,'ylim',[0,15]);
% set(gca,'xlim',[0,2]);
xlabel('Redshift');
ylabel('Seperation/Softening');
% 
% outputdir='/home/kam/BT/data/show';
% print('-depsc',[outputdir,'/overlap_orbit_',runnum]);
