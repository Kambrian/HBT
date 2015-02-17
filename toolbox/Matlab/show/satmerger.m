% runnum='6702DMNSM';Nsnap=45;
% file=['/mnt/A4700/data/',runnum,'/subcat/anal/MatchTo/HBT/',num2str(Nsnap),'_',num2str(Nsnap)];
% NSM2BT=load(file);
% nsm=sortrows(NSM2BT,2);
% nsm(nsm(:,3)==0,:)=[];
% nsm(1:100,:)

partmass=0.0104089;
runnum='6702DM';
scaleF_file=['/mnt/A4700/data/',runnum,'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);  

file=['/mnt/A4700/data/',runnum,'/subcat/anal/history_043_5562.txt'];
data=importdata(file,',',1);
% 
% runnum='6702DMNSM';
% file=['/mnt/A4700/data/',runnum,'/subcat/anal/history_045_3.txt'];
% data2=importdata(file,',',1);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');

axes('position',[0.15,0.5,0.8,0.45]);
plot(a(1+data.data(:,1)),data.data(:,2)/partmass,'-');
hold on;
% plot(data2.data(:,1),data2.data(:,2)/partmass,'-');
grid on;
% set(gca,'yscale','log');
set(gca,'xaxis','top');
set(gca,'xlim',[a(1),a(100)]);
ylabel('Mass/particles');

axes('position',[0.15,0.15,0.8,0.35]);
semilogy(a(1+data.data(:,1)),data.data(:,6));
set(gca,'ylim',[0.01,10],'xlim',[a(1),a(100)]);
grid on;
xlabel('scale factor');
ylabel('$D/R_{vir}$');

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/subhist_043_5562.eps']);