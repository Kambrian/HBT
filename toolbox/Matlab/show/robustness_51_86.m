% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
outputdir='/home/kam/Projects/HBT/code/data/show';
runnum='6702DM';
scaleF_file=['/mnt/A4700/data/',runnum,'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);    

datadir=['/mnt/A4700/data/',runnum,'/subcat/anal/follow'];
BndMean=load([datadir,'/follow_051_86.mean']);
BndMin=load([datadir,'/follow_051_86.minpotH']);
BndW=load([datadir,'/follow_051_86.potW']);
%core
BndC10=load([datadir,'/follow_051_86.1']);
BndC5=load([datadir,'/follow_051_86.4']);
BndC3=load([datadir,'/follow_051_86.11']);
BndC2=load([datadir,'/follow_051_86.25']);
BndC14=load([datadir,'/follow_051_86.51']);
%adaptive; colums: [snapshot,sublen,srclen,srclen2,corefrac]
BndA10=load([datadir,'/follow_051_86.adpt10.0']);
BndA5=load([datadir,'/follow_051_86.adpt5.0']);
BndA3=load([datadir,'/follow_051_86.adpt3.0']);
BndA2=load([datadir,'/follow_051_86.adpt2.0']);
BndA14=load([datadir,'/follow_051_86.adpt1.4']);
%jump
BndJ2=load([datadir,'/follow_051_86.adpt2.0_skip2']);
BndJ5=load([datadir,'/follow_051_86.adpt2.0_skip5']);
BndJ10=load([datadir,'/follow_051_86.adpt2.0_skip10']);
%local; [snap,bndloc,srcloc,sublen,corefrac]
BndL=load([datadir,'/follow_051_86.adpt2.0.loc']);
%%
x=1./a(BndMean(:,1)+1)-1;
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',10,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultAxesYLim',[0,150000]);
h1=subplot(1,3,1);
yticks=0:50000:150000;
plot(x,BndMean(:,2),'k:','displayname','CoM Frame');hold on;
plot(x,BndMin(:,2),'m-.','displayname','MinPot Frame');
plot(x,BndW(:,2),'c--','displayname','CoP Frame');
plot(x,BndC10(:,2),'r-','displayname','CoreFrac=0.01');
% plot(x,BndC5(:,2),'r-','displayname','CoreFrac=0.04');
% plot(x,BndC3(:,2),'- *');
plot(x,BndC2(:,2),'b-','displayname','CoreFrac=0.25');
plot(x,BndC14(:,2),'g-','displayname','CoreFrac=0.5');
set(gca,'xdir','reverse','ytick',yticks);
ylabel('Bound Mass (particles)')

h2=subplot(1,3,2);
plot(x,BndMean(:,3),'k:');hold on;
plot(x,BndMin(:,3),'m-.');
plot(x,BndW(:,3),'c--');
plot(x,BndC10(:,3),'r-');
% plot(x,BndC5(:,3),'r-');
% plot(x,BndC3(:,3),'- *');
plot(x,BndC2(:,3),'b-');
plot(x,BndC14(:,3),'g-');
set(gca,'xdir','reverse','xtick',0:0.2:0.7,'ytick',yticks,'yticklabel','');
xlabel('Redshift');

h3=subplot(1,3,3);
plot(x,BndC2(:,2),'b:','displayname','Infalled,0.25');
hold on;
plot(x,BndC2(:,3),'b--','displayname','Progenitor,0.25');
plot(x,BndA10(:,2),'r-','displayname','Adaptive,0.01');
% plot(x,BndA5(:,2),'r-','displayname','Adaptive,0.04');
% plot(x,BndA3(:,2),'r- .');
plot(x,BndA2(:,2),'b-','displayname','Adaptive,0.25');
plot(x,BndA14(:,2),'g-','displayname','Adaptive,0.5');
set(gca,'xdir','reverse','xtick',0:0.2:0.7,'ytick',yticks,'yaxislocation','right');

% ylabel('Bound Mass (particles)')

set(h1,'position',[0.1,0.1,0.26,0.8]);
set(h2,'position',[0.36,0.1,0.26,0.8]);
set(h3,'position',[0.62,0.1,0.26,0.8]);
% h4=subplot(2,2,4);axis(h4,[0,0.8,0,10000]);
% set(h4,'position',[0.5,0.1,0.4,0.4],'box','on','xtick',0:0.2:0.7,'ytick',0:5000:10000,'xticklabel',[],'yticklabel',[],'xdir','reverse','yaxislocation','right');
legend(h1,'show','location','NorthEast');legend(h1,'boxoff');
legend(h3,'show','location','northeast');legend(h3,'boxoff');

set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,6,20,10]);
print('-depsc',[outputdir,'/robustness_51_86']);
% hgsave([outputdir,'/robustness_51_86.fig']);
%%
% figure;
% set(gcf,...
%     'DefaultLineLineWidth',2,'DefaultAxesLineWidth',.5,...
%     'DefaultAxesFontName','Helvetica',...
%     'DefaultAxesFontSize',20,...
%     'DefaultAxesTickLength',[0.02,0.02],... 
%     'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
myfigure;
ax(1)=subplot(2,1,1);
plot(1./a(BndA2(:,1)+1)-1,BndA2(:,2),'b- .','displayname','skip 0','linewidth',2,'markersize',8);
hold on;
plot(1./a(BndJ2(:,1)+1)-1,BndJ2(:,2),'rx','displayname','skip 1','markersize',8);
plot(1./a(BndJ5(:,1)+1)-1,BndJ5(:,2),'gd','displayname','skip 4','markersize',10);
plot(1./a(BndJ10(:,1)+1)-1,BndJ10(:,2),'ms','displayname','skip 9','markersize',12);
set(gca,'xdir','reverse','ylim',[0,10000]);
% xlabel('Redshift');
% set(gca,'ytick',2000:2000:10000);
% set(gca,'yticklabel',2000:2000:10000);
ylabel('Bound Mass (particles)');legend('show','location','NorthEast');%legend('boxoff');
ax(2)=subplot(2,1,2);
% plot(1./a(BndA2(:,1)+1)-1,BndA2(:,2),'b- .','displayname','skip 0','linewidth',2,'markersize',8);
% hold on;
plot(1./a(BndA2(:,1)+1)-1,zeros(size(BndA2,1),1),'b-'); hold on;
plot(1./a(BndJ2(:,1)+1)-1,BndJ2(:,2)./BndA2(1:2:end,2)-1,'rx','displayname','skip 1','markersize',8);hold on;
plot(1./a(BndJ5(:,1)+1)-1,BndJ5(:,2)./BndA2(1:5:end,2)-1,'gd','displayname','skip 4','markersize',10);
plot(1./a(BndJ10(:,1)+1)-1,BndJ10(:,2)./BndA2(1:10:end,2)-1,'ms','displayname','skip 9','markersize',12);
set(gca,'xdir','reverse');
xlabel('Redshift');
set(gca,'ylim',[-2e-3,2.5e-3])
set(gca,'yticklabel',[-0.002,0,0.002]);
ylabel('Relative Difference');
TightenSubplots(ax',[2,1]);
print('-depsc',[outputdir,'/trace_jump_51_86_res.eps']);
% hgsave([outputdir,'/trace_jump.fig']);

%%
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
plot(1./a(BndL(:,1))-1,BndL(:,4),'b- o','displayname','infalled halo remnant');
hold on;
% plot(1./a(BndL(:,1))-1,BndL(:,4),'k:','displayname','infalled subhalo');
plot(1./a(BndL(:,1))-1,BndL(:,2),'r--','displayname','local accretion');
set(gca,'xdir','reverse');
xlabel('Redshift');ylabel('Bound Mass (particles)');legend('show','location','NorthEast');%legend('boxoff');

print('-depsc',[outputdir,'/local_accrete_51_86.eps']);
% hgsave([outputdir,'/local_accrete.fig']);