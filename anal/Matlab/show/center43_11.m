%CoM MinPot CoP 0.01 0.04 0.1 0.25 0.5 
%for 6702DM snap43 grp11

dat=[
150339,168943,148965,-50.8955,767.657,153.141
150370,168933,148959,-50.8955,767.657,153.141
150357,168939,148961,-61.4487,768.342,166.424
150370,168933,148960,-66.672,753.562,191.246
150370,168933,148960,-72.6299,737.301,188.485
150370,168934,148960,-71.0995,737.365,194.468
150371,168933,148962,-75.9327,775.858,186.975
150371,168930,148960,-70.8517,769.177,180.188
];
Cen=dat(:,1:3);
VCen=dat(:,4:6);

CoM=[150370,168935,148959];
VCoM=[-72.4781,740.146,196.474];
Rvir=617.783;
Vvir=475.051;

d=Cen-repmat(CoM,size(Cen,1),1);
d=sqrt(sum(d.^2,2))/Rvir;
v=VCen-repmat(VCoM,size(VCen,1),1);
v=sqrt(sum(v.^2,2))/Vvir;

myfigure;
plot(d(1),v(1),'k^','markerfacecolor','k','markersize',10,'displayname','CoM Frame');hold on;
plot(d(2),v(2),'ms','markerfacecolor','m','markersize',10,'displayname','MinPot Frame');
plot(d(3),v(3),'cd','markerfacecolor','c','markersize',10,'displayname','CoP Frame');
plot(d(4),v(4),'ro','markerfacecolor','r','markersize',6,'displayname','CoreFrac=0.01');
% plot(x,BndC5(:,2),'r-','displayname','CoreFrac=0.04');
% plot(x,BndC3(:,2),'- *');
plot(d(7),v(7),'bo','markerfacecolor','b','markersize',10,'displayname','CoreFrac=0.25');
plot(d(8),v(8),'go','markerfacecolor','g','markersize',14,'displayname','CoreFrac=0.5');
xlabel('$\Delta R/R_{vir}$','interpreter','latex');
ylabel('$\Delta V/V_{vir}$','interpreter','latex');
legend('show','location','southeast');

% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
outputdir='/home/kam/Projects/HBT/code/data/show';
print('-depsc',[outputdir,'/center_43_11.eps']);