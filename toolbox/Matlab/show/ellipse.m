t=0:0.1:2*pi+0.1;
a=2;b=1;
c=sqrt(a^2-b^2);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultTextFontSize',20,...
    'DefaultTextInterpreter','latex');
   
plot(a*sin(t),b*cos(t),'k');
hold on;
plot(-c,0,'o','markersize',8,'markerfacecolor','k','markeredgecolor','k');
text(-c,-0.2,'M');
plot(a*sin(pi/4),b*cos(pi/4),'o','markersize',4,'markerfacecolor','k','markeredgecolor','k');
plot(a*sin(pi/4+0.1),b*cos(pi/4+0.1),'o','markersize',4,'markerfacecolor','k','markeredgecolor','k');
text(a*sin(pi/4+0.1),b*cos(pi/4+0.1)+0.2,'m');
text(0.5,0.2,'$\overrightarrow{r}$');
% plot([-sqrt(3),2*sin(pi/4+0.1)],[0,cos(pi/4+0.1)]);
axis equal;
box off;
set(gca,'visible','off');

outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/ellipse.eps']);