%% load data
clear;
addpath(genpath('../post'));

markers=['o:';'s:';'d:';'^:';'<:';'x:';'+:';'o:';'s:';'d:';'^:';'<:';'x:';'+:'];
% markers=['.-';'.-';'.-';'.-';'.-';'.-'];
colors=['r';'g';'b';'c';'m';'k';'r';'g';'b';'c';'m';'k'];

snaps=59:-6:35;
% dn=[1,2,4,6,8,10,12,16];
dn=[1,2,4,6,8,10,16];
amp=zeros(numel(snaps),numel(dn)); %relative amplitude, averaged over hostmass bins
eamp=amp; % error of amp
dt=amp;   % dlna
%%
time=1;
Nsnap=snaps(time); %Nsnap=59;

RunNum=8213;name='msfun';
virtype=0;
outputdir=['/home/kam/Projects/HBT/code/data/show/massfun/',num2str(RunNum),'/'];
% outputdir=['/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/',num2str(RunNum),'/'];

scaleF_file=['/mnt/A4700/data/',num2str(RunNum),'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);z=1./a-1; 
[w,dlt,d]=collapse_barrier(a,0.3);

dlna=log(a(Nsnap+1)./a(Nsnap+1-dn));
dt(time,:)=dlna';
x=dlna/(1.3*a(Nsnap+1)+0.508);
dz=z(Nsnap+1-dn)-z(Nsnap+1);
dw=w(Nsnap+1-dn)-w(Nsnap+1);
dlnw=log(w(Nsnap+1-dn)/w(Nsnap+1));
    
data=cell(numel(dn),1);
datadir=['/mnt/A4700/data/',num2str(RunNum),'/subcat/anal/massfun/'];
[data{1},redshift]=read_massfun([datadir,'massfun_',num2str(Nsnap,'%03d'),'.',num2str(virtype,'%d')]);
nfun=numel(data{1});
if time==5, nfun=2; end
Mhost=zeros(nfun,1);
for i=1:nfun
    Mhost(i)=data{1}(i).Mhost;
end

for i=2:numel(dn)
datadir=['/mnt/A4700/data/',num2str(RunNum),'/subcat',num2str(dn(i)),'/anal/massfun/'];
[data{i},redshift]=read_massfun([datadir,'massfun_',num2str(Nsnap,'%03d'),'.',num2str(virtype,'%d')]);
end
%%
ifun=1;
for i=1:numel(dn)
    mfun=data{i}(ifun).mfunspec/Mhost(ifun);
    xmass=data{i}(ifun).xmass(:,2);
    plot(log10(xmass),mfun(:,1),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',[num2str(log10(data{i}(ifun).Mbin(1))+10,'%2.1f'),'$\sim$',num2str(log10(data{i}(ifun).Mbin(2))+10,'%2.1f')]);
    hold on;
%     cftool(log10(xmass), log10(mfun(:,1)))
end
hold off;
set(gca,'yscale','log','xminortick','on','yminortick','on');
xlabel('$\log(M_{sub}/M_{vir})$','interpreter','latex');
ylabel('$dN/d(M_{sub}/M_{vir})/N_{host}$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(redshift,'%2.1f')]);
%% plot logspaced specific mass function
y=zeros(nfun,numel(dn));ey=y;
flagplot=1;
mfun=cell(nfun,numel(dn));
for i=1:nfun
figure;
for j=1:numel(dn)
[y(i,j),ey(i,j),mfun{i,j},err]=fit_massfun_ratio(data{j}(i),data{1}(i),Mhost(i),flagplot);
end
end
y(:,1)=1;ey(:,1)=0.001;   

ym=sum(y./ey.^2)./sum(ey.^-2);
eym=sqrt(1./sum(ey.^-2));
ym(1)=1;eym(1)=min(eym(2:end))/10;
amp(time,:)=ym;
eamp(time,:)=eym;
%%
cftool(dlna,ym,1./eym.^2);
% cftool(x,ym,1./eym.^2);
%%
figure;errorbar(dn,ym,eym);
figure;errorbar(dn,y(1,:),ey(1,:),'r');
hold on;
errorbar(dn,y(2,:),ey(2,:),'g');
errorbar(dn,y(3,:),ey(3,:),'b');
%%
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);
set(gcf,'DefaultTextInterpreter','latex');
xmass=data{1}(1).xmass(:,2);
for i=1:numel(dn)
 plot(log10(xmass)+10,mfun{2,i}(:,1),markers(i,:),...
        'color',colors(i),'markerfacecolor',colors(i),...
        'displayname',num2str(dn(i)));
    hold on;
 plot(log10([xmass(1),xmass(end)])+10,[y(2,i),y(2,i)],'-','color',colors(i));
end
xlabel('$\log(M_{sub}/(M_{\odot}/h)$)','interpreter','latex');
% ylabel('$dN/d\ln M_{sub}/M_{host}\times(10^{10}M_{\odot}/h)$','interpreter','latex');
ylabel('$Relative{\ } Amplitude$','interpreter','latex');
hl=legend('show','location','southwest');set(hl,'interpreter','latex');
% title(['z=',num2str(redshift,'%2.1f')]);

% print('-depsc',[outputdir,fname,'.eps']);
%%
%fun=@(x) A*exp(-(x/b).^2);
A =[1.001,0.001;
    1.001,0.001;
    1.002,0.001;
    1.002,0.001;
    1.003,0.001];%Amplitude and error of A for guassian fit
b=[1.80,0.02;
    1.48,0.03;
    1.26,0.03;
    1.06,0.02;
    0.92,0.03;] %width and error

markers=['o';'s';'d';'^';'<';'x';'+';'o';'s';'d';'^';'<';'x';'+'];
colors=['r';'g';'b';'c';'m';'k';'r';'g';'b';'c';'m';'k'];
myfigure;

for time=1:5
errorbar(dt(time,:),amp(time,:),eamp(time,:),markers(time,:),...
        'color',colors(time),'markerfacecolor',colors(time),...
        'displayname',['z=',num2str(1/a(snaps(time)+1)-1,'%2.1f')]);
hold on;
end
legend('show','location','southwest');
for time=1:5
xplot=linspace(0,dt(time,end)*1.2,20);
plot(xplot,A(time,1)*exp(-(xplot/b(time,1)).^2),'color',colors(time),'displayname','fit');
end

xlabel('$\Delta\ln a$','interpreter','latex');
ylabel('Relative Amplitude','interpreter','latex');

print('-depsc',[outputdir,'msfun_time_resolution.eps']);
%%
t=a(snaps+1);
cftool(t,b(:,1),b(:,2).^-2);
%%
k=[1.29,0.02]; %slope and error
b0=[0.51,0.01]; 
%b=k*t+b0;
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);
set(gcf,'DefaultTextInterpreter','latex');
errorbar(t,b(:,1),b(:,2),'o');
hold on;
plot(t,k(1)*t+b0(1),'-');
xlabel('$a$');
ylabel('$b$');
print('-depsc',[outputdir,'msfun_time_resolution_width.eps']);
%%
myfigure;

errorbar(xa0,y0(1,:),ey0(1,:),'ro','displayname','12~13');
hold on;
errorbar(xa0,y0(2,:),ey0(2,:),'bs','displayname','13~14');
errorbar(xa0,y0(2,:),ey0(2,:),'g<','displayname','14~15');
A=1.001;x0=1.8;
fun=@(x) A*exp(-(x/x0).^2);
xplot=linspace(0,xa0(end)*1.2,20);
plot(xplot,fun(xplot),'k','linewidth',2,'displayname','z=0,x0=1.8');

xlabel('$\Delta\ln a$','interpreter','latex');
ylabel('Relative Amplitude');
%%
% figure;
errorbar(xa,y(1,:),ey(1,:),'ro','markerfacecolor','r','displayname','12~13');
hold on;
errorbar(xa,y(2,:),ey(2,:),'gs','markerfacecolor','g','displayname','13~14');
A=1.003;x0=0.9;
fun=@(x) A*exp(-(x/x0).^2);
xplot=linspace(0,xa(end)*1.2,20);
plot(xplot,fun(xplot),'k','displayname','z=2.1,x0=0.9');
legend('show');
% errorbar(x,ym,eym);

print('-depsc',[outputdir,'msfun_time_resolution.eps']);