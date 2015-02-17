%% data preparation
addpath(genpath('../post'));

Nsnap=99;grpid=0;
runs=[6402,6404,6409;6500,6501,6506;6600,6601,6602;6700,6701,6702];
c=zeros(size(runs));
cerr=c;
M=c;
for i=1:numel(runs)
datadir=['/mnt/A4700/data/',num2str(runs(i)),'/subcat/profile/logbin'];
par=readhalo_param(datadir,Nsnap);
halo=readvirial_size(datadir,Nsnap,'halo');
c(i)=par.c(grpid+1,1);
cerr(i)=par.c(grpid+1,2);
M(i)=halo.Mvir(grpid+1,1);
end

pmass=[5.45271e-6,5.29086e-06,5.23311e-06
       7.44348e-05,5.8378e-5,6.45321e-05
     0.000789798,0.000670876,0.000546019
     0.0130057,0.00834475,0.008848];

M=M.*pmass*1e10;
x=logspace(11,16,5);
x0=1.5e13;
y=9*(x/x0).^-0.13;
yu=y*exp(0.25);
yl=y/exp(0.25);
x=log10(x);

% c=[14.0,16.7,6.0
%     4.1,14.0,11.2
%     5.5,4.3,4.3
%     5.3,4.6,7.3];
% M=[100.7274,99.7226,94.8864;
%     994.0477,1011.4532,1013.251;
%     9483.4238,11126.1074,10329.4365;
%     181948.0156,123550.6562,111621.6641;];
% M=M*1e10;
%%
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data';
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

errorbar(log10(M(:)),c(:),cerr(:),'o','displayname','halos');hold on;plot(x,y,'displayname','$c=9(\frac{M}{M^*})^{\frac{\ }{\ }0.13}$');
hl=legend('show');set(hl,'interpreter','latex');
plot(x,yu,'--');plot(x,yl,'--');set(gca,'yscale','log');
xlabel('$log(M_{vir}/(M_{\odot}/h))$','interpreter','latex');ylabel('$c_{vir}$','interpreter','latex');
title('z=0');

print('-depsc',[outputdir,'/M_c_resim_z0.eps']);
hgsave([outputdir,'/M_c_resim_z0.fig']);