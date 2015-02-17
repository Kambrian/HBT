%% data preparation
addpath(genpath('../post'));

Nsnap=59;
RunNum=8213;

file=['/mnt/A4700/data/',num2str(RunNum),'/subcat/anal/M_c_',num2str(Nsnap,'%03d')];
fid=fopen(file,'r');
Ngrp=fread(fid,1,'int32');
Mc=fread(fid,[2,Ngrp],'float32');
Ngrp2=fread(fid,1,'int32');
if Ngrp~=Ngrp2
    error('file corrupt');
end
M=Mc(1,:)'*1e10;
c=Mc(2,:)';
Mbin=logspace(log10(min(M)),log10(max(M)),10);
cm=zeros(9,1);mm=cm;
for i=1:9
    cm(i)=median(c(M>Mbin(i)&M<Mbin(i+1)));
    mm(i)=mean(M(M>Mbin(i)&M<Mbin(i+1)));
end

x=logspace(11,16,5);
x0=1.5e13;
y=9*(x/x0).^-0.13;
yu=y*exp(0.25);
yl=y/exp(0.25);
x=log10(x);
%%
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data';
figure;
set(gcf,...
    'DefaultLineLineWidth',2,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);

plot(log10(M(:)),c(:),'.','displayname','halos');hold on;
plot(log10(mm),cm,'r-','displayname','median');
plot(x,y,'k-','displayname','$c=9(\frac{M}{M^*})^{\frac{\ }{\ }0.13}$');
hl=legend('show');set(hl,'interpreter','latex');
plot(x,yu,'k--');plot(x,yl,'k--');set(gca,'yscale','log');
xlabel('$log(M_{vir}/(M_{\odot}/h))$','interpreter','latex');ylabel('$c_{vir}$','interpreter','latex');
title('z=0');

print('-depsc',[outputdir,'/M_c_',num2str(RunNum),'S',num2str(Nsnap),'.eps']);
hgsave([outputdir,'/M_c_',num2str(RunNum),'S',num2str(Nsnap),'.fig']);