%% data preparation
addpath(genpath('../post'));

Nsnap=59;RunNum='8213';
pmass=0.0620373e10;


datadir=['/mnt/A4700/data/',RunNum,'/subcat/profile/logbin'];
par=readhalo_param(datadir,Nsnap);
halo=readvirial_size(datadir,Nsnap,'halo');
c=par.c(:,1);
cerr=par.c(:,2);
M=halo.Mvir(:,1);

M=M.*pmass;
%%
x=logspace(11,16,5);
%A=8.4;x0=1.e12;beta=0.108; %Maccio08 WMAP5
% A=10.26;x0=1.e12;beta=0.114; %Maccio08 WMAP1
%A=7.96;x0=2e12;beta=0.091; %Duffy08
A=9;x0=1.3e13;beta=0.13; %Bullock
y=A*(x/x0).^-beta;
yu=y*exp(0.25);
yl=y/exp(0.25);
x=log10(x);
%%
f=(par.stat==0)&(cerr./c<0.5)&(par.Merr<0.2)&(M>1000*pmass);
%%
Mbin=logspace(log10(nanmin(M(f))),log10(nanmax(M(f))),10);
cm=zeros(9,1);mm=cm;
for i=1:9
    cm(i)=median(c(M>Mbin(i)&M<Mbin(i+1)&f));
    mm(i)=mean(M(M>Mbin(i)&M<Mbin(i+1)&f));
end
%%
plot(log10(M(f)),c(f),'.','displayname','halos');hold on;
plot(log10(mm),cm,'r-','displayname','median');
plot(x,y,'k-','displayname','$c=9(\frac{M}{M^*})^{\frac{\ }{\ }0.13}$');
hl=legend('show');set(hl,'interpreter','latex');
plot(x,yu,'k--');plot(x,yl,'k--');set(gca,'yscale','log');
xlabel('$log(M_{vir}/(M_{\odot}/h))$','interpreter','latex');ylabel('$c_{vir}$','interpreter','latex');
title('z=0');