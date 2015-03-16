%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Projects/HBT/code/data/show/massfun/resim';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

M1=134.2184;R1= 179.3841;
virtype=1;

Nsnap=162;grpid=0;
RunName='AqA2';softeninghalo=4.8e-5;mp=1.000007e-06; R1=R1/1e3;
subdir=['/mnt/charon/HBT/data/',RunName,'/subcat/'];

% Nsnap=1023;grpid=0;
% RunName='AqA4';softeninghalo=0.25;mp=2.8679e-5;
% subdir=['/mnt/charon/HBT/data/',RunName,'/subcat/'];

scaleF_file=[subdir,'Redshift.dat'];
a=load_scaleF(scaleF_file);
z=1/a(Nsnap+1)-1;
%%
datadir=[subdir,'profile/logbin/'];
% datadir=['/mnt/A4700/data/',RunName,'/subcat/profile/logbin/'];
halo=readhalo_size(datadir,Nsnap,'halo');
haloprof=readhalo_prof_single(datadir,Nsnap,'halo',halo,grpid+1);
vs=diff(logspace(log10(max(softeninghalo,halo.rmax(grpid+1)*1e-2)),log10(halo.rmax(grpid+1)),halo.nbin(grpid+1)+1).^3)';
rho=(haloprof.ns+haloprof.no+haloprof.nb)./vs;
r=haloprof.rs;
erho=sqrt(rho)./sqrt(vs);
%%
mbdfile='/work/Projects/SubProf/data/A2Mbd.hdf5';
Mbd.m=h5read(mbdfile,'/mass');
Mbd.x=h5read(mbdfile,'/x')';
Mbd.mTVV=h5read(mbdfile,'/massTVV')';
Mbd.snapTVV=h5read(mbdfile, '/snapTVV')';
Mbd.flag=h5read(mbdfile,'/DirectInfall');
Mbd.r=sqrt(sum(Mbd.x.^2,2));
% Mbd.r=Mbd.r*1e3;
%%
% Mass Profile
nmin=100;
iInfall=2;%Mmax
% rbin=logspace(-1.5,0.7,20)*R1;
rbin=linspace(0.01, 3, 30)*R1;
rref=R1;
myfigure;
[~,i1]=min(abs(r-rref));
rhoDM=rho/rho(i1)*1;
h1=loglog(r/R1,rhoDM,'-','linewidth',5,'color',[0.8,0.8,0.8],'displayname','Halo'); hold on;
rhoSat=rhoDM.*(r/r(i1)).^1.3;
rhoSat(r>rref)=rhoDM(r>rref);
% loglog(r(r<rref)/R1, rhoSat(r<rref), 'g-', 'displayname', 'sat Model');

f=Mbd.m>nmin;
[xr1,n1,rhon1]=loghist(Mbd.r(f),rbin,[],[],1,3);
erhon1=sqrt(n1)./diff(rbin.^3);
[~,i1]=min(abs(xr1-rref));
% errorbar(xr1/R1,rhon1/rhon1(i1),erhon1/rhon1(i1),'-','linewidth',3,'color',[0.5,0.5,1], 'displayname','Sub');

f=f&Mbd.flag>0;
[xr11,n11,rhon11]=loghist(Mbd.r(f),rbin,[],[],1,3);
% plot(xr11/R1, rhon11/rhon11(i1),'-','linewidth',3,'color',[1,0.5,0.5],'displayname','Sub-Direct');

% type 1 vs. type 2
f=Mbd.mTVV(:,iInfall)>1000;
% f=f&Mbd.flag>0;
[xr3,n3,rhon3]=loghist(Mbd.r(f),rbin,[],[],1,3);
% plot(xr3/R1, rhon3/rhon3(i1),'k-','linewidth', 3, 'displayname','Ninfall>1e3, all');hold on;
f=Mbd.mTVV(:,iInfall)>1000&Mbd.m<1;
% f=f&Mbd.flag>0;
[xr3,n31,rhon31]=loghist(Mbd.r(f),rbin,[],[],1,3);
% plot(xr3/R1, rhon31/rhon3(i1),'-','color',[1,0.7,0.7],'displayname','Ninfall>1e3, orphan');hold on;
%%
% nmin=500;
% nmin=4000;
nmin=14000;
f=Mbd.mTVV(:,iInfall)>nmin&Mbd.m>0;
% f=f&Mbd.flag>0;
[xr3,n32,rhon32]=loghist(Mbd.r(f),rbin,[],[],1,3);
plot(xr3/R1, rhon32/rhon32(i1),'-','color',[0.7,1,0.7],'displayname','Ninfall>1e3, sub');


xlabel('D/Rvir');ylabel('$n/<n>$');
set(gca,'xscale','log','yscale','log');
legend('show','location','southwest')
% print('-depsc','/work/Projects/SubProf/plots/A4subprof_decompose.eps')
%%
myfigure();
nmin=2000;
nminZ0=20;
iInfall=1;
nbin=80;
f=Mbd.m>nminZ0&Mbd.mTVV(:,iInfall)>nmin;
% f=Mbd.flag>0&f;
loglog(Mbd.r(f)/R1, single(Mbd.m(f))./single(Mbd.mTVV(f,iInfall)),'.')
[xm,ym,yl,xmean,ymean,ysig]=skeleton(Mbd.r(f)/R1, single(Mbd.m(f))./single(Mbd.mTVV(f,iInfall)),logspace(-1,log10(3),nbin));
hold on;
plot(xm,ym,'k-',xm,yl,'g--');
xlabel('D/Rvir');ylabel('$m/m_0$');
xlim([1e-1,3]);
ff=~isnan(xm)&xm<1&xm>0.2;
p=polyfit(log(xm(ff)),log(ym(ff)),1);
h=plot(xm,exp(polyval(p,log(xm))),'m-');
% hh=plot(xm, 10^-0.3*xm.^1.2, 'm-');
plot(xlim(),nminZ0./nmin*[1,1],'r--');
legend(h,['slope=',num2str(p(1),'%.1f')],'location','southeast');
% print('-depsc','/work/Projects/SubProf/plots/A2strip.eps')
%%
[xm,ym,yl,xmean,ymean,ysig,count]=skeleton(log10(Mbd.r(f)/R1), log10(single(Mbd.m(f))./single(Mbd.mTVV(f,iInfall))),linspace(-1,log10(3),nbin));
ff=~isnan(xm);
cftool(xm(ff), ym(ff), [], (ysig(ff)./sqrt(count(ff))).^-1)
%%
f=Mbd.m>nminZ0&Mbd.mTVV(:,iInfall)>nmin&Mbd.r/R1>0.3&Mbd.r/R1<0.4;
data=log10(single(Mbd.m(f))./single(Mbd.mTVV(f,iInfall)));
hist(data, 20)
hold on;
plot([median(data),median(data)], ylim(), 'r-');
plot(log10(nminZ0/nmin)*[1,1], ylim(), 'r--');
%%
myfigure();
[xm,ym,yl,xmean,ymean,ysig]=skeleton(log10(Mbd.r(f)/R1), log10(single(Mbd.m(f))./single(Mbd.mTVV(f,iInfall))),linspace(-1,log10(3),nbin));
plot(xm,ym,'k-',xm,yl,'g--');
hold on;
errorbar(xmean,ymean,ysig,'ro');
ff=~isnan(xm);
p=polyfit((xmean(ff)),(ym(ff)),1);
h=plot(xmean,(polyval(p,xmean)),'b--');
legend(h,['slope=',num2str(p(1),'%.1f')],'location','southeast');
% print('-depsc','/work/Projects/SubProf/plots/A2strip_Err.eps')
%%
ff=~isnan(xm);
p=polyfit(xm(ff),log10(ym(ff)),1);
h=plot(xm, 10.^(polyval(p, xm)), 'r:');
legend(h,['slope=',num2str(p(1),'%.1f')],'location','southeast');
myfigure();
plot(xm,ysig,'-');
myfigure();
[xx,yy,n,s]=densitygrid(log10(Mbd.r(f)/R1),log10(single(Mbd.m(f))./single(Mbd.mTVV(f,iInfall))),[nbin,nbin],[-1,0.8],[-2,0.2]);
contourf(xx,yy,log10((n+1))); hold on;
x=xlim();
h=plot(x,2*x-0.4,'k-')
xlabel('$\log(D/R_{200})$');
ylabel('$\log(m/m_0)$');
legend(h, 'slope=2', 'location','southeast');
% print('-depsc','/work/Projects/SubProf/plots/A2strip_joint.eps')
%%
myfigure();
iInfall=1;
f=Mbd.mTVV(:,iInfall)>0&Mbd.r<R1;
% f=f&Mbd.flag>0;
[xm,ym,dyx]=linhist(log(single(Mbd.mTVV(f,iInfall))),20);
loglog(exp(xm), dyx, '-');
hold on;
dyx(isnan(xm))=[];
xm(isnan(xm))=[];
p=polyfit(xm, log(dyx), 1);
h=plot(exp(xm),exp(polyval(p,xm)),'r--');
legend(h,['slope=',num2str(p(1),'%.1f')],'location','southwest');
xlabel('$M_{\rm max}$(particles)');ylabel('${\rm d}N/{\rm d}\ln M$');
xlim([10, 1e6]);
print('-depsc','/work/Projects/SubProf/plots/A2MFinfall.eps')