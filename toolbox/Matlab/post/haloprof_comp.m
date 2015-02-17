% clear;
% a=load('outputs_zoom.txt');
datadir='/mnt/A4700/data/6700/subcat/profile/logbin';
Nsnap=59;
halo=readhalo_size(datadir,Nsnap,'halo');
ghalo=readhalo_size(datadir,Nsnap,'ghalo');
haloprof=readhalo_prof(datadir,Nsnap,'halo',halo);
ghaloprof=readhalo_prof(datadir,Nsnap,'ghalo',ghalo);
dynamprof=readdynam_prof(datadir,Nsnap,halo);
dynamprof_sph=readdynam_prof_sph(datadir,Nsnap,halo);
%%
figure;subplot(4,2,1);hist(log10(ghalo.mass(logical(ghalo.mass))./halo.mass(logical(ghalo.mass))));xlabel('log10(GasNum/DMNum),raw');
title(['snap',num2str(Nsnap),', z=',num2str(1/a(Nsnap+1)-1,'%2.1f'),', a=',num2str(a(Nsnap+1),'%3.2g')]);
Fbin=0.1:0.2:1.4;
haloflag=repmat([1,0,0],halo.ngrps,1);
% haloflag=haloflag&logical(ghalo.Mvir_this);
haloflag=haloflag&logical(halo.Mvir);
haloflag=haloflag&(repmat(~halo.flag_fakehalo,1,3));
haloflag=haloflag&(~halo.flag_badvir);
haloflag=haloflag&(repmat(halo.mass>100,1,3));
subplot(4,2,2);hist((ghalo.mass(haloflag(:,1))./halo.mass(haloflag(:,1))),Fbin);xlabel('(GasNum/DMNum),filtered');xlim([Fbin(1)-0.1,Fbin(end)+0.1]);
subplot(4,2,3);hist((ghalo.Mvir_all(haloflag)./halo.Mvir(haloflag)),Fbin);xlabel('(GasNum_{vir}/DMNum_{vir}),filtered');xlim([Fbin(1)-0.1,Fbin(end)+0.1]);
subplot(4,2,4);hist((ghalo.Mvir_this(haloflag)./ghalo.Mvir_all(haloflag)));xlabel('(GasNum_{halo,vir}/GasNum_{vir}),filtered');
subplot(4,2,[5,7]);loglog(halo.Mvir(haloflag),ghalo.Mvir_all(haloflag)./halo.Mvir(haloflag),'.');xlabel('halo vir mass');ylabel('vir gas fraction');
ecc=sqrt(sum((halo.CoM-halo.Cen).^2,2))./halo.Rvir(:,1);
subplot(4,2,[6,8]);semilogy(ecc(haloflag(:,1)),ghalo.Mvir_all(haloflag)./halo.Mvir(haloflag),'.');
xlabel('eccentricity=R(haloCoM~mainsubCen)/haloRvir');ylabel('vir gas fraction');

% set(gcf,'paperunits','normalized','paperposition',[0,0,1,1]);
% print('-dpdf',['gasfraction_6702_',num2str(Nsnap,'%03d')]);
% hgsave(['gasfraction_6702_',num2str(Nsnap,'%03d')]);
%%
c=find((ghalo.Mvir_all./halo.Mvir<0.5)&haloflag);
d=ghalo.Mvir_all(c,1)./halo.Mvir(c,1);
d(1:10)
c(1:10)

%%

pmass=0.0088475;
grpind=4;
rvir=halo.Rvir(grpind,1);

ns=haloprof.ns{grpind};
no=haloprof.no{grpind};
nb=haloprof.nb{grpind};
rs=haloprof.rs{grpind};
ro=haloprof.ro{grpind};
rb=haloprof.rb{grpind};

vs=diff(logspace(log10(max(1.5,halo.rmax(grpind)*1e-2)),log10(halo.rmax(grpind)),halo.nbin(grpind)+1).^3)';
vs=vs*(4*pi/3);
svs=cumsum(vs);
na=ns+no+nb;
sna=cumsum(na);
sns=cumsum(ns);
figure; 
loglog(rs/rvir,pmass*ns./vs,'r- .');
hold on;
loglog(ro/rvir,pmass*no./vs,'g: .');
loglog(rb/rvir,pmass*nb./vs,'b-- .');
loglog(rs/rvir,pmass*na./vs,'k-. .');

% % req=2.5;c=4.65;beta=1.15;%for grpind=1
% % req=4;c=4.2;beta=1.2;%for grpind=2
% req=2;c=6.25;beta=1.1;%for grpind=3;
% f=@(x) gas_dens(x,req,beta);
% x=rs/rvir*c;
% % logspace(-1,1,20);
% y=zeros(size(x));
% for i=1:numel(x)
% y(i)=f(x(i));
% end
% hold on;
% loglog(x/c,y.*(ns)*pmass./vs,'-- o');

switch grpind
    case 1
req=1.49;Rhoeq=2.59e-5;c=4.65;beta=0.9;%for grpind=1
    case 2
req=0.295*4.2;Rhoeq=2.094e-5;c=4.2;beta=1.2;%for grpind=2
    case 3
req=0.194*6.25;Rhoeq=6.1e-5;c=6.25;beta=1;%for grpind=3;
end
f=@(x) gas_dens_cosig(x,req,beta,grpind);
x=rs/rvir*c;
% logspace(-1,1,20);
y=zeros(size(x));
for i=1:numel(x)
y(i)=f(x(i));
end
hold on;
loglog(x/c,y*Rhoeq,'r-- o');

ns=ghaloprof.ns{grpind};
no=ghaloprof.no{grpind};
nb=ghaloprof.nb{grpind};
rs=ghaloprof.rs{grpind};
ro=ghaloprof.ro{grpind};
rb=ghaloprof.rb{grpind};

vs=diff(logspace(log10(max(1.5,ghalo.rmax(grpind)*1e-2)),log10(ghalo.rmax(grpind)),ghalo.nbin(grpind)+1).^3)';
vs=vs*(4*pi/3);
svs=cumsum(vs);
na=ns+no+nb;
sna=cumsum(na);
sns=cumsum(ns);
hold on;
loglog(rs/rvir,pmass*ns./vs,'r- *');
loglog(ro/rvir,pmass*no./vs,'g: *');
loglog(rb/rvir,pmass*nb./vs,'b-- *');
loglog(rs/rvir,pmass*na./vs,'k-. *');
hold off;
xlabel('r/rvir');
ylabel('density');
% legend('dm halo','gas model,{\beta}_{eff}=1.12','gas halo','gas all');
legend('halo','other','back','all','gas model,{\beta}_{eff}=1.12','gas model,{\beta}_{eff}=0.88','ghalo','gother','gback','gall');
%%
%isothermal equilibrated density prof for gas
%Rho0:Rho_gas(r=0)
%c: concentration c=rvir/rs
%x: r/rs  (rs is NFW scale radius)
%beta:Tgas=beta*Tvir
f=@(Rho0,x,c,beta) Rho0*exp(3*c/beta/(log(1+c)-c/(1+c))*(log(1+x)./x-1));
Rhos=0.025;c=4.65;xeq=2.5;beta=1;
% Rhos=0.019;c=4.2;xeq=2;beta=1;
% Rhos=0.048;c=6.25;xeq=2;beta=1;
x=rs/rvir*c;
Rho_gas0=Rhos/xeq/(1+xeq)^2/f(1,xeq,c,beta);
hold on; loglog(x/c,f(Rho_gas0,x,c,beta)*pmass,'g-- o');
NFW=@(Rhos,x) Rhos./x./(1+x).^2;
hold on; loglog(x/c,NFW(Rhos,x)*pmass,'r-- o');
% xeq=1;beta=1.7;
xeq=4;beta=1.5;
Rho_gas0=Rhos/xeq/(1+xeq)^2/f(1,xeq,c,beta);
hold on; loglog(x/c,f(Rho_gas0,x,c,beta)*pmass,'b-- o');
%%
legend('halo','other','back','all','ghalo','gother','gback','gall','Tgas=Tvir','Tgas=1.8Tvir');
%%

figure;loglog(rs/rvir,(ghaloprof.ns{grpind}+ghaloprof.no{grpind}+ghaloprof.nb{grpind})./(haloprof.ns{grpind}),'- .')
req=2.5;c=4.65;beta=1.1;%for grpind=1
% req=4;c=4.2;beta=1.1;%for grpind=2
% req=2;c=6.25;beta=1.1;
f=@(x) gas_dens(x,req,beta);
x=rs/rvir*c;
% logspace(-1,1,20);
y=zeros(size(x));
for i=1:numel(x)
y(i)=f(x(i));
end
hold on;
loglog(x/c,y,'-- .');
hold off;
xlabel('R/Rvir');ylabel('{\rho}_{gasall}/{\rho}_{dmhalo}');

%%
grpind=4;
rvir=halo.Rvir(grpind,1);
vm2s=dynamprof.vm2s{grpind};
vm2o=dynamprof.vm2o{grpind};
vm2b=dynamprof.vm2b{grpind};
vsms=dynamprof.vsms{grpind};
vsmo=dynamprof.vsmo{grpind};
vsmb=dynamprof.vsmb{grpind};

ns=haloprof.ns{grpind};
no=haloprof.no{grpind};
nb=haloprof.nb{grpind};
rs=haloprof.rs{grpind};
ro=haloprof.ro{grpind};
rb=haloprof.rb{grpind};
vsma=(vsms.*ns+vsmo.*no+vsmb.*nb)./(ns+no+nb);

vs=diff(logspace(log10(max(1.5,halo.rmax(grpind)*1e-2)),log10(halo.rmax(grpind)),halo.nbin(grpind)+1).^3)';
vs=vs*(4*pi/3);

figure; subplot(2,2,1);
loglog(rs/rvir,vsms-vm2s,'r- .');
hold on;
loglog(ro/rvir,vsmo-vm2o,'g: .');
loglog(rb/rvir,vsmb-vm2b,'b-- .');
loglog(rs/rvir,vsma-vm2s,'k-. .');

s=@(x,sig2_rs,alpha) 4^(2/3)*sig2_rs./(x.^(1-alpha).*(1+x).^2).^(2/3);
c=4.65;sig2_rs=5.2e6;alpha=1.77;
loglog(rs/rvir,s(rs/rvir*c,sig2_rs,alpha),'m- o');

gvm2s=dynamprof.gvm2s{grpind};
gvm2o=dynamprof.gvm2o{grpind};
gvm2b=dynamprof.gvm2b{grpind};
gvsms=dynamprof.gvsms{grpind};
gvsmo=dynamprof.gvsmo{grpind};
gvsmb=dynamprof.gvsmb{grpind};
gvsma=(gvsms.*ns+gvsmo.*no+gvsmb.*nb)./(ns+no+nb);

loglog(rs/rvir,gvsms-gvm2s,'r- s');
hold on;
loglog(ro/rvir,gvsmo-gvm2o,'g: s');
loglog(rb/rvir,gvsmb-gvm2b,'b-- s');
loglog(rs/rvir,gvsma-gvm2s,'k-. s');

Ums=dynamprof.Ums{grpind};
Umo=dynamprof.Umo{grpind};
Umb=dynamprof.Umb{grpind};

Uma=(Ums.*ns+Umo.*no+Umb.*nb)./(ns+no+nb);
loglog(rs/rvir,2*Ums,'r- *');
loglog(ro/rvir,2*Umo,'g: *');
loglog(rb/rvir,2*Umb,'b-- *');
loglog(rs/rvir,2*Uma,'k-. *');
hold off;
legend('halo','other','back','all','model','gvhalo','gvother','gvback','gvall','ghalo','gother','gback','gall','location','southwest');
xlabel('R/Rvir');
ylabel('{\sigma}_{dm}^2 or 2*Ugas');


subplot(2,2,2);
loglog(rs/rvir,haloprof.ns{grpind}./vs./(Ums).^(3/2),'r- *');hold on;loglog(rs/rvir,haloprof.ns{grpind}./vs./(vsms-vm2s).^(3/2),'- .');
loglog(rs/rvir,(haloprof.ns{grpind}+haloprof.nb{grpind}+haloprof.no{grpind})./vs./(dynamprof_sph.vd_rtf{grpind}(1,:)').^(3/2),'g- o');
rflag=logical(rs<rvir);
x=log(rs(rflag)/rvir);
y=log(haloprof.ns{grpind}./vs./(Ums).^(3/2));y=y(rflag);
p=polyfit(x,y,1);hold on;loglog(rs/rvir,exp(polyval(p,log(rs/rvir))),'r--');agas=p(1);
y=log(haloprof.ns{grpind}./vs./(vsms-vm2s).^(3/2));y=y(rflag);
p=polyfit(x,y,1);hold on;loglog(rs/rvir,exp(polyval(p,log(rs/rvir))),'--');adm=p(1);
y=log((haloprof.ns{grpind}+haloprof.nb{grpind}+haloprof.no{grpind})./vs./(dynamprof_sph.vd_rtf{grpind}(1,:)').^(3/2));y=y(rflag);
p=polyfit(x,y,1);hold on;loglog(rs/rvir,exp(polyval(p,log(rs/rvir))),'g--');admr=p(1);
legend('gas','dm','dmr',['{\alpha}_{gas}=',num2str(agas)],['{\alpha}_{dm}=',num2str(adm)],['{\alpha}_{dmr}=',num2str(admr)],'location','southwest');
xlabel('R/Rvir');
ylabel('{\rho}_{dm}/{\sigma}_{dm}^3 or {\rho}_{dm}/(Ugas)^{3/2}');

subplot(2,2,[3,4]);
semilogx(rs/rvir,(2*Ums)./(vsms-vm2s),'- .');
hold on;
% semilogx(rs/rvir,(2*Ums+(gvsms-gvm2s))./(vsms-vm2s),'- s');
semilogx(rs/rvir,(2*Ums+dynamprof_sph.gvd_rtf{grpind}(1,:)'*3)./(vsms-vm2s),'- s');
semilogx(rs/rvir,dynamprof_sph.vd_rtf{grpind}(1,:)'*3./(vsms-vm2s),': o');
rflag=logical(rs<rvir);
x=log(rs(rflag)/rvir);
y=2*Ums./(vsms-vm2s);y=y(rflag);
rat=polyfit(x,y,0);
hold on;loglog([0.05,5],[rat,rat],'r--');
y=(2*Ums+(gvsms-gvm2s))./(vsms-vm2s);y=y(rflag);
rat2=polyfit(x,y,0);
hold on;loglog([0.05,5],[rat2,rat2],'r--');
y=dynamprof_sph.vd_rtf{grpind}(1,:)'*3./(vsms-vm2s);y=y(rflag);
rat3=polyfit(x,y,0);
hold on;loglog([0.05,5],[rat3,rat3],'r--');
xlabel('R/Rvir');
ylabel('${ 2U_{gas} \over {\sigma}_{DM}^2}={{\sigma}_{gas\_ion}^2  \over {\sigma}_{DM}^2}$','interpreter','latex','fontsize',15);
legend('Tgas/Tdyn','Teff/Tdyn','Tdyn_r/Tdyn',['{\beta}=',num2str(rat)],['{\beta}_{eff}=',num2str(rat2)],['{\beta}_{r}=',num2str(rat3)]);

%%
grpind=2;
pmass=0.008848;
switch grpind
    case 1
Rhos=0.025*pmass;Rs=331;c=4.65;sig2_rs=1.7e6;alpha=1.92;
    case 2
Rhos=0.019*pmass;Rs=183;c=4.2;sig2_rs=3e5;alpha=1.7;
    case 3
Rhos=0.048*pmass;Rs=103;c=6.25;sig2_rs=2.5e5;alpha=1.88;
end
s=@(x,sig2_rs,alpha) 4^(2/3)*sig2_rs./(x.^(1-alpha).*(1+x).^2).^(2/3);
g=@(x,Rhos,Rs) 43007.1*4*pi*Rhos*Rs^3*(log(1+x)-x./(1+x))./x.^2/Rs^2;
rvir=halo.Rvir(grpind,1);
rs=haloprof.rs{grpind};
figure;
subplot(2,2,1);semilogx(rs/rvir,dynamprof_sph.vm_xyz{grpind}','- .');
subplot(2,2,2);semilogx(rs/rvir,dynamprof_sph.vm_rtf{grpind}','- .');
subplot(2,2,3);loglog(rs/rvir,dynamprof_sph.vd_xyz{grpind}','- .');hold on;loglog(rs/rvir,dynamprof_sph.Um{grpind}*2/3,': o');
subplot(2,2,4);loglog(rs/rvir,dynamprof_sph.vd_rtf{grpind}','- .');hold on;loglog(rs/rvir,dynamprof_sph.Um{grpind}*2/3,': o');
loglog(rs/rvir,dynamprof_sph.Um{grpind}*2/3+dynamprof_sph.gvd_rtf{grpind}(1,:)','r: s');
loglog(rs/rvir,s(rs/rvir*c,sig2_rs,alpha),'m- o');
figure;semilogx(rs/331,1-(dynamprof_sph.vd_rtf{grpind}(2,:)+dynamprof_sph.vd_rtf{grpind}(3,:))/2./dynamprof_sph.vd_rtf{grpind}(1,:),'- .');
hold on;
fplot(@(x) (-(2.452405824333459*pi*(x+1)^(4/3)*(log(x+1)-x/(x+1)))/(4^(2/3)*x^1.613333333333333)-(.002659291564727156*x^1.386666666666667*(x+1)^(10/3)*(-(145.4021333333334*4^(2/3))/(x^1.386666666666667*(x+1)^(10/3))-(1253.466666666667*4^(2/3))/(x^.3866666666666667*(x+1)^(13/3))))/4^(2/3))/(2),[0.1,20]);
figure;
subplot(2,2,1);semilogx(rs/rvir,dynamprof_sph.gvm_xyz{grpind}','- .');
subplot(2,2,2);semilogx(rs/rvir,dynamprof_sph.gvm_rtf{grpind}','- .');
subplot(2,2,3);loglog(rs/rvir,dynamprof_sph.gvd_xyz{grpind}','- .');hold on;loglog(rs/rvir,dynamprof_sph.Um{grpind}*2/3,': o');
subplot(2,2,4);loglog(rs/rvir,dynamprof_sph.gvd_rtf{grpind}','- .');hold on;loglog(rs/rvir,dynamprof_sph.Um{grpind}*2/3,': o');
figure;semilogx(rs/Rs,1-(dynamprof_sph.gvd_rtf{grpind}(2,:)+dynamprof_sph.gvd_rtf{grpind}(3,:))/2./dynamprof_sph.gvd_rtf{grpind}(1,:),'- .');
figure;loglog(rs/Rs,(2*dynamprof_sph.gvd_rtf{grpind}(1,:)-(dynamprof_sph.gvd_rtf{grpind}(2,:)+dynamprof_sph.gvd_rtf{grpind}(3,:))),'- .');
hold on;
g=@(x,Rhos,Rs) 43007.1*4*pi*Rhos*Rs^3*(log(1+x)-x./(1+x))./x.^2/Rs^2;
loglog(rs/Rs,g(rs/Rs,Rhos,Rs).*rs,'m- o');
loglog(rs/Rs,rs'.*(g(rs/Rs,Rhos,Rs)'+(2*dynamprof_sph.gvd_rtf{grpind}(1,:)-(dynamprof_sph.gvd_rtf{grpind}(2,:)+dynamprof_sph.gvd_rtf{grpind}(3,:)))./rs'),'b- o');

% figure;loglog(rs/Rs,(2*dynamprof_sph.vd_rtf{grpind}(1,:)-(dynamprof_sph.vd_rtf{grpind}(2,:)+dynamprof_sph.vd_rtf{grpind}(3,:))),'- .');
% hold on;
% loglog(rs/Rs,rs.*g(rs/Rs,Rhos,Rs),'m- o');
% loglog(rs/Rs,rs'.*(g(rs/Rs,Rhos,Rs)'+(2*dynamprof_sph.vd_rtf{grpind}(1,:)-(dynamprof_sph.vd_rtf{grpind}(2,:)+dynamprof_sph.vd_rtf{grpind}(3,:)))./rs'),'b- o');

%%
beta=zeros(100,1);
for grpind=1:100
    beta(grpind)=polyfit(haloprof.rs{grpind},2*dynamprof.Ums{grpind}./(dynamprof.vsms{grpind}-dynamprof.vm2s{grpind}),0);
end
hist(beta(haloflag(1:100,1)),0:0.1:2);
