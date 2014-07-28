
RunNum='HY300';SnapLoad=59;
datadir=['/mnt/uv/HBT/data/',RunNum,'/subcat/anal/image/'];
halo=readhalo_size([datadir,'../../profile/logbin'],SnapLoad,'halo');
%
RunNum2='HY300Nd20012';SnapLoad=59;
datadir=['/mnt/uv/HBT/data/',RunNum2,'/subcat/anal/image/'];
halo2=readhalo_size([datadir,'../../profile/logbin'],SnapLoad,'halo');
simu300.halo=halo;
simu300N.halo=halo2;
%%
RunNum='MilliMill';SnapLoad=63;
datadir=['/mnt/charon/HBT/data/',RunNum,'/subcat/'];
halo=readhalo_size([datadir,'/profile/logbin'],SnapLoad,'halo');
%%

%%
save /work/Projects/HBT/code/data/ELUCID.mat simu300 simu300N simu3002 simu3002N 
%%
load /work/Projects/HBT/code/data/ELUCID.mat
%%
dndM=importdata('/data/jxhan/Downloads/mass_functions_010.dat',' ',1);%tinker
dndM2=importdata('/data/jxhan/Downloads/mass_functions_006.dat',' ',1);%ST
dndM3=importdata('/data/jxhan/Downloads/mass_functions_009.dat',' ',1);%tinker 200c
%%
x=logspace(1,6,40);
dlnm=(log(x(2))-log(x(1)));
myfigure;
set(gcf,'DefaultLineMarkerSize',10,'DefaultLineLineWidth',1);
[xm,ym]=loghist(simu300.halo.Mvir(:,3)*mp,x);
xm=xm*1e10; cm=1/(300^3)/dlnm*xm;
h1=plot(xm,ym.*cm,'r-','linewidth',2);
ym0=ym;
hold on;
[xm,ym]=loghist(simu300N.halo.Mvir(:,3)*mp,x);
xm=xm*1e10; cm=1/(300^3)/dlnm*xm;
h2=ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'ro','logxy','hhy',0.5);
yratN=ym./ym0; eratN=sqrt(1./ym+1./ym0);

[xm,ym]=loghist(simu3002.halo.Mvir(:,3)*mp,x);
xm=xm*1e10; cm=1/(300^3)/dlnm*xm;
h3=plot(xm,ym.*cm,'g-','linewidth',2);
yrat2=ym./ym0;erat2=sqrt(1./ym+1./ym0);
hold on;
[xm,ym]=loghist(simu3002N.halo.Mvir(:,3)*mp,x);
xm=xm*1e10; cm=1/(300^3)/dlnm*xm;
h4=ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'gs','logxy','hhy',0.5);
yrat2N=ym./ym0;erat2N=sqrt(1./ym+1./ym0);
h5=plot(dndM.data(:,1),dndM.data(:,2),'k--','linewidth',1);
h6=plot(dndM2.data(:,1),dndM2.data(:,2),'k-','linewidth',1);
% plot(dndM3.data(:,1),dndM3.data(:,2),'k-');
xlim([1e12,2e15]);
ylim([2e8,1e10]);
xscale('log');yscale('log');
l=legend([h1,h2(1),h3,h4(1),h5,h6],'300','300N','3002','3002N','Tinker08','SMT01');
set(l,'location','southwest');
xlabel('$M_{200b}/[M_\odot/h]$');ylabel('$Mdn/d\ln M$');
outputdir='/work/Projects/HBT/code/data/show/images/show';
print('-depsc',[outputdir,'/msfunCmp.eps']);
%%
x=logspace(1,6,50);
dlnm=(log(x(2))-log(x(1)));
figure;
[xm,ym]=loghist(simu300.halo.Mvir(:,3),x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'ro','logxy');
ym0=ym;
hold on;
[xm,ym]=loghist(simu300N.halo.Mvir(:,3),x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'gs','logxy');
yratN=ym./ym0; eratN=sqrt(1./ym+1./ym0);
plot(dndM.data(:,1),dndM.data(:,3),'k-');
plot(dndM2.data(:,1),dndM2.data(:,3),'k-');
xlim([1e12,1e15]);

[xm,ym]=loghist(simu3002.halo.Mvir(:,3),x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'r-','logxy');
yrat2=ym./ym0;erat2=sqrt(1./ym+1./ym0);
hold on;
[xm,ym]=loghist(simu3002N.halo.Mvir(:,3),x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'g-','logxy');
yrat2N=ym./ym0;erat2N=sqrt(1./ym+1./ym0);
%%
mp=1.440505;
getmass=@(n) n*mp; %n.*(1-n.^-0.6)*mp;
x=logspace(1,6,50);
dlnm=(log(x(2))-log(x(1)));
figure;
[xm,ym]=loghist(getmass(simu300.halo.mass),x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'r-','logxy');
ym0=ym;
hold on;
[xm,ym]=loghist(getmass(simu300N.halo.mass),x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'ro','logxy');
yratN=ym./ym0; eratN=sqrt(1./ym+1./ym0);
plot(dndM.data(:,1),dndM.data(:,3),'k-')
plot(dndM2.data(:,1),dndM2.data(:,3),'k-');
xlim([1e12,1e15]);

[xm,ym]=loghist(getmass(simu3002.halo.mass),x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'g-','logxy');
yrat2=ym./ym0;erat2=sqrt(1./ym+1./ym0);
hold on;
[xm,ym]=loghist(getmass(simu3002N.halo.mass),x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'gs','logxy');
yrat2N=ym./ym0;erat2N=sqrt(1./ym+1./ym0);
%%
mp=1.440505;
x=logspace(1,6,50);
dlnm=(log(x(2))-log(x(1)));
figure;
[xm,ym]=loghist(simu300.halo.mass.*(1-simu300.halo.mass.^-0.6)*mp,x);
xm=xm*1e10; cm=1/(300^3)/dlnm*xm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'r-','logxy');
ym0=ym;
hold on;
[xm,ym]=loghist(simu300N.halo.mass.*(1-simu300N.halo.mass.^-0.6)*mp,x);
xm=xm*1e10; cm=1/(300^3)/dlnm*xm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'ro','logxy');
yratN=ym./ym0; eratN=sqrt(1./ym+1./ym0);
plot(dndM.data(:,1),dndM.data(:,2),'k--')
plot(dndM2.data(:,1),dndM2.data(:,2),'k-');
xlim([1e12,1e15]);

[xm,ym]=loghist(simu3002.halo.mass.*(1-simu3002.halo.mass.^-0.6)*mp,x);
xm=xm*1e10; cm=1/(300^3)/dlnm*xm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'g-','logxy');
yrat2=ym./ym0;erat2=sqrt(1./ym+1./ym0);
hold on;
[xm,ym]=loghist(simu3002N.halo.mass.*(1-simu3002N.halo.mass.^-0.6)*mp,x);
xm=xm*1e10; cm=1/(300^3)/dlnm*xm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'gs','logxy');
yrat2N=ym./ym0;erat2N=sqrt(1./ym+1./ym0);
%%
myfigure;
h1=ploterr(xm,yrat2,[],[],'k-','logx');hold on;
h2=ploterr(xm,yratN,[],eratN,'ro','logx');set(h2(1),'markersize',8);
h3=ploterr(xm,yrat2N,[],[],'g--','logx');set(h3(1),'markersize',15);
plot(xm,ones(size(xm)),'k:');
l=legend([h1(1),h2(1),h3(1)],'3002/300','300N/300','3002N/300');set(l,'location','northwest');
xlabel('$M_{200b}/[M_\odot/h]$');ylabel('$dn/dn_0$');
xlim([1e11,1e15]);
outputdir='/work/Projects/HBT/code/data/show/images/show';
print('-depsc',[outputdir,'/msfunRat.eps']);
%%
x=logspace(1,6,50);
dlnm=log(x(2))-log(x(1));
figure;
[xm,ym]=loghist(simu300.halo.mass,x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'ro','logxy');
hold on;
[xm,ym]=loghist(simu300N.halo.mass,x);
xm=xm*1e10; cm=1/(300^3)/dlnm;
ploterr(xm,ym.*cm,[],sqrt(ym).*cm,'gs','logxy');
plot(dndM2.data(:,1),dndM2.data(:,3),'k-');
%%
myfigure;
for i=1:sum(halo.Rvir(:,1)>1)
h=plot_circle(halo.Cen(i,:),halo.Rvir(i,1)^3);set(h,'color','r'); hold on;
end


for i=1:sum(halo2.Rvir(:,1)>1)
h=plot_circle(halo2.Cen(i,:),halo2.Rvir(i,1)^3);set(h,'line','--','color','g'); hold on;
end

xlim([0,300]);ylim([0,300]);
title('L300(red) vs. L300Nd(green)');
outputdir='/work/Projects/HBT/code/data/show/images';
print('-depsc',[outputdir,'/halos',RunNum,'.eps']);