clear;
addpath('../post');

runnum='8213';Nsnap=59;
% runnum='6702DM';Nsnap=99;

file=['/mnt/A4700/data/',runnum,'/subcat/anal/sub_shape_',num2str(Nsnap)];
Q=load(file);
Q=sort(Q,2);
Q=Q(logical(sum(abs(Q),2)),:);%exclude Q=0;
file=['/mnt/A4700/data/',runnum,'/subcatS/anal/sub_shape_',num2str(Nsnap)];
Qs=load(file);
Qs=sort(Qs,2);
Qs=Qs(logical(sum(abs(Qs),2)),:);%exclude Q=0;

% Q=sqrt(Q);
% Qs=sqrt(Qs);

rc=20;
fc=1.0;

r=Q(:,3)./Q(:,1);
ir=1./r;
f=std(Q,0,2)./mean(Q,2);
[nf,xf]=hist(f,0:0.1:1.5);
[nr,xr]=hist(r,1:0.4:10);
[nir,xir]=hist(ir,0:0.1:1);
% partmass=0.008848;
% load ../massfun/v66/sub_mass_099
% idr=[sub_mass_099(find(r>rc))/partmass,r(find(r>rc)),find(r>rc)];
% idf=[sub_mass_099(find(f>fc))/partmass,f(find(f>fc)),r(find(f>fc)),find(f>fc)];

rs=Qs(:,3)./Qs(:,1);
irs=1./rs;
fs=std(Qs,0,2)./mean(Qs,2);
[nfs,xfs]=hist(fs,0:0.1:1.5);
[nrs,xrs]=hist(rs,1:0.4:10);
[nirs,xirs]=hist(irs,0:0.1:1);
% load ../massfun/subfind/SubFind_mass_099
% idrs=[SubFind_mass_099(find(rs>rc)), rs(find(rs>rc)),find(rs>rc)];
% idfs=[SubFind_mass_099(find(fs>fc)), fs(find(fs>fc)),rs(find(fs>fc)),find(fs>fc)];
%%
% fun=@plot;
fun=@stairs;

myfigure;
axes('position',[.2,.15,0.32,0.8]);
fun(xf,nf,'r--');hold on;fun(xfs,nfs,'k-');
xlabel('$\zeta$');
ylabel('Counts');
xlim([0,1.9]);
set(gca,'yscale','log');
l=legend('HBT','SUBFIND');
set(l,'fontsize',15);
% title('differential distribution in principle axes of moments of inertia');

axes('position',[0.53,0.15,0.32,0.8]);
% fun(xr,nr,'r--');hold on;fun(xrs,nrs,'k-');
% xlabel('$\gamma$');
% ylabel('Counts');set(gca,'yaxislocation','right');
% xlim([0,10]);
% set(gca,'yscale','log');

fun(xir,nir,'r--');hold on;fun(xirs,nirs,'k-');
xlabel('$\gamma$');
ylabel('Counts');set(gca,'yaxislocation','right');
xlim([0,1]);
set(gca,'yscale','log');

% legend('HBT','SUBFIND');
% set(gcf,'paperposition',[0.6,6,20,20]);
outputdir='/home/kam/Projects/HBT/code/data/show';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/subshape_',runnum,'_',num2str(Nsnap),'.new.eps']);