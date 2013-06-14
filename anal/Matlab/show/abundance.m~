% RunNum=6113;
% Nsnap=99:-10:59;

RunNum=6120;
Nsnap=99:-10:49;
MinSubMass=1e-4;
nbin=5;
outputdir=['/home/kam/Projects/HBT/code/data/show/massfun/',num2str(RunNum),'/'];
datadir=['/mnt/A4700/data/',num2str(RunNum),'/subcat/anal/massfun/'];
n=numel(Nsnap);
scaleF=load_scaleF(num2str(RunNum));
A=cell(n,1);
for i=1:n
file=[datadir,'abundanceN_',num2str(Nsnap(i),'%03d'),'.',num2str(MinSubMass,'%e')];
% file=[datadir,'abundance_',num2str(Nsnap(i),'%03d')];
A{i}=load_abundance(file);
A{i}.rs=A{i}.rs;%*scaleF(Nsnap(i)+1);
A{i}.rhos=A{i}.rhos;%*scaleF(Nsnap(i)+1)^-3;
end

figure;
x=[];y=[];z=[];
for i=1:n
% errorbar(A{i}.Mvir,A{i}.Ns,sqrt(A{i}.Ns),colors{i});set(gca,'xscale','log
% ','yscale','log');hold on;
f=A{i}.Ns>50;
x=[x;log10(A{i}.Mvir(f))];
y=[y;log10(A{i}.rs(f))];
z=[z;log10(A{i}.Ns(f))];
end
plot3(x,y,z,'.');hold on;
xlabel('M');ylabel('rs');
[xx,yy]=meshgrid([3:0.1:5],[1:0.1:3]);
co=[x,y,ones(size(x))]\z;
zz=co(1)*xx+co(2)*yy+co(3);
mesh(xx,yy,zz);

colors={'r.';'g.';'b.';'k.';'c.';'m.'};
colorl={'r.-';'g.-';'b.-';'k.-';'c.-';'m.-'};
figure;
for i=1:n
% errorbar(A{i}.Mvir,A{i}.Ns,sqrt(A{i}.Ns),colors{i});set(gca,'xscale','log','yscale','log');
% plot(A{i}.Mvir,A{i}.Ns,colors{i});hold on;
[xm,ym]=logbin(A{i}.Mvir,A{i}.Ns,nbin);
semilogx(xm,ym,colorl{i});
hold on;
end
set(gca,'xscale','log');
% set(gca,'yscale','log');

figure;
for i=1:n
% errorbar(A{i}.rhos,A{i}.Ns,sqrt(A{i}.Ns),colors{i});set(gca,'xscale','log','yscale','log');
% plot(A{i}.rhos,A{i}.Ns,colors{i});hold on;
[xm,ym]=logbin(A{i}.rhos,A{i}.Ns,nbin);
semilogx(xm,ym,colorl{i});
hold on;
end
set(gca,'xscale','log');
% set(gca,'yscale','log');

figure;
for i=1:n
% errorbar(A{i}.rs,A{i}.Ns,sqrt(A{i}.Ns),colors{i});set(gca,'xscale','log','yscale','log');
% plot(A{i}.rs,A{i}.Ns,colors{i});hold on;
[xm,ym]=logbin(A{i}.rs,A{i}.Ns,nbin);
semilogx(xm,ym,colorl{i});
hold on;
end
set(gca,'xscale','log');
% set(gca,'yscale','log');

figure;
for i=1:n
% errorbar(A{i}.c,A{i}.Ns,sqrt(A{i}.Ns),colors{i});set(gca,'xscale','log','yscale','log');
% plot(A{i}.c,A{i}.Ns,colors{i});hold on;
[xm,ym]=logbin(A{i}.c,A{i}.Ns,nbin);
semilogx(xm,ym,colorl{i});
hold on;
end
set(gca,'xscale','log');
% set(gca,'yscale','log');

figure;
for i=1:n
% errorbar(A{i}.Rvir,A{i}.Ns,sqrt(A{i}.Ns),colors{i});set(gca,'xscale','log','yscale','log');
% plot(A{i}.Rvir,A{i}.Ns,colors{i});hold on;
[xm,ym]=logbin(A{i}.Rvir,A{i}.Ns,nbin);
semilogx(xm,ym,colorl{i});
hold on;
end
set(gca,'xscale','log');
% set(gca,'yscale','log');