%% data preparation
% to add: 6702low, cosmo

outputdir='/home/kam/Projects/HBT/code/data/show/massfun/resim';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show/massfun/resim';
addpath(genpath('../post'));

virtype=1;
rmin=0.01;rmax=1;
Nsnap=162;grpid=0;
RunName='AqA2';softeninghalo=4.8e-5;mp=1.000007e-06;
% softeninghalo=1.5;
% Nsnap=1023;grpid=0;
% RunName='AqA4';softeninghalo=0.25;mp=2.8679e-5;

scaleF_file=['/mnt/charon/HBT/data/',RunName,'/subcatmore/Redshift.dat'];
% scaleF_file=['/mnt/A4700/data/',RunName,'/subcat/Redshift.dat'];
a=load_scaleF(scaleF_file);
z=1/a(Nsnap+1)-1;
% mp=0.0104089; %6702DM
% mp=2.8679e-5;%AqA4

% datadir=['/mnt/A4700/data/',RunName,'/subcat/anal/massfun/'];
datadir=['/mnt/charon/HBT/data/',RunName,'/subcatmore/anal/massfun/'];
[Mlist1,M1,R1,Rlist1,Sublist1]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','');
% datadir=['/mnt/A4700/data/',RunName,'/subcatS/anal/massfun/'];
% datadir=['/mnt/charon/HBT/data/',RunName,'/subcat/anal/subfind/massfun/'];
% [Mlist2,M2,R2,Rlist2]=Msublist_in_radii(datadir,virtype,Nsnap,grpid,rmin,rmax,'bindata','');
%%
datadir=['/mnt/charon/HBT/data/',RunName,'/subcatmore/profile/logbin/'];
% datadir=['/mnt/A4700/data/',RunName,'/subcat/profile/logbin/'];
halo=readhalo_size(datadir,Nsnap,'halo');
haloprof=readhalo_prof_single(datadir,Nsnap,'halo',halo,grpid+1);
vs=diff(logspace(log10(max(softeninghalo,halo.rmax(grpid+1)*1e-2)),log10(halo.rmax(grpid+1)),halo.nbin(grpid+1)+1).^3)';
rho=(haloprof.ns+haloprof.no+haloprof.nb)./vs;
r=haloprof.rs;
erho=sqrt(rho)./sqrt(vs);
%%
% Mass Profile
myfigure;
h1=loglog(r/R1,rho/(M1/mp/R1^3),'k-','displayname','Halo'); hold on;
loglog(r/R1, rho/(M1/mp/R1^3).*(r/R1).^(1.5), 'g-', 'displayname', 'sat Model');
nbin=15;
nmin=2000;
rbin=logspace(-2,0.7,nbin)*R1;
% f=Mlist1>nmin*mp;
f=Mlist1>10000;%&Mlist1<10000000;
[xr1,n1,rhon1]=loghist(Rlist1(f),rbin,[],[],1,3);
erhon1=sqrt(n1)./diff(rbin.^3);
rhomean=sum(n1)/R1^3;
h4=errorbar(xr1/R1,rhon1/rhomean,erhon1/rhomean,'-','linewidth',3,'color',[0.8,0.8,1], 'displayname','sat Infall');
% f=Mlist2>nmin*mp;
% [xr2,n2,rhon2]=loghist(Rlist2(f),rbin,[],[],1,3);
% hold on;
% h5=plot(xr2/R1,rhon2/rhomean,'r-','linewidth',3,'color',[1,0.8,0.8],'displayname','SUBFIND');
xlabel('D/Rvir');ylabel('$n/<n>$');
set(gca,'xscale','log','yscale','log');
legend('show')
% print('-depsc','/work/Projects/SubProf/plots/A2subprof_infall.eps')
%%
hold on;
reinasto=199e-3;
alpha=0.678;
neinasto=einasto(xr1/reinasto,alpha);
nnorm=quad(@(r) einasto(r/reinasto,alpha).*r.^2*4*pi,0,R1)/(R1^3);
plot(xr1/R1,neinasto/nnorm,'r');
%%
figure();
loglog(r/R1, einasto(r/reinasto,alpha)./rho);
%%
tmp=importdata('/mnt/charon/HBT/data/AqA2/subcatmore/anal/steller/SnapInfall_162',',',1);
submass=tmp.data(:,5);%current mass
directinfall=tmp.data(:,end);
submass=submass(Sublist1+1);
directinfall=directinfall(Sublist1+1);
f=Mlist1>1000;
semilogx(Rlist1(f)/R1, submass(f)./Mlist1(f),'r.');

%%

nbin=10;
nmin=1000;
rbin=logspace(-2,0.7,nbin)*R1;
f=Mlist1>nmin*mp;
[xr1,m1,rho1]=loghist(Rlist1(f),rbin,[],[],Mlist1(f),3);
[xr1,n1,rhon1]=loghist(Rlist1(f),rbin,[],[],1,3);
erhon1=sqrt(n1)./diff(rbin.^3);
f=Mlist2>nmin*mp;
[xr2,m2,rho2]=loghist(Rlist2(f),rbin,[],[],Mlist2(f),3);
[xr2,n2,rhon2]=loghist(Rlist2(f),rbin,[],[],1,3);
erhon2=sqrt(n2)./diff(rbin.^3);

hold on;
rhomean=sum(n1)/R1^3;
h2=ploterr(xr1/R1,rhon1/rhomean,[],erhon1/rhomean,'o-','logxy');
hold on;
h3=ploterr(xr2/R1,rhon2/rhomean,[],erhon2/rhomean,'ro-','logxy');

f=@(x) exp(-2/0.678*((x/0.7).^0.678-1));
plot(xr2/R1,f(xr2/R1),'r-');




xlabel('D/Rvir');ylabel('$n/<n>$');
legend([h1,h2(1),h3(1),h4,h5],'\rho/<\rho>','HBT,m>100','SUBFIND,m>100','HBT,m>20','SUBFIND,m>20');
% hl=legend('show','location','southwest');set(hl,'interpreter','latex');
% title(['z=',num2str(z,'%2.1f')]);
% print('-depsc',[outputdir,'/numberprof_',RunName,'.eps']);
% hgsave([outputdir,'/msfun_',RunName,'.fig']);

%% average mass
myfigure;
% 
loglog(xr1/R1,m1./n1/M1,':','displayname','HBT,<m>');
hold on;
loglog(xr2/R1,m2./n2/M1,'r:','displayname','SUBFIND,<m>');
    
set(gca,'yscale','linear','yminortick','on');
xlabel('D/Rvir');ylabel('$<m>/Mvir$');
% hl=legend('show','location','southwest');set(hl,'interpreter','latex');
title(['z=',num2str(z,'%2.1f')]);
% print('-depsc',[outputdir,'/msfun_',RunName,'.eps']);
% hgsave([outputdir,'/msfun_',RunName,'.fig']);
figure;
semilogx(xr1/R1,m1./n1./(m2./n2),'o-');
%% radial distr for different mass bins
mbin=10.^([log10(100*mp),-3,-2,-1,0]);
[~,ibin]=histc(Mlist2,mbin);
nbin=6;
rbin=logspace(-2,0.7,nbin)*R1;
lines=['r-';'g-';'b-';'c-';'m-';'k-';];
h1=[];
figure;
for i=1:numel(mbin)-1
f=(ibin==i);
% [xr1,m1,rho1]=loghist(Rlist1(f),rbin,[],[],Mlist1(f),3);
[xr1,n1,rhon1]=loghist(Rlist2(f),rbin,[],[],1,3);
erhon1=sqrt(n1)./diff(rbin.^3);
rhomean=sum(n1)/R1^3;
h=ploterr(xr1/R1,rhon1/rhomean,[],erhon1/rhomean,lines(i,:),'logxy');
h1=[h1,h];
hold on;
end
legend(h1(1,:),num2str((1:numel(mbin)-1)'))
% print('-depsc',[outputdir,'/ra
%%
mbin=10.^([log10(100*mp),-3.5,0]);
[~,ibin]=histc(Mlist2,mbin);
nbin=6;
rbin=logspace(-2,0,nbin)*R1;
lines=['r--';'g--';'b--';'c--';'m--';'k--';];
h2=[];
% figure;
for i=1:numel(mbin)-1
f=(ibin==i);
% [xr1,m1,rho1]=loghist(Rlist2(f),rbin,[],[],Mlist1(f),3);
[xr1,n1,rhon1]=loghist(Rlist2(f),rbin,[],[],1,3);
erhon1=sqrt(n1)./diff(rbin.^3);
rhomean=sum(n1)/R1^3;
h=ploterr(xr1/R1,rhon1/rhomean,[],erhon1/rhomean,lines(i,:),'logxy');
h2=[h2,h];
hold on;
end
legend(h2(1,:),num2str((1:numel(mbin-1))'))
%% shape profile
RunName='6402DM';
Nsnap=99;
file=['/mnt/A4700/data/',RunName,'/subcat/anal/sub_shape_',num2str(Nsnap)];
Q=load(file);
Q=sort(Q,2);
% Q=Q(logical(sum(abs(Q),2)),:);%exclude Q=0;
file=['/mnt/A4700/data/',RunName,'/subcatS/anal/sub_shape_',num2str(Nsnap)];
Qs=load(file);
Qs=sort(Qs,2);
% Qs=Qs(logical(sum(abs(Qs),2)),:);%exclude Q=0;

rQ=Q(:,end)./Q(:,1);
irQ=1./rQ;
fQ=std(Q,0,2)./mean(Q,2);

rQs=Qs(:,end)./Qs(:,1);
irQs=1./rQs;
fQs=std(Qs,0,2)./mean(Qs,2);
%%
virtype=0;
grpid=0;

  datadir=['/mnt/A4700/data/',RunName,'/subcat/anal/massfun/'];
        submass=read_submass([datadir,'submass_',num2str(Nsnap,'%03d')]);
        subcom=read_subcom([datadir,'subcom_',num2str(Nsnap,'%03d')]);
        massdata=[submass,subcom];
        clear submass subcom
        cenid=read_cid([datadir,'cid_',num2str(Nsnap,'%03d')]);
        switch virtype
            case 0
                grpsize=read_grpsize([datadir,'grpsizeVIR_',num2str(Nsnap,'%03d')]);
            case 1
                grpsize=read_grpsize([datadir,'grpsizeC200_',num2str(Nsnap,'%03d')]);
            case 2
                grpsize=read_grpsize([datadir,'grpsizeB200_',num2str(Nsnap,'%03d')]);
        end
        Mhost=grpsize(grpid+1,1);
        Rvir=grpsize(grpid+1,2);
        clear grpsize
        center_id=cenid(grpid+1)+1;
        Rsub=sqrt(sum((massdata(:,2:4)-repmat(massdata(center_id,2:4),size(massdata,1),1)).^2,2));
        Msub=massdata(:,1);
  
  %SUBFIND
  datadir=['/mnt/A4700/data/',RunName,'/subcatS/anal/massfun/'];
        submass=read_submass([datadir,'submass_',num2str(Nsnap,'%03d')]);
        subcom=read_subcom([datadir,'subcom_',num2str(Nsnap,'%03d')]);
        massdata=[submass,subcom];
        clear submass subcom
        cenid=read_cid([datadir,'cid_',num2str(Nsnap,'%03d')]);
        center_id=cenid(grpid+1)+1;
        RsubS=sqrt(sum((massdata(:,2:4)-repmat(massdata(center_id,2:4),size(massdata,1),1)).^2,2));
        MsubS=massdata(:,1);
%%        
nbin=10;
rbin=logspace(-2,0.7,nbin)*Rvir;
f1=logical(sum(abs(Q),2));
f2=logical(sum(abs(Qs),2));
val=sqrt(irQ);
valS=sqrt(irQs);
[xr1,m1,rho1]=loghist(Rsub(f1),rbin,[],[],val(f1),3);
[xr1,n1,rhon1]=loghist(Rsub(f1),rbin,[],[],1,3);
[xr2,m2,rho2]=loghist(RsubS(f2),rbin,[],[],valS(f2),3);
[xr2,n2,rhon2]=loghist(RsubS(f2),rbin,[],[],1,3);

figure;
for i=1:numel(f1)
    if f1(i)
    semilogx(Rsub(i)/Rvir,val(i),'o','markersize',Msub(i).^(1/2)*20); hold on;
    end
end
h1=plot(-1,-1,'o');
% semilogx(Rsub(f1)/Rvir,val(f1),'.');
hold on;
h2=semilogx(xr1/Rvir,m1./n1,'r');
h3=semilogx(xr2/Rvir,m2./n2,'g');
l=legend([h1,h2,h3],'HBT','<HBT>','<SUBFIND>');
set(l,'location','southwest');
xlabel('R/Rvir');ylabel('sqrt(I_{min}/I_{max})');
print('-depsc',['/work/Projects/HBT/code/data/show/shape',num2str(size(Q,2)),'.',RunName,'_HBT.eps']);

figure;
for i=1:numel(f2)
    if f2(i)
    semilogx(RsubS(i)/Rvir,val(i),'o','markersize',MsubS(i).^(1/2)*20); hold on;
    end
end
h1=plot(-1,-1,'o');
% semilogx(RsubS(f2)/Rvir,valS(f2),'gs');
hold on;
h2=semilogx(xr1/Rvir,m1./n1,'r');
h3=semilogx(xr2/Rvir,m2./n2,'g');
l=legend([h1,h2,h3],'SUBFIND','<HBT>','<SUBFIND>');
xlabel('R/Rvir');ylabel('sqrt(I_{min}/I_{max})');
set(l,'location','southwest');
print('-depsc',['/work/Projects/HBT/code/data/show/shape',num2str(size(Q,2)),'.',RunName,'_SUBFIND.eps']);
