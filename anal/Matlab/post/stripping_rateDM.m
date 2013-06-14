clear;
global G a HUBBLE0 Omega0 OmegaLambda
G=43007.1;HUBBLE0=0.1;Omega0=0.3;OmegaLambda=0.7;

runnum=6114;
DMpmass=zeros(size(runnum));
% runnum=[6700,6701,6702,6500,6501,6506,6402,6404,6409,6600,6601,6602];
% DMpmass=[0.0130057,0.00834475,0.008848,7.44348e-05,5.8378e-5,6.45321e-05,...
%     5.45271e-6,5.29086e-06,5.23311e-06,0.000789798,0.000670876,0.000546019];

param_gas=[];snap_gas=[];Mrate=[];Mhost=[];K=[];j2=[];Chost=[];Csat=[];hist_gas=[];pmass_gas=[];
param_dm=[];snap_dm=[];MrateD=[];MhostD=[];KD=[];j2D=[];ChostD=[];CsatD=[];hist_dm=[];pmass_dm=[];
param_both=[];
rmin=[];
for i=1:numel(runnum)
disp(['Simulation',num2str(runnum(i))])
   
history_file=['/mnt/A4700/data/',num2str(runnum(i)),'/subcat/anal/history_000_005.dat'];
history_par_file=['/mnt/A4700/data/',num2str(runnum(i)),'/subcat/anal/historypar_000_005.dat'];

disp('loading files...')
scaleF_file=['/mnt/A4700/data/',num2str(runnum(i)),'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);    
[hist_bin,DMpmass(i)]=load_hist_bin(history_file);
pmass_bin=repmat(DMpmass(i),size(hist_bin));
[Snappar,Mratepar,Mhostpar,Rhostpar,Kpar,j2par,Kerrpar,j2errpar,Chostpar,Csatpar]=load_hist_par(history_par_file);

% rmin=[rmin;radial_range(hist_bin,pmass_bin,Snappar)];

disp('fitting DM...')
par_dm=DM_stripping_rate_stat(hist_bin,pmass_bin,Snappar,Chostpar,Csatpar);
filterD=logical(par_dm(:,3)>0&Csatpar(:,1)>0&Csatpar(:,2)>0&Chostpar(:,1)>0&Chostpar(:,2)>0);
param_dm=[param_dm;par_dm(filterD,:)];%alpha,beta,xi
snap_dm=[snap_dm;Snappar(filterD,:)];
MrateD=[MrateD;Mratepar(filterD,:)];
MhostD=[MhostD;Mhostpar(filterD,:)*DMpmass(i)];
KD=[KD;Kpar(filterD,:)];
j2D=[j2D;j2par(filterD,:)];
ChostD=[ChostD;Chostpar(filterD,:)];
CsatD=[CsatD;Csatpar(filterD,:)];
hist_dm=[hist_dm;hist_bin(filterD)];
pmass_dm=[pmass_dm;pmass_bin(filterD)];
disp(['N=',num2str(numel(find(filterD)))]);
% param_both=[param_both;par_gas(filter&filterD,:),par_dm(filter&filterD,:)];

clear hist_bin Snappar Mratepar Mhostpar Rhostpar Kpar j2par Kerrpar j2errpar Chostpar Csatpar par_gas filter
end
par=param_dm(:,3);
figure;hist(par,20);xlabel('<\Delta f_{dm}>');ylabel('Counts');
x=linspace(0,0.4,40);
figure;h1=subplot(2,1,1);h2=subplot(2,1,2);
filter=logical(par<0.1);
y=histc(par(filter),x);
subplot(h1);plot(x,y,'k- .');hold on;
subplot(h2);plot(x,cumsum(y),'k- .');

% x=logspace(-2,1,50);
x=0.1:0.05:5;
figure;h1=subplot(2,1,1);h2=subplot(2,1,2);
% filter=logical(param_dm(:,3)<0.1);
filter=logical(param_dm(:,1)<1.2&param_dm(:,1)>0.8);
y=histc(param_dm(filter,1),x);
subplot(h1);plot(x+0.025,y,'k- .');xlabel('\alpha_{dm}');ylabel('Counts');set(gca,'xlim',[0.1,3]);
x=0.1:0.04:1;
% x=logspace(-1,1,50);
% y=histc(param_dm(filter,2).*param_dm(filter,1),x);
y=histc(param_dm(filter,2),x);
subplot(h2);plot(x+0.02,y,'k- .');xlabel('\beta_{dm}');ylabel('Counts');

Msat=Mrate.*Mhost;MsatD=MrateD.*MhostD;
e2=1-j2;e2D=1-j2D;
disp('all done')
figure('name','parameter distributions');
subplot(2,2,3);hist(param_dm(:,1));xlabel('\alpha_{dm}');subplot(2,2,4);hist(param_dm(:,2));xlabel('\beta_{dm}');

figure('name','Mass ratio dependence');
subplot(2,2,3);loglog(MrateD(:,2),param_dm(:,1),'.');ylabel('\alpha_{dm}');xlabel('m/M');
subplot(2,2,4);loglog(MrateD(:,2),param_dm(:,2),'.');ylabel('\beta_{dm}');xlabel('m/M');

figure('name','Host Mass dependence');
subplot(2,2,3);loglog(MhostD(:,2),param_dm(:,1),'.');ylabel('\alpha_{dm}');xlabel('M');
subplot(2,2,4);loglog(MhostD(:,2),param_dm(:,2),'.');ylabel('\beta_{dm}');xlabel('M');

figure('name','Sat Mass dependence');
subplot(2,2,3);loglog(MsatD(:,2),param_dm(:,1),'.');ylabel('\alpha_{dm}');xlabel('m');
subplot(2,2,4);loglog(MsatD(:,2),param_dm(:,2),'.');ylabel('\beta_{dm}');xlabel('m');

figure('name','Host concentration dependence');
subplot(2,2,3);loglog(ChostD(:,2),param_dm(:,1),'.');ylabel('\alpha_{dm}');xlabel('Chost');
subplot(2,2,4);loglog(ChostD(:,2),param_dm(:,2),'.');ylabel('\beta_{dm}');xlabel('Chost');

figure('name','Sat concentration dependence');
subplot(2,2,3);loglog(CsatD(:,1),param_dm(:,1),'.');ylabel('\alpha_{dm}');xlabel('Csat');
subplot(2,2,4);loglog(CsatD(:,1),param_dm(:,2),'.');ylabel('\beta_{dm}');xlabel('Csat');
x=CsatD(:,1);y=param_dm(:,2);x=x(y>0);y=y(y>0);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

figure('name','orbit shape dependence');
subplot(2,2,3);loglog(e2D(:,2),param_dm(:,1),'.');ylabel('\alpha_{dm}');xlabel('e2');
subplot(2,2,4);loglog(e2D(:,2),param_dm(:,2),'.');ylabel('\beta_{dm}');xlabel('e2');

figure('name','orbital energy dependence');
subplot(2,2,3);loglog(KD(:,2),param_dm(:,1),'.');ylabel('\alpha_{dm}');xlabel('K');
subplot(2,2,4);loglog(KD(:,2),param_dm(:,2),'.');ylabel('\beta_{dm}');xlabel('K');

% figure;hist(log10(rmin));

%%
Ndm=size(hist_dm,1);
hand=figure();hand2=figure();

for h=1:Ndm
[Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_dm{h},pmass_dm(h),snap_dm(h,2),1,1,1);
if Ninfall>5
% alpha=0;beta=CsatD(h,1)^-0.4;
par=strip_coeff(ChostD(h,2),CsatD(h,2));
% flag_redistr=0;alpha=1;beta=0.45;p=[par(1)/alpha,par(2)/alpha/beta];
% flag_resitr=0;alpha=0;beta=0;p=[par(1)/param_dm(h,1),par(2)/param_dm(h,1)/param_dm(h,2)];
flag_redistr=-6;alpha=0.6;beta=0.4;p=[alpha,beta];
% flag_redistr=-6;alpha=0;beta=0;p=param_dm(h,1:2);
% flag_redistr=-5;alpha=1;beta=0.5;p=[alpha,beta];
% flag_redistr=-5;alpha=0;beta=0;p=param_dm(h,1:2);
figure(hand);
xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,p,3,flag_redistr);
figure(hand2);
xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,p,2,flag_redistr);
end
end 

figure(hand);
subplot(2,1,1);axis([0,0.3,0,1]);hold off;xlabel('\Deltaln(a)');ylabel('Mdm/Mdm(0)');title('simulation,z_{infall}>1,m/M~[0.01,0.1]')%set(gca,'yscale','log');
subplot(2,1,2);axis([0,0.3,0,1]);hold off;xlabel('\Deltaln(a)');ylabel('Mdm/Mdm(0)');title('Best fit Model: $df/dln(a)=\frac{\sqrt\Delta}{2\pi}[\alpha{-}\beta(f\frac{Rvir}{R})^{2/3}]$','interpreter','latex','FontSize',16);%set(gca,'yscale','log');
figure(hand2);
plot([0,1],[0,1],'r');axis([0,1,0,1]);hold off;xlabel('simulation');ylabel('model');title(['dm,\alpha=',num2str(alpha),',\beta=',num2str(beta)]);
%%
Ngas=size(hist_gas,1);
hand3=figure();hand4=figure();
for h=1:Ngas
[Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_gas{h},pmass_gas(h),snap_gas(h,2),1,1,2);
if Ninfall>5
% alpha=0;beta=Csat(h,1)^-0.4;
par=strip_coeff(Chost(h,2),Csat(h,2));
% flag_redistr=2;alpha=2;beta=1;p=[par(3)/alpha,par(4)/alpha/beta];
% flag_resitr=2;alpha=0;beta=0;p=[par(3)/param_gas(h,1),par(4)/param_gas(h,1)/param_gas(h,2)];
flag_redistr=-2;alpha=1.5;beta=1;p=[alpha,beta];
% flag_redistr=-2;alpha=0;beta=0;p=param_gas(h,1:2);
% flag_redistr=-1;alpha=2;beta=1;p=[alpha,beta];
% flag_redistr=-1;alpha=0;beta=0;p=param_gas(h,1:2);
figure(hand3);
xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,p,3,flag_redistr);
figure(hand4);
xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,p,2,flag_redistr);
end
end 

figure(hand3);
subplot(2,1,1);axis([0,0.3,0,1]);hold off;xlabel('\Deltaln(a)');ylabel('Mgas/Mgas(0)');title('simulation,z_{strp}>1,m/M~[0.01,0.1]')%set(gca,'yscale','log');
subplot(2,1,2);axis([0,0.3,0,1]);hold off;xlabel('\Deltaln(a)');ylabel('Mgas/Mgas(0)');title('Best fit Model: $df/dln(a)=\frac{\alpha\sqrt\Delta{-}\beta f/\sqrt{\Omega_{sat}/\Omega_{host}}(M/m)^{1/3}(v/v_H)}{2\pi}$','interpreter','latex','FontSize',16);%set(gca,'yscale','log');
figure(hand4);
plot([0,1],[0,1],'r');axis([0,1,0,1]);hold off;xlabel('simulation');ylabel('model');title(['gas,\alpha=',num2str(alpha),',\beta=',num2str(beta)]);%set(gca,'yscale','log');

%% gas fitting
% idlist=find(pard>0&pard<0.03&par>0.01&par<0.03&par_dm(:,4)<0.66);
% idlist=find(pard>0&par>0&par_dm(:,4)<0.66);
% idlist2=find(pard>0&par>0&par_dm(:,4)>=0.66);
% h=idlist(1)
% h=idlist2(18)
% close all;
h=2;
Mrate(h,:)
[Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_gas{h},pmass_gas(h),snap_gas(h,2),1,1,2);
virsub

flag_plot=0;
flag_redistr=-1;
% if Mratepar(h,1)>0.02, flag_redistr=4; end 
% alpha=1;
% alpha=0:0.05:3;
% alpha=1;betag=0.8;betad=1.9;
alpha=0.01:0.2:3;
betag=0.1:0.2:2;
% alpha=logspace(0,2,10);
% betag=alpha;
xig=zeros(numel(alpha),numel(betag));
for i=1:numel(alpha)
    for j=1:numel(betag)
xig(i,j)=xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[alpha(i),betag(j)],flag_plot,flag_redistr);
    end
end
[ximin,indmin]=min(abs(xig(:)));
[ig,jg]=ind2sub(size(xig),indmin);
level=[0.001,0.01:0.01:0.05,0.1:0.1:0.5];
figure;imagesc(alpha,betag,sqrt(abs(xig(:,:)')));set(gca,'ydir','normal','clim',[min(abs(xig(:))),0.5]);hold on;
[c,handel]=contour(alpha,betag,sqrt(abs(xig(:,:)')),level,'linewidth',2,'linecolor','k');clabel(c);
% figure;c=contourf(alpha,betag,sqrt(abs(xig(:,:)')),level);clabel(c);%set(gca,'xscale','log');
figure;subplot(2,1,1);plot(alpha,sqrt(xig(:,jg)));subplot(2,1,2);plot(betag,sqrt(xig(ig,:)));
% figure;contourf(alpha,betad,log10(abs(xid'))/2,20);
flag_plot=1;
% clc;
xi2=xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[alpha(ig),betag(jg)],flag_plot,flag_redistr);
xi2=xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[1.5,1],1,-2);
xi2=xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[2,1],1,-1);
par=strip_coeff(Chost(h,2),Csat(h,2));
xi2=xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[par(3)/1.5,par(4)/1.5],flag_plot,2);
if flag_redistr==-1
xi2=xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[par(3)/alpha(ig),par(4)/betag(jg)/alpha(ig)],flag_plot,2);
end
%% whole gas history
[Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_gas{h},pmass_gas(h),1,0,0,0);
figure;
% plot(mgas./mdm);
hold on;
plot(mgas/max(mgas),'b:');plot(mdm/max(mdm),'k-- o');plot(mhost/max(mhost),'r-. .');
plot([snap_gas(h,2)-hist_gas{h}.node(1).Nsnap+1,snap_gas(h,2)-hist_gas{h}.node(1).Nsnap+1],[0.01,1]);
plot([snap_gas(h,3)-hist_gas{h}.node(1).Nsnap+1,snap_gas(h,3)-hist_gas{h}.node(1).Nsnap+1],[0.01,1],'--');
set(gca,'yscale','log');
% legend('\Omega','gas','DM','host','Dstrp','Merge');
legend('gas','DM','host','Dstrp','Merge');
plot(r./Rvir,'b- *');

%% dm fitting
% close all;
% idlist=find(param_dm(:,1)>5);
% idlist=find(dominance>2);
% h=idlist(12);
h=1;
% [Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_dm{h},pmass_dm(h),snap_dm(h,2),1,1,1);
[Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_bin{h},DMpmass,Snappar(h,2),1,1,1);
virsub
% alpha=[1,2];
% beta=[1,1];
flag_plot=0;
flag_redistr=-5;
kt=ones(size(t));
% [x,xi2]=fmincon(@(params) xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,params,0,flag_redistr),[1,0.5],[],[],[],[],[0,0],[2,1]);
[x,xi2]=fmincon(@(params) xi2_DMstrp_history(Hz,Csatpar(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,params,0,flag_redistr),[1,0.5],[],[],[],[],[0,0],[2,1]);
% [x,xi2]=fminsearch(@(params) xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,params,0,flag_redistr),[0.3,3])
% x=[0.5,1.5];
% xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,x,1,fl
% ag_redistr);hold on;title(['e2=',num2str(1-j2D(h,:))]);
xi2_DMstrp_history(Hz,Csatpar(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,x,1,flag_redistr);hold on;title(['e2=',num2str(1-j2D(h,:))]);
% par=strip_coeff(ChostD(h,2),CsatD(h,2));
% xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,par(1:2),1,flag_redistr);hold on;title(['e2=',num2str(1-j2D(h,:))]);
%%
% alpha=0.36;betad=1.5;
% alpha=0.085;betad=0.93;
% alpha=0.1:0.03:3;
alpha=0.1:0.1:2;
% alpha=0.8:0.1:1.2;
% betad=2;
betad=0.3:0.1:1;
% alpha=1;betag=0.8;betad=1.9;
xid=zeros(numel(alpha),numel(betad));
for i=1:numel(alpha)
    for j=1:numel(betad)
xid(i,j)=xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,[alpha(i),betad(j)],0,flag_redistr);
    end
end
level=[0.001,0.01:0.01:0.05,0.1:0.1:0.5];
figure;imagesc(alpha,betad,sqrt(abs(xid(:,:)')));set(gca,'ydir','normal','clim',[min(abs(xid(:))),0.5]);hold on;
[c,handel]=contour(alpha,betad,sqrt(abs(xid(:,:)')),level,'linewidth',2,'linecolor','k');clabel(c);
% figure;
% imagesc(alpha,betad,log(sqrt(abs(xid(:,:)'))));set(gca,'ydir','normal');hold on;
% [c,handel]=contour(alpha,betad,sqrt(abs(xid(:,:)')),level,'linewidth',2,'linecolor','k');clabel(c);
% [c,handel]=contourf(alpha,betad,sqrt(abs(xid(:,:)')),level);clabel(c);%set(gca,'xscale','log');
[ximin,indmin]=min(abs(xid(:)));
[id,jd]=ind2sub(size(xid),indmin);
figure;subplot(2,1,1);plot(alpha,sqrt(xid(:,jd)));subplot(2,1,2);plot(betad,sqrt(xid(id,:)));
flag_plot=1;
xi2=xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,[alpha(id),betad(jd)],flag_plot,flag_redistr);
if flag_redistr==-5
par=strip_coeff(ChostD(h,2),CsatD(h,2));
xi2=xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,[par(1)/alpha(id),par(2)/alpha(id)/betad(jd)],1,0);
end
%% whole dm history
% [Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_dm{h},pmass_dm(h),1,0,0,0);
[Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_bin{h},DMpmass,Snappar(h,2),1,1,1);
figure;
% plot(mgas./mdm);
hold on;
plot(mgas/max(mgas),'b:');plot(mdm/max(mdm),'k-- o');plot(mhost/max(mhost),'r-. .');
plot([snap_dm(h,2)-hist_dm{h}.node(1).Nsnap+1,snap_dm(h,2)-hist_dm{h}.node(1).Nsnap+1],[0.01,1]);
plot([snap_dm(h,3)-hist_dm{h}.node(1).Nsnap+1,snap_dm(h,3)-hist_dm{h}.node(1).Nsnap+1],[0.01,1],'--');
set(gca,'yscale','log');
% legend('\Omega','gas','DM','host','Dstrp','Merge');
legend('gas','DM','host','Dstrp','Merge');
plot(r./Rvir,'b- *');

%%
f=mdm/mdm(1);y=diff(f)./diff(t)*2*pi./sqrt(virialF(1:end-1));x=f(1:end-1).*Rvir(1:end-1)./r(1:end-1);
p=polyfit(x,y,1);
cftool(x,y)
%%
%% DM stat
Nhist=size(hist_bin,1);
hlist=[];
yy=[];xx=[];
Mrange=[0,0.001];
for h=1:Nhist
    [Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_bin{h},Snappar(h,2),1,1,1);
    if Ninfall>20&&Mratepar(h,2)>Mrange(1)&&Mratepar(h,1)<Mrange(2)
    f=mdm/mdm(1);y=diff(f)./diff(t)*2*pi./sqrt(virialF(1:end-1));x=f.*Rvir./r;x=(x(1:end-1)+x(2:end))/2;    
    yy=[yy;y];xx=[xx;x];
    end
end
p=polyfit(xx,yy,1)
cftool(xx,yy)
%% plot concentration dependence
cs=logspace(-1,1,10);
ch=logspace(-1,1,5);
par1=zeros(numel(ch),numel(cs));
par2=par1;par4=par1;
for i=1:numel(ch)
    for j=1:numel(cs)
        coeff=strip_coeff(ch(i),cs(j));
        par1(i,j)=coeff(1);par2(i,j)=coeff(2);par4(i,j)=coeff(4);
    end
end
figure;subplot(1,3,1);semilogx(cs,par1(1,:));xlabel('sat concentration');ylabel('A');
subplot(1,3,2);semilogx(repmat(cs',size(ch)),par2');legend(cellstr(num2str(ch','%2.1f')));xlabel('sat concentration');ylabel('B_d');title('coefficients dependence on sat and host (shown in legend) concentrations');
subplot(1,3,3);semilogx(repmat(cs',size(ch)),par4');legend(cellstr(num2str(ch','%2.1f')));xlabel('sat concentration');ylabel('B_g');


    
%% Model Concentration Dependence
mass=@(x) log(1+x)-x./(1+x);
mgrad=@(x) x.^2./(1+x).^2./mass(x); %dlnM/dlnx, x=r/rs
vc=@(x) sqrt(mass(x)./x);
rhox=@(x) x./(1+x).^2;
rhoc=@(c) c./mass(c); %F=rhoc*rhox; F=rho*r^2/rv^2;
fdmh=@(x) vc(x).*sqrt(3-mgrad(x)); %dependence of tidal radius on host profile
% fgass=@(x) 1./vc(x)./sqrt(rhox(x));
fgass=@(x) 1./sqrt(rhox(x));
fgash=@(x) sqrt(rhox(x));
x=logspace(-2,1,20);
figure;plot(x,1./vc(x).^2);
figure;loglog(x,vc(x).*sqrt(rhox(x)));
c=logspace(-1,1,20);
Adm=zeros(size(c));
Bdms=Adm;
Bdmh=Adm;
Agas=sqrt(5/12)*ones(size(c));
Bgass=Adm;
Bgash=Adm;
for i=1:20
    Adm(i)=quad(vc,0.01*c(i),c(i))/vc(c(i))/c(i)/sqrt(2);
    Bdms(i)=quad(mgrad,0.01*c(i),c(i))/c(i);
    Bdmh(i)=quad(fdmh,0.01*c(i),c(i))/vc(c(i))/c(i);
%     Bgass(i)=Bdms(i)*quad(fgass,0.01*c(i)c(i))*vc(c(i))/sqrt(rhoc(c(i)))/c(i);
    Bgass(i)=Bdms(i)*quad(fgass,0.01*c(i),c(i))/sqrt(rhoc(c(i)))/c(i);
    Bgash(i)=quad(fgash,0.01*c(i),c(i))*sqrt(rhoc(c(i)))/c(i);
end
Bdm=Bdms'*Bdmh/sqrt(2); %colum vec: same host c; line vec: same sat c;
Bgas=Bgass'*Bgash/sqrt(6/5);

betag=0.7;alphag=1.5;
betad=3;alphad=0.2;
% figure(2);

subplot(2,2,1);hold on;plot(c,Adm/betag);
% subplot(2,2,1);hold on;plot(c,Agas/betag);
subplot(2,2,3);hold on;plot(c,Adm/betad);
subplot(2,2,4);hold on;plot(repmat(c',size(c)),Bdm/alphad/betad);
subplot(2,2,2);hold on;plot(repmat(c',size(c)),Bgas/alphag/betag);title('gas');
figure;plot(c,Bdmh);
figure;plot(c,Bgash);

nfw=@(x) log(1+x)-x./(1+x);
mr=@(x) x.^2./(1+x).^2./nfw(x); %dlnM/dlnx, x=r/rs
vc=@(x) sqrt(nfw(x)./x);
fun=@(x) mr(x).*vc(x).^3;
fung=@(x) mr(x).*vc(x).^2;
funb=@(x) mr(x)./vc(x);
funf=@(x) x./(1+x).^2;

x=logspace(-2,1,20);
figure;semilogx(x,mr(x),'- .',x,mr(x).*vc(x).^2,': .');
figure;semilogx(x,vc(x).*sqrt(3-mr(x)));
figure;plot(x,mr(x)./vc(x));
figure;plot(x,mr(x).*vc(x).^2);

c=logspace(-1,1,10);
y=zeros(size(c));
z=y;
yy=y;
u=y;v=y;w=y;
for i=1:10
    y(i)=quad(fun,0.01*c(i),c(i))/vc(c(i)).^3/c(i);
    z(i)=quad(mr,0.01*c(i),c(i))/c(i);
    yy(i)=quad(vc,0.01*c(i),c(i))/vc(c(i))/c(i);
    u(i)=quad(fung,0.01*c(i),c(i))/vc(c(i))^2/c(i);
    v(i)=quad(funf,0.01*c(i),c(i))/c(i);
    w(i)=quad(funb,0.01*c(i),c(i))*vc(c(i))/c(i);
end
figure;semilogx(c,u);
figure;semilogx(c,w./(v.*c./nfw(c)));
figure(2);hold on;semilogx(c,y);figure(4);hold on;semilogx(c,z,'r');hold on;semilogx(c,yy.^3.*z,'k');

p=@(x) (3*(log(1+x)-x./(1+x))-(x./(1+x)).^2).^0.5./(sqrt(x)).*(x.^2./(1+x).^2./(log(1+x)-x./(1+x)));
q=@(c) sqrt(c./(log(1+c)-c./(1+c)));
c=logspace(-1,1,10);
y=zeros(size(c));
for i=1:10
    y(i)=quad(p,0.01*c(i),c(i))/c(i)*q(c(i));
end
figure;semilogx(c,y);

x=logspace(-1,1,10);
figure;semilogx(x,sqrt(vc2(x)./vc2(5)));
x=logspace(-1,0,10);
figure;
for c=[0.1,0.5,1,2,3,5,8]
plot(x,vc2(x*c)./vc2(c).*x);hold on;
end

rho=@(x) x./(1+x).^2;
x=logspace(-2,1,10);
figure;plot(x,rho(x),'- .');
%%  sound speed
alp=2;
nfw=@(x) 1./x./(1+x).^2;
f=@(x) log(1+x)-x./(1+x);
g=@(x) x.*(1+x).^2;
T=@(x) f(x).^1.3./g(x).^(2/3);
TT=@(x) x.^(2/3*(alp-1)).*(1+x).^(-4/3);
K=@(x) TT(x)./nfw(x).^(2/3);
vc2=@(x) (log(1+x)-x./(1+x))./x;
G=43007.1;Mv=10000;Rv=comoving_virial_radius(Mv,1);
coef=G*Mv/Rv/(5/3)/T(1);
coefTT=G*Mv/Rv/(5/3)/TT(1);
coefvc=G*Mv/Rv/(5/3)/vc2(1);
boltzman=1.38e-16;
c=3e5;
proton=938*1.6e-6/c^2;
km=boltzman/proton;
x=logspace(-2,1,20);
figure;loglog(f(x),K(x),'.');xlabel('M(>r)');ylabel('entropy (r)');
hold on;p=polyfit(log(f(x)),log(K(x)),1);loglog(f(x),exp(p(2))*f(x).^p(1));title(['K=M^{',num2str(p(1),'%.1f'),'}']);%set alp=3 to get K=M^1.3
figure;loglog(x*Rv/4,coef*T(x)/km,'r',x*Rv/4,coefTT*TT(x)/km,'g')
y=x;
yy=x;
for i=1:numel(x)
    y(i)=quad(@(x) 1./sqrt(T(x)),0.001,x(i));
    yy(i)=quad(@(x) 1./sqrt(TT(x)),0.001,x(i));
end
y=y./sqrt(coef);yy=yy/sqrt(coefTT);
figure;loglog(x,y,'r-',x,yy,'g-',x,2*pi*x./sqrt(T(x)*coef),'r:',x,2*pi*x./sqrt(TT(x)*coefTT),'g:',x,2*pi*x./sqrt(vc2(x)*coefvc),'k');
figure;loglog(x,y./x,'r-',x,1./sqrt(T(x)*coef),'r:',x,yy./x,'g-',x,1./sqrt(TT(x)*coefTT),'g:',x,1./sqrt(vc2(x)*coefvc),'k');
figure;semilogx(x,yy./x.*sqrt(vc2(x))./(yy(end)./x(end).*sqrt(vc2(x(end)))),'r')
% figure;loglog(x,nfw(x).*vc2(x));
%%
Ndm=size(hist_dm,1);
betad=zeros(Ndm,1);alphad=zeros(Ndm,1);gammad=zeros(Ndm,1);
for h=1:Ndm
    par=strip_coeff(ChostD(h,2),CsatD(h,2));
    betad(h)=par(1)/param_dm(h,1);
    gammad(h)=par(2)/param_dm(h,2);
end
alphad=gammad./betad;

figure('name','parameter distributions');
subplot(3,1,1);hist(betad(betad<5),20);xlabel('\beta_d');subplot(3,1,2);hist(gammad,20);xlabel('\alpha_d*\beta_d');
subplot(3,1,3);hist(alphad(betad<5),20);xlabel('\alpha_d');

figure('name','Mass ratio dependence');
subplot(3,1,1);semilogx(MrateD(:,2),betad,'.');ylabel('\alpha_{dm}');xlabel('m/M');set(gca,'ylim',[0,5]);
subplot(3,1,2);semilogx(MrateD(:,2),gammad,'.');ylabel('\beta_{dm}');xlabel('m/M');
subplot(3,1,3);semilogx(MrateD(:,2),alphad,'.');ylabel('\beta_{dm}');xlabel('m/M');

figure('name','Host Mass dependence');
subplot(3,1,1);semilogx(MhostD(:,2),betad,'.');ylabel('\alpha_{dm}');xlabel('M');set(gca,'ylim',[0,5]);
subplot(3,1,2);semilogx(MhostD(:,2),gammad,'.');ylabel('\beta_{dm}');xlabel('M');
subplot(3,1,3);semilogx(MhostD(:,2),alphad,'.');ylabel('\beta_{dm}');xlabel('M');

figure('name','Sat Mass dependence');
subplot(3,1,1);semilogx(MsatD(:,2),betad,'.');ylabel('\alpha_{dm}');xlabel('m');set(gca,'ylim',[0,5]);
subplot(3,1,2);semilogx(MsatD(:,2),gammad,'.');ylabel('\beta_{dm}');xlabel('m');
subplot(3,1,3);semilogx(MsatD(:,2),alphad,'.');ylabel('\beta_{dm}');xlabel('m');

figure('name','Host concentration dependence');
subplot(3,1,1);plot(ChostD(:,2),betad,'.');ylabel('\alpha_{dm}');xlabel('Chost');set(gca,'ylim',[0,5]);
subplot(3,1,2);plot(ChostD(:,2),gammad,'.');ylabel('\beta_{dm}');xlabel('Chost');
subplot(3,1,3);plot(ChostD(:,2),alphad,'.');ylabel('\beta_{dm}');xlabel('Chost');

figure('name','Sat concentration dependence');
subplot(3,1,1);plot(CsatD(:,1),betad,'.');ylabel('\alpha_{dm}');xlabel('Csat');set(gca,'ylim',[0,5]);
subplot(3,1,2);plot(CsatD(:,1),gammad,'.');ylabel('\beta_{dm}');xlabel('Csat');
x=CsatD(:,1);y=gammad;x=x(y>0);y=y(y>0);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;plot(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);
subplot(3,1,3);plot(CsatD(:,1),alphad,'.');ylabel('\beta_{dm}');xlabel('Csat');

figure('name','orbit shape dependence');
subplot(3,1,1);semilogx(e2D(:,2),betad,'.');ylabel('\alpha_{dm}');xlabel('e2');set(gca,'ylim',[0,5]);
subplot(3,1,2);semilogx(e2D(:,2),gammad,'.');ylabel('\beta_{dm}');xlabel('e2');
subplot(3,1,3);semilogx(e2D(:,2),alphad,'.');ylabel('\beta_{dm}');xlabel('e2');

figure('name','orbital energy dependence');
subplot(3,1,1);semilogx(KD(:,2),betad,'.');ylabel('\alpha_{dm}');xlabel('K');set(gca,'ylim',[0,5]);
subplot(3,1,2);semilogx(KD(:,2),gammad,'.');ylabel('\beta_{dm}');xlabel('K');
subplot(3,1,3);semilogx(KD(:,2),alphad,'.');ylabel('\beta_{dm}');xlabel('K');    
%%
Ndm=size(hist_dm,1);par=zeros(Ndm,4);
for h=1:Ndm
    par(h,:)=strip_coeff(ChostD(h,2),CsatD(h,2));
end

alpha=0.5:0.1:1.5;
beta=0.1:0.1:1;
% alpha=1;beta=0.8;
xid=zeros(numel(alpha),numel(beta));
for i=1:numel(alpha)
    for j=1:numel(beta)
        for h=1:Ndm
        [Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_dm{h},pmass_dm(h),snap_dm(h,2),1,1,1);
%         xi2=xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,[par(h,1)/alpha(i),par(h,2)/alpha(i)/beta(j)],0,0);
        xi2=xi2_DMstrp_history(Hz,CsatD(h,2),r,vt2,kt,t,mdm,chost,Rvir,virialF,[alpha(i),beta(j)],0,-5);
        xid(i,j)=xid(i,j)+xi2;
        end
    end
end
xid=sqrt(xid/Ndm);
level=[0.001,0.01:0.01:0.05,0.1:0.05:0.5];
figure;
imagesc(alpha,beta,(abs(xid(:,:)')));set(gca,'ydir','normal','clim',[min(abs(xid(:))),0.5]);hold on;
[c,handel]=contour(alpha,beta,abs(xid(:,:)'),level,'linewidth',2,'linecolor','k');clabel(c);
% [c,handel]=contourf(alpha,beta,sqrt(abs(xid(:,:)')),level);clabel(c);%set(gca,'xscale','log');
[ximin,indmin]=min(abs(xid(:)));
[id,jd]=ind2sub(size(xid),indmin);
level=[0.1:0.02:0.5];
[c,handel]=contour(alpha,beta,abs(xid(:,:)'),level);clabel(c);xlabel('\alpha_{dm}');ylabel('\beta_{dm}');title('corrected isothermal model, overall fit goodness');


Ngas=size(hist_gas,1);par=zeros(Ngas,4);
for h=1:Ngas
    par(h,:)=strip_coeff(Chost(h,2),Csat(h,2));
end

alphag=1.5:0.1:2.5;
betag=0.5:0.1:1.5;
% alphag=2;beta=1;
xig=zeros(numel(alphag),numel(betag));
for i=1:numel(alphag)
    for j=1:numel(betag)
        for h=1:Ngas
        [Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_gas{h},pmass_gas(h),snap_gas(h,2),1,1,2);
%         xi2=xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[par(h,3)/alphag(i),par(h,4)/alphag(i)/betag(j)],0,2);
        xi2=xi2_GASstrp_history(pmass_gas(h),Csat(h,2),r,v,t,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,[alphag(i),betag(j)],0,-1);
        xig(i,j)=xig(i,j)+xi2;
        end
    end
end
xig=sqrt(xig/Ngas);
level=[0.001,0.01:0.01:0.05,0.1:0.05:0.5];
figure;
imagesc(alphag,betag,(abs(xig(:,:)')));set(gca,'ydir','normal','clim',[min(abs(xid(:))),0.5]);hold on;
[c,handel]=contour(alphag,betag,abs(xig(:,:)'),level,'linewidth',2,'linecolor','k');clabel(c);
% [c,handel]=contourf(alpha,beta,sqrt(abs(xid(:,:)')),level);clabel(c);%set(gca,'xscale','log');
[ximin,indmin]=min(abs(xig(:)));
[ig,jg]=ind2sub(size(xig),indmin);
level=[0.1:0.02:0.5];
[c,handel]=contour(alphag,betag,abs(xig(:,:)'),level);clabel(c);xlabel('\alpha_{gas}');ylabel('\beta_{gas}');title('NFW stripping model, overall fit goodness');

%%
dominance=zeros(Ndm,1);
for h=1:Ndm
    [Ninfall,t,r,v,vt2,kt,mdm,mgas,mhost,fhost,chost,Umsat,Hz,Rvir,virialF,virsub]=get_strp_history(hist_dm{h},pmass_dm(h),snap_dm(h,2),1,1,1);
    dominance(h)=max(virsub);
end
hist(dominance);
    