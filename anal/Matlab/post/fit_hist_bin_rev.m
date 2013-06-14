clear;
global G a pmass HUBBLE0 Omega0 OmegaLambda
G=43007.1;HUBBLE0=0.1;Omega0=0.3;OmegaLambda=0.7;

runnum=6702;
pmass=0.008848;

scaleF_file='/home/kam/Documents/research/Galaxy/code/BoundTracing/v7.9/anal/Matlab/outputs_zoom.txt';
a=load(scaleF_file);
history_file=['/mnt/A4700/data/',num2str(runnum),'/subcat/anal/history_000_005.dat'];
history_par_file=['/mnt/A4700/data/',num2str(runnum),'/subcat/anal/historypar_000_005.dat'];

disp('loading files...')
hist_bin=load_hist_bin(history_file);
[Snappar,Mratepar,Mhostpar,Rhostpar,Kpar,j2par,Kerrpar,j2errpar,Chostpar,Csatpar]=load_hist_par(history_par_file);
%% GAS stat
glist=[];
par_gas=zeros(Nhist,4);%alpha,beta,xi2,Kself
Dstrp=zeros(Nhist,1);
virsubs=zeros(Nhist,1);
flag_plot=0;
% alpha=1;
alpha=0.1:0.05:5;
betag=0.1:0.05:3;
% betag=1;
xig=zeros(numel(alpha),numel(betag));
for h=1:Nhist
    [Ninfall,t,r,v,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,virsubs(h),Kself]=get_strp_history(hist_bin{h},Snappar(h,2),1,1,2);
    Dstrp(h)=r(1)/Rvir(1);
    if Ninfall>5&mgas(1)>500
        h
        glist=[glist,h];
        for i=1:numel(alpha)
            for j=1:numel(betag)
                    xig(i,j)=xi2_GASstrp_history(r,v,t,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,[alpha(i),betag(j)],flag_plot,2);
            end
        end
        [ximin,indmin]=min(abs(xig(:)));
        if isreal(xig(indmin))&&(~isnan(xig(indmin)))
            [ig,jg]=ind2sub(size(xig),indmin);
            par_gas(h,:)=[alpha(ig),betag(jg),sqrt(ximin),Kself(1)];
        end
    end
end
par=par_gas(:,3);
% figure;subplot(2,1,1);hist((par(logical(par))),20);xlabel('<\Delta f_{gas}>');ylabel('Counts');
% subplot(2,1,2);hist((par_dm(logical(par),3)),20);xlabel('<\Delta f_{dm}>');ylabel('Counts');
x=linspace(0,0.4,40);
figure;h1=subplot(2,1,1);h2=subplot(2,1,2);
filter=logical(par>0&par_gas(:,4)<0.7);
    y=histc(par_gas(filter,3),x);
    subplot(h1);plot(x,y,'k- .');hold on;
    subplot(h2);plot(x,cumsum(y),'k- .');hold on;
%% GAS corr
infall_type=2;
par_id=1;
% filter=logical(par);
filter=logical(par>0&par<0.05&par_gas(:,4)<0.66);
% figure;plot(virsubs(filterd),par_dm(filterd,4),'.');xlabel('sub prominence');ylabel('dynamical compactness');
% figure;subplot(2,1,1);plot(virsubs(filter),par(filter),'.');subplot(2,1,2);plot(virsubs(filterd),pard(filterd),'.');xlabel('sub dominence');ylabel('goodness');
figure;plot(par_gas(filter,4),par(filter),'.');xlabel('compactness');ylabel('goodness');
% figure;subplot(2,1,1);plot(par_gas(filter,4),par_gas(filter,2),'.');subplot(2,1,2);plot(par_dm(filterd,4),par_dm(filterd,2),'.');xlabel('compactness');ylabel('beta');
% figure;subplot(2,1,1);plot(par_gas(filter,1),par_gas(filter,2),'.');xlabel('\alpha');ylabel('\beta');title('gas');
% subplot(2,1,2);plot(par_dm(filter,1),par_dm(filter,2),'.');xlabel('\alpha');ylabel('\beta');title('dm');
figure;hist(par_gas(filter,par_id),20);title('gas');

% orbital parameters=======
% figure;subplot(2,2,1);semilogy(par_gas(filter,1),Kpar(filter,infall_type),'.');xlabel('\alpha');ylabel('K');
% subplot(2,2,3);semilogy(par_gas(filter,2),Kpar(filter,infall_type),'.');xlabel('\beta');ylabel('K');

% subplot(2,2,2);semilogy(par_gas(filter,1),1-j2par(filter,infall_type),'.');xlabel('\alpha');ylabel('e2');
% subplot(2,2,4);semilogy(par_gas(filter,2),1-j2par(filter,infall_type),'.');xlabel('\beta');ylabel('e2');

figure;subplot(2,1,1);loglog(Mratepar(filter,infall_type),par_gas(filter,par_id),'.');ylabel('\betag');xlabel('m/M');
x=Mratepar(filter,infall_type);y=par_gas(filter,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);
ylabel('\betad');xlabel('m/M');

subplot(2,1,2);loglog(Mhostpar(filter,infall_type).*Mratepar(filter,infall_type),par_gas(filter,par_id),'.');xlabel('Msat');ylabel('\betag');
x=Mhostpar(filter,infall_type).*Mratepar(filter,infall_type);y=par_gas(filter,par_id);p=polyfit(log(x),log(y),1); 
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

figure;subplot(2,1,1);loglog(Mhostpar(filter,infall_type),par_gas(filter,par_id),'.');ylabel('\betag');xlabel('Mhost');
x=Mhostpar(filter,infall_type);y=par_gas(filter,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

subplot(2,1,2);loglog(Rhostpar(filter,infall_type),par_gas(filter,par_id),'.');ylabel('\betag');xlabel('Rhost');
x=Rhostpar(filter,infall_type);y=par_gas(filter,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

figure;subplot(2,1,1);loglog(Chostpar(filter,infall_type),par_dm(filter,par_id),'.');ylabel('\betad');xlabel('Chost');
x=Chostpar(filter,infall_type);y=par_dm(filter,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

filterc=logical(Csatpar(:,infall_type)<inf)&filter;
subplot(2,1,2);loglog(Csatpar(filterc,infall_type),par_dm(filterc,par_id),'.');ylabel('\betad');xlabel('Csat');
x=Csatpar(filterc,infall_type);y=par_dm(filterc,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

figure;plot3(log(Mhostpar(filter,infall_type).*Mratepar(filter,infall_type)),log(Mhostpar(filter,infall_type)),log(par_gas(filter,par_id)),'.');
x1=log(Mhostpar(filter,infall_type).*Mratepar(filter,infall_type));x2=log(Mhostpar(filter,infall_type));X=[x1,x2,ones(size(x1))];y=log(par_gas(filter,par_id));
[p,se_p,mse] = lscov(X,y);
x1grid=linspace(min(x1),max(x1),5)';x2grid=linspace(min(x2),max(x2),5)';X=[x1grid,x2grid,ones(5,1)];
[x11,x22]=meshgrid(x1grid,x2grid);yy=p(1)*x11+p(2)*x22+p(3);
hold on;mesh(x11,x22,yy);
hold on;plot3(x1grid,x2grid,X*p,'--');
% hold on;plot3(x1grid,x2grid,X*[0;-0.1;log(3)],'r-');yy=-0.1*x22+log(3);mesh(x11,x22,yy);
title(['\beta=',num2str(exp(p(3)),'%.1f'),'m^{',num2str(p(1),'%.2f'),'}','M^{',num2str(p(2),'%.2f'),'}',sprintf('\n'),'stderr:',num2str(exp(se_p(3))),',',num2str(se_p(1)),',',num2str(se_p(2))]);
x=logspace(-3,-0,20);
figure;
    y=histc(par_gas(filter,3),x);
    subplot(2,1,1);plot(x,y,'- .');set(gca,'xscale','log');
    subplot(2,1,2);plot(x,cumsum(y),'- .');set(gca,'xscale','log');
%  y=histc(goodness{3},x);figure(1);hold on;
%     subplot(2,1,1);plot(x,y,[colors(3),': .']);hold on;
%     subplot(2,1,2);plot(x,cumsum(y),[colors(3),': .']);hold on;

%% DM stat
Nhist=size(hist_bin,1);
hlist=[];
par_dm=zeros(Nhist,4);%alpha,beta,xi2,Kself
Dstrp=zeros(Nhist,1);
virsubs=zeros(Nhist,1);
flag_plot=0;
flag_redistr=0;
% alpha=1;
alpha=0.1:0.05:3;
betad=0.5:0.05:4;
% betad=1.7;
xid=zeros(numel(alpha),numel(betad));
for h=1:Nhist
    [Ninfall,t,r,v,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,virsubs(h),Kself]=get_strp_history(hist_bin{h},Snappar(h,2),1,1,1);
    Dstrp(h)=r(1)/Rvir(1);
    if Ninfall>15
        h
        hlist=[hlist,h];
        for i=1:numel(alpha)
            for j=1:numel(betad)
                xid(i,j)=xi2_DMstrp_history(r,t,mdm,Rvir,virialF,[alpha(i),betad(j)],flag_plot,flag_redistr);
            end
        end
        [ximin,indmin]=min(abs(xid(:)));
        if isreal(xid(indmin))&&(~isnan(xid(indmin)))
            [id,jd]=ind2sub(size(xid),indmin);
            par_dm(h,:)=[alpha(id),betad(jd),sqrt(ximin),Kself(1)];
        end
    end
end
pard=par_dm(:,3);
% figure;subplot(2,1,1);hist((par(logical(par))),20);xlabel('<\Delta f_{gas}>');ylabel('Counts');
% subplot(2,1,2);hist((par_dm(logical(par),3)),20);xlabel('<\Delta f_{dm}>');ylabel('Counts');
x=linspace(0,0.4,40);
figure;h1=subplot(2,1,1);h2=subplot(2,1,2);
filter=logical(pard>0&pard<0.05&par_dm(:,4)<0.7);
y=histc(par_dm(filter,3),x);
subplot(h1);plot(x,y,'k- .');hold on;
subplot(h2);plot(x,cumsum(y),'k- .');
figure;subplot(2,1,1);hist(par_dm(filter,1),numel(alpha));subplot(2,1,2);hist(par_dm(filter,2),numel(betad));
%% DM corr
infall_type=2;
par_id=1;
% filterd=logical(pard);
filterd=logical(par_dm(:,par_id)<40&pard>0&pard<0.03&par_dm(:,4)<0.9);
% figure;plot(virsubs(filterd),par_dm(filterd,4),'.');xlabel('sub prominence');ylabel('dynamical compactness');
% figure;subplot(2,1,1);plot(virsubs(filter),par(filter),'.');subplot(2,1,2);plot(virsubs(filterd),pard(filterd),'.');xlabel('sub dominence');ylabel('goodness');
figure;plot(par_dm(filterd,4),pard(filterd),'.');xlabel('compactness');ylabel('goodness');
% figure;subplot(2,1,1);plot(par_gas(filter,4),par_gas(filter,2),'.');subplot(2,1,2);plot(par_dm(filterd,4),par_dm(filterd,2),'.');xlabel('compactness');ylabel('beta');
figure;plot(par_dm(pard>0.0001&pard<0.02,par_id),Dstrp(pard>0.0001&pard<0.02),'.');xlabel('\beta_{DM}');ylabel('Dstrp');
% figure;subplot(2,1,1);plot(par_gas(filter,1),par_gas(filter,2),'.');xlabel('\alpha');ylabel('\beta');title('gas');
% subplot(2,1,2);plot(par_dm(filter,1),par_dm(filter,2),'.');xlabel('\alpha');ylabel('\beta');title('dm');
figure;hist(par_dm(filter,par_id),20);title('dm');

% orbital parameters=======
% figure;subplot(2,2,1);semilogy(par_gas(filter,1),Kpar(filter,infall_type),'.');xlabel('\alpha');ylabel('K');
% subplot(2,2,3);semilogy(par_gas(filter,2),Kpar(filter,infall_type),'.');xlabel('\beta');ylabel('K');

% subplot(2,2,2);semilogy(par_gas(filter,1),1-j2par(filter,infall_type),'.');xlabel('\alpha');ylabel('e2');
% subplot(2,2,4);semilogy(par_gas(filter,2),1-j2par(filter,infall_type),'.');xlabel('\beta');ylabel('e2');

figure;subplot(2,1,1);loglog(Mratepar(filterd,infall_type),par_dm(filterd,par_id),'.');ylabel('\betad');xlabel('m/M');
x=Mratepar(filterd,infall_type);y=par_dm(filterd,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

subplot(2,1,2);loglog(Mhostpar(filterd,infall_type).*Mratepar(filterd,infall_type),par_dm(filterd,par_id),'.');xlabel('Msat');ylabel('\betad');
x=Mhostpar(filterd,infall_type).*Mratepar(filterd,infall_type);y=par_dm(filterd,par_id);p=polyfit(log(x),log(y),1); 
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

figure;subplot(2,1,1);loglog(Mhostpar(filterd,infall_type),par_dm(filterd,par_id),'.');ylabel('\betad');xlabel('Mhost');
x=Mhostpar(filterd,infall_type);y=par_dm(filterd,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

subplot(2,1,2);loglog(Rhostpar(filterd,infall_type),par_dm(filterd,par_id),'.');ylabel('\betad');xlabel('Rhost');
x=Rhostpar(filterd,infall_type);y=par_dm(filterd,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

figure;subplot(2,1,1);loglog(Chostpar(filterd,infall_type),par_dm(filterd,par_id),'.');ylabel('\betad');xlabel('Chost');
x=Chostpar(filterd,infall_type);y=par_dm(filterd,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);

filterc=logical(Csatpar(:,infall_type)<inf)&filterd;
subplot(2,1,2);loglog(Csatpar(filterc,infall_type),par_dm(filterc,par_id),'.');ylabel('\betad');xlabel('Csat');
x=Csatpar(filterc,infall_type);y=par_dm(filterc,par_id);p=polyfit(log(x),log(y),1);
x=logspace(log10(min(x)),log10(max(x)),5);hold on;loglog(x,exp(p(2))*x.^p(1),'--');hold off;
legend('data',['\beta=',num2str(exp(p(2)),'%.1f'),'x^{',num2str(p(1),'%.2f'),'}']);



figure;plot3(log(Mhostpar(filterd,infall_type).*Mratepar(filterd,infall_type)),log(Mhostpar(filterd,infall_type)),log(par_dm(filterd,par_id)),'.');
x1=log(Mhostpar(filterd,infall_type).*Mratepar(filterd,infall_type));x2=log(Mhostpar(filterd,infall_type));X=[x1,x2,ones(size(x1))];y=log(par_dm(filterd,par_id));
[p,se_p,mse] = lscov(X,y);
x1grid=linspace(min(x1),max(x1),5)';x2grid=linspace(min(x2),max(x2),5)';X=[x1grid,x2grid,ones(5,1)];
[x11,x22]=meshgrid(x1grid,x2grid);yy=p(1)*x11+p(2)*x22+p(3);
hold on;mesh(x11,x22,yy);
hold on;plot3(x1grid,x2grid,X*p,'--');
% hold on;plot3(x1grid,x2grid,X*[0;-0.1;log(3)],'r-');yy=-0.1*x22+log(3);mesh(x11,x22,yy);
% hold on;plot3(x1grid,x2grid,X*[-0.1;0;log(2.4)],'r-');yy=-0.1*x11+log(2.4);mesh(x11,x22,yy);
title(['\beta=',num2str(exp(p(3)),'%.1f'),'m^{',num2str(p(1),'%.2f'),'}','M^{',num2str(p(2),'%.2f'),'}',sprintf('\n'),'stderr:',num2str(exp(se_p(3))),',',num2str(se_p(1)),',',num2str(se_p(2))]);
x=logspace(-3,-0,20);
figure;
    y=histc(par_dm(filterd,3),x);
    subplot(2,1,1);plot(x,y,'- .');set(gca,'xscale','log');
    subplot(2,1,2);plot(x,cumsum(y),'- .');set(gca,'xscale','log');
%  y=histc(goodness{3},x);figure(1);hold on;
%     subplot(2,1,1);plot(x,y,[colors(3),': .']);hold on;
%     subplot(2,1,2);plot(x,cumsum(y),[colors(3),': .']);hold on;
%%
for h=1:Nhist
[Ninfall,t,r,v,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,virsub,Kself]=get_strp_history(hist_bin{h},Snappar(h,2),1,1);
if Ninfall>5&&Kself(1)<0.7
df_dt_gas=zeros(Ninfall,1);
ft_gas=ones(Ninfall,1);
df_dt_dm=zeros(Ninfall,1);
ft_dm=ones(Ninfall,1);
for i=1:Ninfall
    if i>1
        ft_gas(i)=ft_gas(i-1)+df_dt_gas(i-1)*0.0116;
        ft_dm(i)=ft_dm(i-1)+df_dt_dm(i-1)*0.0116;
    end
	mrate=mdm(i)/mhost(i);
    grate=mgas(i)/mdm(i)/fhost(i);
%     df_dt_gas(i)=(sqrt(virialF(i))*par_gas(h,1)-par_gas(h,2)*(ft_gas(i))/sqrt(grate)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
%     df_dt_dm(i)=sqrt(virialF(i))*(par_dm(h,1)-par_dm(h,2)*(ft_dm(i)*Rvir(i)/r(i))^(2/3))/2/pi;  %tidal stripping
    df_dt_gas(i)=(sqrt(virialF(i))*0.8-0.8*(ft_gas(i))/sqrt(grate)*mrate^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
% df_dt_gas(i)=(sqrt(virialF(i))-4*(ft_gas(i))/sqrt(grate)*mhost(i)^(0.2)*mdm(i)^(-0.3)*(v(i)/r(i)/Hz(i)))/2/pi;
% df_dt_gas(i)=(sqrt(virialF(i))-5*(ft_gas(i))/sqrt(grate)*mhost(i)^(0.2)*mdm(i)^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
%     df_dt_gas(i)=(sqrt(virialF(i))*2.4*mdm(i)^(-0.1)-5*(ft_gas(i))/sqrt(grate)*mhost(i)^(0.2)*mdm(i)^(-1./3)*(v(i)/r(i)/Hz(i)))/2/pi;
    df_dt_dm(i)=sqrt(virialF(i))*(1-1.7*(ft_dm(i)*Rvir(i)/r(i))^(3/3))/2/pi;  %tidal stripping
end
frac_gas=mgas/mgas(1);
frac_dm=mdm/mdm(1);
         figure(1);
            subplot(2,1,1); plot((t-t(1)),frac_gas,'- .','userdata',h);hold on;
            subplot(2,1,2); plot((t-t(1)),ft_gas,'r-- o','userdata',h);hold on;
          figure(2);
            subplot(2,1,1); plot((t-t(1)),frac_dm,'- .','userdata',h);hold on;
            subplot(2,1,2); plot((t-t(1)),ft_dm,'r-- o','userdata',h);hold on;
          figure(3);
            subplot(2,1,1); plot(frac_gas,ft_gas,'.','userdata',h);hold on;
            subplot(2,1,2); plot(frac_dm,ft_dm,'.','userdata',h);hold on;  
end
end 

figure(1);
subplot(2,1,1);axis([0,0.3,0,1]);hold off;xlabel('\Deltaln(a)');ylabel('Mgas/Mgas(0)');title('simulation,z_{strp}>1,m/M~[0.01,0.1]')%set(gca,'yscale','log');
subplot(2,1,2);axis([0,0.3,0,1]);hold off;xlabel('\Deltaln(a)');ylabel('Mgas/Mgas(0)');title('Best fit Model: $df/dln(a)=\frac{\alpha\sqrt\Delta{-}\beta f/\sqrt{\Omega_{sat}/\Omega_{host}}(M/m)^{1/3}(v/v_H)}{2\pi}$','interpreter','latex','FontSize',16);%set(gca,'yscale','log');
figure(2);
subplot(2,1,1);axis([0,0.3,0,1]);hold off;xlabel('\Deltaln(a)');ylabel('Mdm/Mdm(0)');title('simulation,z_{infall}>1,m/M~[0.01,0.1]')%set(gca,'yscale','log');
subplot(2,1,2);axis([0,0.3,0,1]);hold off;xlabel('\Deltaln(a)');ylabel('Mdm/Mdm(0)');title('Best fit Model: $df/dln(a)=\frac{\sqrt\Delta}{2\pi}[\alpha{-}\beta(f\frac{Rvir}{R})^{2/3}]$','interpreter','latex','FontSize',16);%set(gca,'yscale','log');
figure(3);
subplot(2,1,1);plot([0,1],[0,1]);axis([0,1,0,1]);hold off;xlabel('simulation');ylabel('model');title('gas');%set(gca,'yscale','log');
subplot(2,1,2);plot([0,1],[0,1]);axis([0,1,0,1]);hold off;xlabel('simulation');ylabel('model');title('dm');
%% gas fitting
% idlist=find(pard>0&pard<0.03&par>0.01&par<0.03&par_dm(:,4)<0.66);
% idlist=find(pard>0&par>0&par_dm(:,4)<0.66);
% idlist2=find(pard>0&par>0&par_dm(:,4)>=0.66);
% h=idlist(26)
% h=idlist2(18)
close all;
 h=1;
Mratepar(h,:)
[Ninfall,t,r,v,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,virsub,Kself]=get_strp_history(hist_bin{h},Snappar(h,2),1,1,2);
virsub,'Kself=',Kself(1)

flag_plot=0;
flag_redistr=2;
% if Mratepar(h,1)>0.02, flag_redistr=4; end 
% alpha=1;
% alpha=0:0.05:3;
alpha=0.8:0.1:1.2;
% alpha=1;betag=0.8;betad=1.9;
betag=0.1:0.05:2;
xig=zeros(numel(alpha),numel(betag));
for i=1:numel(alpha)
    for j=1:numel(betag)
xig(i,j)=xi2_GASstrp_history(r,v,t,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,[alpha(i),betag(j)],flag_plot,flag_redistr);
    end
end
[ximin,indmin]=min(abs(xig(:)));
[ig,jg]=ind2sub(size(xig),indmin);
% figure;[c,h]=contourf(alpha,betag,log10(abs(xig(:,:)'))/2,20);
level=[0.001,0.01:0.01:0.05,0.1:0.1:0.5];
figure;c=contourf(alpha,betag,sqrt(abs(xig(:,:)')),level);clabel(c);
% figure;contourf(alpha,betad,log10(abs(xid'))/2,20);
flag_plot=1;
xi2=xi2_GASstrp_history(r,v,t,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,[alpha(ig),betag(jg)],flag_plot,flag_redistr);

%% dm fitting
close all;
h=1;
[Ninfall,t,r,v,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,virsub,Kself]=get_strp_history(hist_bin{h},Snappar(h,2),1,1,1);
% alpha=[1,2];
% beta=[1,1];
flag_plot=0;
flag_redistr=0;
% alpha=0.36;betad=1.5;
% alpha=0.085;betad=0.93;
alpha=0.1:0.03:3;
% alpha=0.8:0.1:1.2;
betad=0.5:0.03:3;
% alpha=1;betag=0.8;betad=1.9;
xid=zeros(numel(alpha),numel(betad));
for i=1:numel(alpha)
    for j=1:numel(betad)
xid(i,j)=xi2_DMstrp_history(r,t,mdm,Rvir,virialF,[alpha(i),betad(j)],flag_plot,flag_redistr);
    end
end
% level=[0.001,0.01:0.01:0.05,0.1:0.1:0.5];
% figure;
% imagesc(alpha,betad,log(sqrt(abs(xid(:,:)'))));set(gca,'ydir','normal');hold on;
% [c,handel]=contour(alpha,betad,sqrt(abs(xid(:,:)')),level,'linewidth',2,'linecolor','k');clabel(c);
[ximin,indmin]=min(abs(xid(:)));
[id,jd]=ind2sub(size(xid),indmin);
% adm=[adm;alpha(id),mdm(1)];
flag_plot=1;
xi2=xi2_DMstrp_history(r,t,mdm,Rvir,virialF,[alpha(id),betad(jd)],flag_plot,flag_redistr);
%% whole history
[Ninfall,t,r,v,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,virsub,Kself]=get_strp_history(hist_bin{h},1,0,0,0);
figure;
% plot(mgas./mdm);
hold on;
plot(mgas/max(mgas),'b:');plot(mdm/max(mdm),'k-- o');plot(mhost/max(mhost),'r-. .');
plot([Snappar(h,2)-hist_bin{h}.node(1).Nsnap+1,Snappar(h,2)-hist_bin{h}.node(1).Nsnap+1],[0.01,1]);
plot([Snappar(h,3)-hist_bin{h}.node(1).Nsnap+1,Snappar(h,3)-hist_bin{h}.node(1).Nsnap+1],[0.01,1],'--');
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
    [Ninfall,t,r,v,mdm,mgas,mhost,fhost,Hz,Rvir,virialF,virsubs(h),Kself]=get_strp_history(hist_bin{h},Snappar(h,2),1,1,1);
    if Ninfall>20&&Mratepar(h,2)>Mrange(1)&&Mratepar(h,1)<Mrange(2)
    f=mdm/mdm(1);y=diff(f)./diff(t)*2*pi./sqrt(virialF(1:end-1));x=f.*Rvir./r;x=(x(1:end-1)+x(2:end))/2;    
    yy=[yy;y];xx=[xx;x];
    end
end
p=polyfit(xx,yy,1)
cftool(xx,yy)

