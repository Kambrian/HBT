%%
clear;
sim='/home/kam/Documents/research/Galaxy/code/BoundTracing/v6/sim6702/anal/data';
load /home/kam/Documents/research/Galaxy/code/gas/outputs_zoom.txt
uiimport([sim,'/steller/SnapInfall'])
load([sim,'/prof/LxTidal/CoM/subprof_099'])
load([sim,'/massfun/v60/group_offset_099'])
load([sim,'/massfun/v60/sub_mass_099'])
load([sim,'/prof/subprof_dm_099.3vir'])
% load /mnt/A4700/data/6702/subcat/anal/LxVir/CoM/subprof_099.2
%%
nbin=60;
fact_relax=1.5;
L=subprof_099(:,1:nbin);
gm=group_offset_099(:,1);
Ro=subprof_099(:,nbin+1:2*nbin);
T=subprof_099(:,2*nbin+1:3*nbin);
r=subprof_099(:,3*nbin+1:4*nbin);
n=subprof_099(:,4*nbin+1:5*nbin);
s=subprof_099(:,5*nbin+1:5*nbin+2);
Lt=L.*n;
x=0:nbin;
V=diff(x.^3);
V=repmat(V,size(L,1),1);
Lm=Lt./V;
id=subprof_099(:,end);
rcen=subprof_dm_099(id+1,722);
rvir=subprof_dm_099(gm+1,721);
rratio=[];
for j=1:size(gm,1)
    if group_offset_099(j,2)>1
        for k=(gm(j)+2):(gm(j)+group_offset_099(j,2))
            rratio=[rratio;subprof_dm_099(k,722)/subprof_dm_099(gm(j)+1,721)];
        end
    end
end
clear subprof_dm_099
stellermass;
%%
figure;
subplot(3,1,1);loglog(s(:,2),Ksat(id+1),'.');xlabel('Mdm');ylabel('Lk');set(gca,'xminortick','on','yminortick','on');title('steller vs. DM');
subplot(3,1,2);loglog(s(:,2),s(:,1),'.');xlabel('Mdm');ylabel('Rtidal');set(gca,'xminortick','on','yminortick','on');
subplot(3,1,3);loglog(Ksat(id+1),s(:,1),'.');xlabel('Lk');ylabel('Rtidal');set(gca,'xminortick','on','yminortick','on');
Tav=sum(T(:,1:10).*Lt(:,1:10),2)./sum(Lt(:,1:10),2);
Roav=sum(Ro(:,1:10).*V(:,1:10),2)./sum(V(:,1:10),2);%Roav2=sum(n(:,1:10),2)./s(:,1).^3;
Mgas=sum(n(:,1:10),2);
% figure;plot3(Mgas,Roav,Tav,'.');set(gca,'Xscale','log','Yscale','log','zscale','log');xlabel('m');ylabel('Rho');zlabel('T');

figure;
h1=subplot(3,1,1);loglog(s(:,2),Mgas,'.');ylabel('M_{gas}');set(gca,'xminortick','on','yminortick','on');title('gas vs. DM');
h2=subplot(3,1,2);loglog(s(:,2),Tav,'.');ylabel('<T>/keV');set(gca,'xminortick','on','yminortick','on');
h3=subplot(3,1,3);loglog(s(:,2),Roav,'.');ylabel('<\rho_{gas}>');set(gca,'xminortick','on','yminortick','on');
xlabel('Mdm');
align([h1,h2,h3],'horizontalalignment','fixed',30);

figure;
h1=subplot(3,1,1);loglog(Mgas,Roav,'.');xlabel('M_{gas}');ylabel('<\rho_{gas}>');set(gca,'xminortick','on','yminortick','on');title('gas properties');
h2=subplot(3,1,2);loglog(Mgas,Tav,'.');xlabel('M_{gas}');ylabel('<T>/keV');set(gca,'xminortick','on','yminortick','on');
h3=subplot(3,1,3);loglog(Mgas.*Roav,Tav,'.');xlabel('M_{gas}*<\rho_{gas}>');set(gca,'xminortick','on','yminortick','on');ylabel('<T>/keV');

figure;
h1=subplot(3,1,1);loglog(s(:,1),Mgas,'.');ylabel('M_{gas}');set(gca,'xminortick','on','yminortick','on');title('gas vs. size');
h2=subplot(3,1,2);loglog(s(:,1),Tav,'.');ylabel('<T>/keV');set(gca,'xminortick','on','yminortick','on');
h3=subplot(3,1,3);loglog(s(:,1),Roav,'.');ylabel('<\rho_{gas}>');set(gca,'xminortick','on','yminortick','on');
xlabel('Rtidal');
align([h1,h2,h3],'horizontalalignment','fixed',30);

Lav=sum(Lt(:,1:10),2)-sum(Lt(:,11:13),2)/(1.3^3-1);
figure;
h1=subplot(3,1,1);loglog(Ksat(id+1),Mgas,'.');ylabel('M_{gas}');set(gca,'xminortick','on','yminortick','on');title('steller vs. gas');
h2=subplot(3,1,2);loglog(Ksat(id+1),Tav,'.');ylabel('<T>/keV');set(gca,'xminortick','on','yminortick','on');
h3=subplot(3,1,3);loglog(Ksat(id+1),Roav,'.');ylabel('<\rho_{gas}>');set(gca,'xminortick','on','yminortick','on');
xlabel('Lk');
align([h1,h2,h3],'horizontalalignment','fixed',30);

figure;hist(Mgas(1:100)./s(1:100,2),0:0.1:3);xlabel('Ngas/Ndm');ylabel('count');

% figure;plot3(Lav*1e40,s(:,1),s(:,2),'.');set(gca,'Xscale','log','Yscale','log','zscale','log');xlabel('L');ylabel('r');zlabel('m');
%%
nplot=1:size(Lt,1);
bi=1:size(Lt,1);
bid=id(bi);
LL=sum(Lt(bi(nplot),1:10),2);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LL=LL(find(LL>0));
rr=rratio(find(LL>0));
figure;
for j=1:size(LL,1)
    clr=min(rr(j)/5,1);
    plot(log10(KK(j)),log10(LL(j))+40,'o','markerfacecolor',[clr,clr,clr]);
    hold on;
end
Lav=sum(Lt(:,1:10),2)-sum(Lt(:,11:13),2)/(1.3^3-1);
Lbk=sum(Lt(:,11:13),2)./sum(V(:,11:13),2);%total background emissivity
Lerr=sqrt(sum((Lm(:,11:13)-repmat(Lbk,1,3)).^2,2)/3/2)*(10^3);%the std error for L=(L_inner-L_ith_backgrd)
% bi=find(Lbk>Lerr);
%bi=find(Lav>2*Lerr);
rratiobin=[0,0.3,0.6,0.9,1.5,5,inf];colors=['k','b','g','c','r','m','y'];
for j=1:size(rratiobin,2)-1
    bi=find(rratio>rratiobin(j)&rratio<rratiobin(j+1));
    bid=id(bi);
    KK=Ksat(bid+1);LL=Lav(bi);LLerr=Lerr(bi);
    LLu=LL+LLerr;LLl=LL-LLerr;
    LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
    KK=KK(find(LL>0));
    LL=LL(find(LL>0));
    le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
    errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),[colors(mod(j,7)),'.'])
    hold on;
end
% plot([10,11.5],[36.8,39]+1,':r',[11,12],[39.1,43.65]+1,':b');
plot([10,12],[38.1,41.27],'--k',[10,12],[36.89,42.37],'-k');
set(gca,'xminortick','on','yminortick','on');

Lav=sum(Lt(:,1:10),2)-sum(Lt(:,11:13),2)/(1.3^3-1);
Lbk=sum(Lt(:,11:13),2)./sum(V(:,11:13),2);%total background emissivity
Lerr=sqrt(sum((Lm(:,11:13)-repmat(Lbk,1,3)).^2,2)/3/2)*(10^3);%the std error for L=(L_inner-L_ith_backgrd)
bi=find(id<gm(2)&rratio<0.4&rratio>0.2);
bid=id(bi);
KK=Ksat(bid+1);LL=Lav(bi);LLerr=Lerr(bi);
LLu=LL+LLerr;LLl=LL-LLerr;
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
KK=KK(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'.')
hold on;
% plot([10,11.5],[36.8,39]+1,':r',[11,12],[39.1,43.65]+1,':b');
plot([10,12],[38.1,41.27],'--k',[10,12],[36.89,42.37],'-k');
set(gca,'xminortick','on','yminortick','on');
%%
nplot=1:size(Lt,1);
bi=1:size(Lt,1);
bid=id(bi);
LL=sum(Lt(:,1:10),2)-sum(Lt(:,11:15),2)/(1.5^3-1);
LLl=sum(Lt(:,1:10),2)-sum(Lt(:,11:13),2)/(1.3^3-1);
LLu=sum(Lt(:,1:10),2)-sum(Lt(:,14:15),2)/(1.5^3-1.3^3);
Larr=[LL,LLl,LLu];
LL=Larr(bi(nplot),1);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=Larr(bi(nplot),2);LLu=Larr(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
figure;errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'g.')
hold on;
% plot([10,11.5],[36.8,39]+1,':r',[11,12],[39.1,43.65]+1,':b');
plot([10,12],[38.1,41.27],'--k',[10,12],[36.89,42.37],'-k');
set(gca,'xminortick','on','yminortick','on');
%%
bi=find(sum(Lt(:,1:10),2)./(sum(Lt(:,11:15),2)/(1.5^3-1))>2);
bid=id(bi);
nplot=1:size(bid,1);
% nplot=1:92;
% nplot=1:89;
% nplot=1:43;
% LL=subLx_099(bi(nplot),4);
LL=Larr(bi(nplot),1);LLl=Larr(bi(nplot),2);LLu=Larr(bi(nplot),3);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));

% LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
figure;errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'b.')
hold on;
plot([10,11.5],[36.8,39],':r',[11,12],[39.1,43.65],':b');
plot([10,12],[38.1,41.27],'--k',[10,12],[36.89,42.37],'--k');
%%
[i,j]=find(~(repmat(id,1,size(gm,1))-repmat(gm'+1,size(id,1),1)));%find main subs [subid,grpid]
% mr=[];
% for j=1:size(gm,1)
%     if group_offset_099(j,2)>1
%         for k=(gm(j)+2):(gm(j)+group_offset_099(j,2))
%             mr=[mr;sub_mass_099(k)/sub_mass_099(gm(j)+1)];
%         end
%     end
% end

p=[];
x=log10(1:8);
for j=1:size(Lm,1)
    y=log10(Lm(j,1:8));
    if(size(find(isfinite(y)),2)<3)
        p=[p;NaN,0];
    else
    p=[p;polyfit(x(isfinite(y)),y(isfinite(y)),1)];
    end
end
pp=p(:,1);
FWHM=0.5.^(1./pp).*s(:,1)/10*2;

z=1./outputs_zoom(SnapInfall(id+1,3)+1)-1;
% plot(z,p(:,1),'.',z(isnan(pp)),1,'ro');
% figure;plot(1-SnapInfall(id+1,5)./SnapInfall(id+1,2),p(:,1),'.')
frac_lost=1-SnapInfall(id+1,5)./SnapInfall(id+1,2);
figure;semilogx(frac_lost,p(:,1),'.',frac_lost(isnan(pp)),1,'r.');
 figure;hist(p(:,1),100)

% bi=find(pp<0);pp<1&FWHM>0&
bi=find(rcen<800);%chandra psf fwhm 0.5 arcsec=0.2kpc at z=0.03; FOV (diameter) 17arcmin=500kpc at z=0.03 (factor of 0.4 smaller when z=0.01)
                  % ACIS-S FOV 51arcmin*8arcmin=1500kpc*250kpc
bid=id(bi);
ppp=pp(bi);
zz=1./outputs_zoom(SnapInfall(bid+1,3)+1)-1;
% hold on;plot(zz(1:end),ppp(1:end),'rd')
figure;semilogy(Larr(bi(1:end),:),'.-')
% figure;semilogy(subLx_099(bi(1:end),2:4),'.-')
% errorbar(Ksat(bid(1:28)+1),subLx_099_000(bid(1:28),4),subLx_099_000(bid(1:28),4)-subLx_099_000(bid(1:28),2),subLx_099_000(bid(1:28),3)-subLx_099_000(bid(1:28),4),'.');
% set(gca,'xscale','log','yscale','log');
% figure;errorbar(log10(Ksat(bid(1:15)+1)),log10(subLx_099_000(bi(1:15),4))+40,log10(subLx_099_000(bi(1:15),4))-log10(subLx_099_000(bi(1:15),2)),log10(subLx_099_000(bi(1:15),3))-log10(subLx_099_000(bi(1:15),4)),'.');

%%
rmax=500;
bi=find(rcen<rmax);%chandra psf fwhm 0.5 arcsec=0.2kpc at z=0.03; FOV (diameter) 17arcmin=500kpc at z=0.03 (factor of 0.4 smaller when z=0.01)
                  % ACIS-S FOV 51arcmin*8arcmin=1500kpc*250kpc
bid=id(bi);
colors=['k','b','g','c','r','m','y'];
% nplot=1:92;
% nplot=1:89;
figure;
h1=subplot(3,1,1);
for j=1:6
nplot=find(bid>gm(j)&bid<gm(j+1));
nplot([1,end])'
LL=Larr(bi(nplot),1);LLl=Larr(bi(nplot),2);LLu=Larr(bi(nplot),3);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),[colors(j),'.'])
hold on;
end
plot([10,12],[38.1,41.27],'-k',[10,12],[36.89,42.37],'--k');
hold on;
legend([mat2cell([num2str(group_offset_099(1:j,3)*10^10,'%1.1e'),repmat('M_{sun}/h',j,1)],ones(6,1),16);'cluster';'group'],...
        'location','northwest','color','none');
legend boxoff;
% nplot=find(bid>gm(j+1));
% nplot([1,end])'
% LL=Larr(bi(nplot),1);
% KK=Ksat(bid(nplot)+1);
% KK=KK(find(LL>0));
% LLl=Larr(bi(nplot),2);LLu=Larr(bi(nplot),3);
% LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
% LL=LL(find(LL>0));
% le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
% errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),[colors(j),'.'])
set(gca,'xminortick','on','yminortick','on');
% xlabel('log10(Lk/Lk_{sun})');
ylabel('log10(Lx/(erg/s))');title(['sim6702,within ',num2str(rmax),'kpc/h of central gals'])
hold off;

h2=subplot(3,1,2);
for j=1:6
nplot=find(bid>gm(j)&bid<gm(j+1));
nplot([1,end])'
LL=Larr(bi(nplot),1);LLl=Larr(bi(nplot),2);LLu=Larr(bi(nplot),3);
KK=Ksat(bid(nplot)+1);
ind=find(log10(KK)>10.45&LL>0&LLl>0&LLu>0);
KK=KK(ind);LL=LL(ind);
LLl=LLl(ind);LLu=LLu(ind);
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);
errorbar(log10(KK),log10(LL)+40,le,ue,[colors(j),'.'])
hold on;
end
plot([10,12],[38.1,41.27],'-k',[10,12],[36.89,42.37],'--k');
hold on;
% legend([mat2cell([num2str(group_offset_099(1:j,3)*10^10,'%1.1e'),repmat('M_{sun}/h',j,1)],ones(6,1),16);'cluster';'group']);
set(gca,'xminortick','on','yminortick','on');
axis([10,12,36.5,43])
% xlabel('log10(Lk/Lk_{sun})');
ylabel('log10(Lx/(erg/s))');text(10.2,42,{'K band luminosity limited';'             +';['central dist<',num2str(rmax),'kpc/h']});
hold off;

h3=subplot(3,1,3);
nplot=1:size(Lt,1);
bi=1:size(Lt,1);
bid=id(bi);
LL=sum(Lt(bi(nplot),1:10),2);
KK=Ksat(bid(nplot)+1);
plot(log10(KK),log10(LL)+40,'.');
hold on;
plot([10,12],[38.1,41.27],'-k',[10,12],[36.89,42.37],'--k');
hold on;
set(gca,'xminortick','on','yminortick','on');
xlabel('log10(Lk/Lk_{sun})');ylabel('log10(Lx/(erg/s))');text(6.5,43,'raw data, all satellites');
hold off;
align([h1,h2,h3],'horizontalalignment','fixed',35);
%%

nplot=1:96;
LL=Larr(bi(nplot),4);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=Larr(bi(nplot),2);LLu=Larr(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
cftool(log10(KK(ri)),log10(LL(ri))+40)
% errorbar(KK,LL*10^40,LL-LLl,LLu-LL,'.');set(gca,'xscale','log','yscale','log');

%%
nplot=1:9;
figure;
loglog(r(nplot,:)',Ro(nplot,:)','-');title('rho');
figure;
loglog(r(nplot,:)',T(nplot,:)','-');title('T');
figure;
loglog(r(nplot,:)',Lm(nplot,:)','-');title('L');
