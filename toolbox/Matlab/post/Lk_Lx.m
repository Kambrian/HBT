%%
clear;
% sim='6702/v5';
% load /home/kam/Documents/research/Galaxy/code/gas/outputs_zoom.txt
% uiimport([sim,'/SnapInfall'])
% load([sim,'/subprof_099'])
% load([sim,'/group_offset_099'])
% load([sim,'/sub_mass_099'])
% load([sim,'/subLx_099'])
sim='/home/kam/Documents/research/Galaxy/code/BoundTracing/v6/sim6702/anal/data';
load /home/kam/Documents/research/Galaxy/code/gas/outputs_zoom.txt
uiimport([sim,'/steller/SnapInfall'])
load([sim,'/prof/LxTidal/subprof_099'])
load([sim,'/massfun/group_offset_099'])
load([sim,'/massfun/sub_mass_099'])
load([sim,'/prof/LxTidal/subLx_099'])
%%
subLx_099=subLx_099*10;
L=subprof_099(:,1:10)*10;

gm=group_offset_099(:,1);
Ro=subprof_099(:,11:20);
T=subprof_099(:,21:30);
r=subprof_099(:,31:40);
n=subprof_099(:,41:50);
s=subprof_099(:,51:52);

Lt=L.*n;
x=0:10;
V=diff(x.^3);
V=repmat(V,size(L,1),1);
Lm=Lt./V;
id=subprof_099(:,53);

steller;
%%
nplot=1:size(Lt,1);
bi=1:size(Lt,1);
bid=id(bi);
LL=subLx_099(bi(nplot),4);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
figure;errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'g.')
hold on;
plot([10,11.5],[36.8,39]+1,':r',[11,12],[39.1,43.65]+1,':b');
plot([10,12],[38.1,41.27],'--k',[10,12],[36.89,42.37],'--k');
set(gca,'xminortick','on','yminortick','on');
%%
[i,j]=find(~(repmat(id,1,size(gm,1))-repmat(gm'+1,size(id,1),1)));%find main subs [subid,grpid]

mr=[];
for j=1:size(gm,1)
    if group_offset_099(j,2)>1
        for k=(gm(j)+2):(gm(j)+group_offset_099(j,2))
            mr=[mr;sub_mass_099(k)/sub_mass_099(gm(j)+1)];
        end
    end
end

p=[];
for j=1:size(L,1)
    p=[p;polyfit(log10(1:8),log10(Lm(j,1:8)),1)];
end
pp=p(:,1);

z=1./outputs_zoom(SnapInfall(id+1,3)+1)-1;
% plot(z,p(:,1),'.',z(isnan(pp)),1,'ro');
% figure;plot(1-SnapInfall(id+1,5)./SnapInfall(id+1,2),p(:,1),'.')
frac_lost=1-SnapInfall(id+1,5)./SnapInfall(id+1,2);
figure;semilogx(frac_lost,p(:,1),'.',frac_lost(isnan(pp)),1,'r.');
 figure;hist(p(:,1),50)

bi=find(pp<0);
bid=id(bi);
ppp=pp(bi);
zz=1./outputs_zoom(SnapInfall(bid+1,3)+1)-1;
% hold on;plot(zz(1:end),ppp(1:end),'rd')

figure;semilogy(subLx_099(bi(1:end),2:4),'.-')
% errorbar(Ksat(bid(1:28)+1),subLx_099_000(bid(1:28),4),subLx_099_000(bid(1:28),4)-subLx_099_000(bid(1:28),2),subLx_099_000(bid(1:28),3)-subLx_099_000(bid(1:28),4),'.');
% set(gca,'xscale','log','yscale','log');
% figure;errorbar(log10(Ksat(bid(1:15)+1)),log10(subLx_099_000(bi(1:15),4))+40,log10(subLx_099_000(bi(1:15),4))-log10(subLx_099_000(bi(1:15),2)),log10(subLx_099_000(bi(1:15),3))-log10(subLx_099_000(bi(1:15),4)),'.');

% nplot=1:92;
nplot=1:89;
% nplot=1:51;
LL=subLx_099(bi(nplot),4);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
figure;errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'k.')
hold on;
% nplot=93:107;
nplot=90:103;
% nplot=52:63;
LL=subLx_099(bi(nplot),4);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'b.')
% nplot=108:120;
nplot=104:115
% nplot=64:71;
LL=subLx_099(bi(nplot),4);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'r.')
% nplot=121:124;
nplot=116:117;
% nplot=72:74;
LL=subLx_099(bi(nplot),4);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'m.')
% nplot=125:141;
nplot=118:139;
% nplot=75:77;
LL=subLx_099(bi(nplot),4);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
LL=LL(find(LL>0));
le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'g.')
% 
% nplot=78:82;
% LL=subLx_099(bi(nplot),4);
% KK=Ksat(bid(nplot)+1);
% KK=KK(find(LL>0));
% LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
% LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
% LL=LL(find(LL>0));
% le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
% errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'c.')
% 
% nplot=83:96;
% LL=subLx_099(bi(nplot),4);
% KK=Ksat(bid(nplot)+1);
% KK=KK(find(LL>0));
% LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
% LLl=LLl(find(LL>0));LLu=LLu(find(LL>0));
% LL=LL(find(LL>0));
% le=log10(LL)-log10(LLl);ue=log10(LLu)-log10(LL);ri=find(~imag(le)&~imag(ue));
% errorbar(log10(KK(ri)),log10(LL(ri))+40,le(ri),ue(ri),'y.')
% hold on;
plot([10,11.5],[36.8,39],':r',[11,12],[39.1,43.65],':b');
plot([10,12],[38.1,41.27],'--k',[10,12],[36.89,42.37],'--k');
set(gca,'xminortick','on','yminortick','on');
%%

nplot=1:96;
LL=subLx_099(bi(nplot),4);
KK=Ksat(bid(nplot)+1);
KK=KK(find(LL>0));
LLl=subLx_099(bi(nplot),2);LLu=subLx_099(bi(nplot),3);
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
