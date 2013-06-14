%%
clear;
sim='6702/v5';
load outputs_zoom.txt
uiimport([sim,'/SnapInfall'])
load([sim,'/gasprof_099'])
load([sim,'/group_offset_099'])
load([sim,'/sub_mass_099'])
%%
nbin=10;
gm=group_offset_099(:,1);
L=gasprof_099(:,1:nbin);
Ro=gasprof_099(:,nbin+1:2*nbin);
T=gasprof_099(:,2*nbin+1:3*nbin);
r=gasprof_099(:,3*nbin+1:4*nbin);
n=gasprof_099(:,4*nbin+1:5*nbin);
m=gasprof_099(:,5*nbin+1:5*nbin+2);
Lt=L.*n;
x=0:nbin;
V=diff(x.^3);
V=repmat(V,size(L,1),1);
Lm=Lt./V;
rho=n./V;
Lx=sum(Lt,2);

steller;
%%

nplot=1:10;
figure;
loglog(r(nplot,:)',rho(nplot,:)','-');title('rho_cat');
figure;
loglog(r(nplot,:)',T(nplot,:)','-');title('T_cat');
figure;
loglog(r(nplot,:)',Lm(nplot,:)','-');title('L-cat');

figure;
LL=Lx(1:length(Ksat))*1e40;
KK=Ksat(find(LL>0));
LL=LL(find(LL>0));
plot(log10(KK),log10(LL),'.');
cftool(log10(KK),log10(LL));

LLs=[];
KKs=[];
Ngroups=size(group_offset_099,1);
for grpid=1:Ngroups
    mainsub=group_offset_099(grpid,1)+1;
    group_len=group_offset_099(grpid,2);
for i=mainsub+1:mainsub+group_len-1
    LLs=[LLs;Lx(i)];
    KKs=[KKs;Ksat(i)];
end
end
KK=KKs(find(LLs>0));
LL=LLs(find(LLs>0));
plot(log10(KK),log10(LL)+40,'.');
cftool(log10(KK),log10(LL)+40);

