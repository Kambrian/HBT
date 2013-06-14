runnum='6702DM';
Nsnap=99;
fact_relax=3;

file=['/mnt/A4700/data/',runnum,'/subcat/anal/subprof_dm_',num2str(Nsnap,'%03d'),'.',num2str(fact_relax),'vir'];
load(file);
% scaleF_file=['/mnt/A4700/data/',runnum,'/subcat/Redshift.dat'];
% tmp=load(scaleF_file);a=tmp(:,2);z=1./a-1;

%%

nbin=120;
fact_subr=0.02;

rs=subprof_dm_099(:,1:nbin);
ns=subprof_dm_099(:,nbin+1:nbin*2);
vs=1:nbin;
vs=vs.^3-(vs-1).^3;
vs=repmat(vs,size(ns,1),1);
ro=subprof_dm_099(:,2*nbin+1:3*nbin);
no=subprof_dm_099(:,3*nbin+1:4*nbin);
rb=subprof_dm_099(:,4*nbin+1:5*nbin);
nb=subprof_dm_099(:,5*nbin+1:6*nbin);
nplot0=[2,6,7,9,15,16,18,22,32];
figure(1);
loglog(rs(nplot0,:)'*fact_relax,ns(nplot0,:)'./vs(nplot0,:)');
% legend('HBT');
hold on;

nplot=1:10;
figure;
loglog(rs(nplot,:)'*fact_relax,ns(nplot,:)'./vs(nplot,:)')
hold on;
loglog(rb(nplot,:)'*fact_relax,nb(nplot,:)'./vs(nplot,:)','--')
loglog(ro(nplot,:)'*fact_relax,no(nplot,:)'./vs(nplot,:)',':')
nr=abs(ns-fact_subr*(nb+no));
I=[];
for i=1:length(rs)
    [v,j]=max(ns(i,:)./vs(i,:),[],2);
    [v,j2]=min(nr(i,j:end)./vs(i,j:end)./ns(i,j:end),[],2);
    I=[I;j+j2-1];
end

% [v,I]=min(abs(ns-nb-no),[],2);%this is problematic since the inner counts can be 0 and identical
for i=nplot
loglog(rs(i,I(i))*fact_relax,ns(i,I(i))/vs(i,I(i)),'o');
end
figure(1);hold on;
for i=nplot0
loglog(rs(i,I(i))*fact_relax,ns(i,I(i))/vs(i,I(i)),'o');
end

req=zeros(length(rs),1);
for i=1:length(rs)
    req(i)=rs(i,I(i));
end
rvir=subprof_dm_099(:,6*nbin+1)*fact_relax;
rcen=subprof_dm_099(:,6*nbin+2);
M=subprof_dm_099(:,6*nbin+3);
SubRank=subprof_dm_099(:,6*nbin+5);

% nbin=10;
rtidal=[];
na=ns+no+nb;
sna=cumsum(na(1,:));
% dr=rvir(1)/nbin;
for i=2:group_offset_099(1,2)
    icen=floor(rcen(i)/rvir(1)*nbin)+1;
    if(icen>nbin)
        icen=nbin
    end
    mr=sna(icen);
%     dmr=ns(1,icen)/mr/dr*rcen(i);
    dmr=na(1,icen)/vs(1,icen)*icen^3/mr;
    rtidal=[rtidal;rcen(i)*(M(i)/(3-dmr)/mr)^(1.0/3)];
end

figure;
nplot=5:600;
loglog(req(nplot).*rvir(nplot),rtidal(nplot-1),'.');
axis([1,1000,1,1000]);
hold on;loglog([1,1000],[1,1000]);hold off;
xlabel(['Req_{',num2str(fact_subr),'}']);ylabel('Rtidal');

figure;plot(rtidal(1:50),'.');
hold on;plot(req(2:51).*rvir(2:51),'ro');
figure; plot(rtidal(1:50)./req(2:51)./rvir(2:51),'.');
%%
clear;
basedir='/home/kam/Documents/research/Galaxy/code/BoundTracing/v6/sim6702/anal/data/prof/'
load([basedir,'subfindprof_dm_099'])
load([basedir,'SubFindgroup_offset_099']);
% load subprof_dm_099
% load ../../BoundTracing/v3.0/sim6702/anal/massfun/group_offset_099
subprof_dm_099=subfindprof_dm_099;

nbin=120;
fact_relax=3;
fact_subr=0.02;
% fact_subr=1;

rs=subprof_dm_099(:,1:nbin);
ns=subprof_dm_099(:,nbin+1:nbin*2);
vs=1:nbin;
vs=vs.^3-(vs-1).^3;
vs=repmat(vs,size(ns,1),1);
ro=subprof_dm_099(:,2*nbin+1:3*nbin);
no=subprof_dm_099(:,3*nbin+1:4*nbin);
rb=subprof_dm_099(:,4*nbin+1:5*nbin);
nb=subprof_dm_099(:,5*nbin+1:6*nbin);
nplot0=[3,8,51,24,18,19,20,25];
figure(1);hold on;
loglog(rs(nplot0,:)'*fact_relax,ns(nplot0,:)'./vs(nplot0,:)','--');
% legend('SubFind');

nplot=1:10;
figure;
loglog(rs(nplot,:)'*fact_relax,ns(nplot,:)'./vs(nplot,:)')
hold on;
loglog(rb(nplot,:)'*fact_relax,nb(nplot,:)'./vs(nplot,:)','--')
loglog(ro(nplot,:)'*fact_relax,no(nplot,:)'./vs(nplot,:)',':')
% snb=cumsum(nb+no,2);
% svs=cumsum(vs,2);
% nr=abs(ns-fact_subr*snb./svs.*vs);
nr=abs(ns-fact_subr*(nb+no));
I=[];
for i=1:length(rs)
    [v,j]=max(ns(i,:)./vs(i,:),[],2);
    [v,j2]=min(nr(i,j:end)./vs(i,j:end)./ns(i,j:end),[],2);
    I=[I;j+j2-1];
end

% [v,I]=min(abs(ns-nb-no),[],2);%this is problematic since the inner counts can be 0 and identical
for i=nplot
loglog(rs(i,I(i))*fact_relax,ns(i,I(i))/vs(i,I(i)),'o');
end

figure(1);hold on;
for i=nplot0
loglog(rs(i,I(i))*fact_relax,ns(i,I(i))/vs(i,I(i)),'o');
end

req=zeros(length(rs),1);
for i=1:length(rs)
    req(i)=rs(i,I(i));
end
rvir=subprof_dm_099(:,6*nbin+1)*fact_relax;
rcen=subprof_dm_099(:,6*nbin+2);
M=subprof_dm_099(:,6*nbin+3);
SubRank=subprof_dm_099(:,6*nbin+5);

% nbin=10;
rtidal=[];
na=ns+no+nb;
sna=cumsum(na(1,:));
% dr=rvir(1)/nbin;
for i=2:SubFindgroup_offset_099(1,1)
    icen=floor(rcen(i)/rvir(1)*nbin)+1;
    if(icen>nbin)
        icen=nbin
    end
    mr=sna(icen);
%     dmr=ns(1,icen)/mr/dr*rcen(i);
    dmr=na(1,icen)/vs(1,icen)*icen^3/mr;
    rtidal=[rtidal;rcen(i)*(M(i)/(3-dmr)/mr)^(1.0/3)];
end

figure;
nplot=2:SubFindgroup_offset_099(1,1);
loglog(rtidal(nplot-1),req(nplot).*rvir(nplot),'.');
axis([1,10000,1,10000]);
hold on;loglog([1,10000],[1,10000]);hold off;
ylabel(['Req_{',num2str(fact_subr),'}']);xlabel('Rtidal');

figure;plot(rtidal(1:50),'.');
hold on;plot(req(2:51).*rvir(2:51),'ro');
figure; plot(rtidal(1:50)./req(2:51)./rvir(2:51),'.');
%%
load subprof_dm_099
load ../../BoundTracing/v3.0/sim6702/anal/massfun/group_offset_099

nbin=20;
fact_relax=1;
fact_subr=0.02;

rs=subprof_dm_099(:,1:nbin);
ns=subprof_dm_099(:,nbin+1:nbin*2);
vs=1:nbin;
vs=vs.^3-(vs-1).^3;
vs=repmat(vs,size(ns,1),1);
ro=subprof_dm_099(:,2*nbin+1:3*nbin);
no=subprof_dm_099(:,3*nbin+1:4*nbin);
rb=subprof_dm_099(:,4*nbin+1:5*nbin);
nb=subprof_dm_099(:,5*nbin+1:6*nbin);
nplot0=[2,6,30,7,14,24,15,20,32];
figure(1);hold on;
loglog(rs(nplot0,:)'*fact_relax,ns(nplot0,:)'./vs(nplot0,:)'/8,':');
% legend('BT3');
nplot=1:10;
figure;
loglog(rs(nplot,:)'*fact_relax,ns(nplot,:)'./vs(nplot,:)')
hold on;
loglog(rb(nplot,:)'*fact_relax,nb(nplot,:)'./vs(nplot,:)','--')
loglog(ro(nplot,:)'*fact_relax,no(nplot,:)'./vs(nplot,:)',':')
nr=abs(ns-fact_subr*(nb+no));
I=[];
for i=1:length(rs)
    [v,j]=max(ns(i,:)./vs(i,:),[],2);
    [v,j2]=min(nr(i,j:end)./vs(i,j:end)./ns(i,j:end),[],2);
    I=[I;j+j2-1];
end

% [v,I]=min(abs(ns-nb-no),[],2);%this is problematic since the inner counts can be 0 and identical
for i=nplot
loglog(rs(i,I(i))*fact_relax,ns(i,I(i))/vs(i,I(i)),'o');
end
figure(1);hold on;
for i=nplot0
loglog(rs(i,I(i))*fact_relax,ns(i,I(i))/vs(i,I(i))/8,'o');
end

req=zeros(length(rs),1);
for i=1:length(rs)
    req(i)=rs(i,I(i));
end
rvir=subprof_dm_099(:,6*nbin+1)*fact_relax;
rcen=subprof_dm_099(:,6*nbin+2);
M=subprof_dm_099(:,6*nbin+3);
SubRank=subprof_dm_099(:,6*nbin+5);

% nbin=10;
rtidal=[];
na=ns+no+nb;
sna=cumsum(na(1,:));
% dr=rvir(1)/nbin;
for i=2:group_offset_099(1,2)
    icen=floor(rcen(i)/rvir(1)*nbin)+1;
    if(icen>nbin)
        icen=nbin
    end
    mr=sna(icen);
%     dmr=ns(1,icen)/mr/dr*rcen(i);
    dmr=na(1,icen)/vs(1,icen)*icen^3/mr;
    rtidal=[rtidal;rcen(i)*(M(i)/(3-dmr)/mr)^(1.0/3)];
end

figure;
nplot=5:600;
loglog(req(nplot).*rvir(nplot),rtidal(nplot-1),'.');
axis([1,1000,1,1000]);
hold on;loglog([1,1000],[1,1000]);hold off;
xlabel(['Req_{',num2str(fact_subr),'}']);ylabel('Rtidal');

figure;plot(rtidal(1:50),'.');
hold on;plot(req(2:51).*rvir(2:51),'ro');
figure; plot(rtidal(1:50)./req(2:51)./rvir(2:51),'.');
%%
load /home/kam/Documents/research/Galaxy/code/BoundTracing/v7.1/data/prof/subprof_dm_099.3vir
load /home/kam/Documents/research/Galaxy/code/BoundTracing/v7.1/data/massfun/v66/group_offset_099

nbin=120;
fact_relax=3;
fact_subr=0.02;

rs=subprof_dm_099(:,1:nbin);
ns=subprof_dm_099(:,nbin+1:nbin*2);
vs=1:nbin;
vs=vs.^3-(vs-1).^3;
vs=repmat(vs,size(ns,1),1);
ro=subprof_dm_099(:,2*nbin+1:3*nbin);
no=subprof_dm_099(:,3*nbin+1:4*nbin);
rb=subprof_dm_099(:,4*nbin+1:5*nbin);
nb=subprof_dm_099(:,5*nbin+1:6*nbin);
nplot0=[2,6,30,7,14,24,15,20,32];
figure(1);
loglog(rs(nplot0,:)'*fact_relax,ns(nplot0,:)'./vs(nplot0,:)'/8,':');hold on;
% legend('BT3');
nplot=1:10;%[42,86,104,109,118,165];
figure;
loglog(rs(nplot,:)'*fact_relax,ns(nplot,:)'./vs(nplot,:)','- ')
hold on;
loglog(rb(nplot,:)'*fact_relax,nb(nplot,:)'./vs(nplot,:)','-- ')
loglog(ro(nplot,:)'*fact_relax,no(nplot,:)'./vs(nplot,:)',': ')
nr=abs(ns-fact_subr*(nb+no));
I=[];
for i=1:length(rs)
    [v,j]=max(ns(i,:)./vs(i,:),[],2);
    [v,j2]=min(nr(i,j:end)./vs(i,j:end)./ns(i,j:end),[],2);
    I=[I;j+j2-1];
end

% [v,I]=min(abs(ns-nb-no),[],2);%this is problematic since the inner counts can be 0 and identical
for i=nplot
loglog(rs(i,I(i))*fact_relax,ns(i,I(i))/vs(i,I(i)),'o');
end
figure(1);hold on;
for i=nplot0
loglog(rs(i,I(i))*fact_relax,ns(i,I(i))/vs(i,I(i))/8,'o');
end

req=zeros(length(rs),1);
for i=1:length(rs)
    req(i)=rs(i,I(i));
end
rvir=subprof_dm_099(:,6*nbin+1)*fact_relax;
rcen=subprof_dm_099(:,6*nbin+2);
M=subprof_dm_099(:,6*nbin+3);
SubRank=subprof_dm_099(:,6*nbin+5);

% nbin=10;
rtidal=[];
na=ns+no+nb;
sna=cumsum(na(1,:));
% dr=rvir(1)/nbin;
for i=2:group_offset_099(1,2)
    icen=floor(rcen(i)/rvir(1)*nbin)+1;
    if(icen>nbin)
        icen=nbin
    end
    mr=sna(icen);
%     dmr=ns(1,icen)/mr/dr*rcen(i);
    dmr=na(1,icen)/vs(1,icen)*icen^3/mr;
    rtidal=[rtidal;rcen(i)*(M(i)/(3-dmr)/mr)^(1.0/3)];
end

figure;
nplot=5:600;
loglog(req(nplot).*rvir(nplot),rtidal(nplot-1),'.');
axis([1,1000,1,1000]);
hold on;loglog([1,1000],[1,1000]);hold off;
xlabel(['Req_{',num2str(fact_subr),'}']);ylabel('Rtidal');

figure;plot(rtidal(1:50),'.');
hold on;plot(req(2:51).*rvir(2:51),'ro');
figure; plot(rtidal(1:50)./req(2:51)./rvir(2:51),'.');

