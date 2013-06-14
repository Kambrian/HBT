load subprof_dm_099.3vir
load ../massfun/v66/group_offset_099
nbin=120;
fact_relax=3;
fact_subr=1;

rs=subprof_dm_099(:,1:nbin);
ns=subprof_dm_099(:,nbin+1:nbin*2);
vs=1:nbin;
vs=vs.^3-(vs-1).^3;
vs=repmat(vs,size(ns,1),1);
ro=subprof_dm_099(:,2*nbin+1:3*nbin);
no=subprof_dm_099(:,3*nbin+1:4*nbin);
rb=subprof_dm_099(:,4*nbin+1:5*nbin);
nb=subprof_dm_099(:,5*nbin+1:6*nbin);
rvir=subprof_dm_099(:,6*nbin+1);
rscale=rvir*fact_relax;
rcen=subprof_dm_099(:,6*nbin+2);
M=subprof_dm_099(:,6*nbin+3);
SubRank=subprof_dm_099(:,6*nbin+5);

rtidal=[];outlier=[];badsmp=[];
na=ns+no+nb;
Ngroups=size(group_offset_099,1);
for grpid=1:Ngroups
    mainsub=group_offset_099(grpid,1)+1;
    group_len=group_offset_099(grpid,2);
    sna=cumsum(na(mainsub,:));
    rtidal=[rtidal; rvir(mainsub)];
for i=mainsub+1:mainsub+group_len-1
    icen=floor(rcen(i)/rscale(mainsub)*nbin)+1;
    if(icen>nbin)
        icen=nbin;
        outlier=[outlier;i,mainsub];
    end
    mr=sna(icen);
    dmr=na(mainsub,icen)/vs(mainsub,icen)*icen^3/mr;
    if (dmr>3)%probably a badly sampled host or a quite loose host
        rtidal=[rtidal;rvir(i)];
        badsmp=[badsmp;i,mainsub];
    else    
        rtidal=[rtidal;rcen(i)*(M(i)/(3-dmr)/mr)^(1.0/3)];
    end
end
end

%%

rtidal=[rtidal;rvir(length(rtidal)+1:end)];
save tidal_radius_099 rtidal -ascii

%%

% plot(rtidal(outlier(:,1)),rvir(outlier(:,1)),'.');
nplot=1:10;
% nplot=[2,6,7,9,15,16,18,22,32];%[7,8,10];%[2,7];%,621,893,915,951];
ncen=[];incen=[];
for i=nplot
icen=floor(rtidal(i)/rscale(i)*nbin)+1;
if(icen>nbin)
    icen=nbin;
end
incen=[incen,icen];
ncen=[ncen,sub2ind(size(rs),i,icen)];
% ncen=[ncen;icen];
% loglog([rtidal(i)/rvir(i),rtidal(i)/rvir(i)],[0.1,1]);
end
sns=cumsum(ns,2);
% figure;
% loglog(rs(nplot,:)'*fact_relax.*repmat(rvir(nplot)'./rtidal(nplot)',size(rs,2),1),ns(nplot,:)'./vs(nplot,:)'./repmat(sns(ncen)./incen.^3,size(rs,2),1),'-')
% hold on;
% loglog(rb(nplot,:)'*fact_relax.*repmat(rvir(nplot)'./rtidal(nplot)',size(rs,2),1),nb(nplot,:)'./vs(nplot,:)'./repmat(sns(ncen)./incen.^3,size(rs,2),1),'--')
% loglog(ro(nplot,:)'*fact_relax.*repmat(rvir(nplot)'./rtidal(nplot)',size(rs,2),1),no(nplot,:)'./vs(nplot,:)'./repmat(sns(ncen)./incen.^3,size(rs,2),1),':')

% figure;
% loglog(rs(nplot,:)'*fact_relax,ns(nplot,:)'./vs(nplot,:)'./repmat(sns(ncen)./incen.^3,size(rs,2),1),'-')
% hold on;
% loglog(rb(nplot,:)'*fact_relax,nb(nplot,:)'./vs(nplot,:)'./repmat(sns(ncen)./incen.^3,size(rs,2),1),'--')
% loglog(ro(nplot,:)'*fact_relax,no(nplot,:)'./vs(nplot,:)'./repmat(sns(ncen)./incen.^3,size(rs,2),1),':')

figure;
loglog(rs(nplot,:)'*fact_relax,ns(nplot,:)'./vs(nplot,:)','-')
hold on;
loglog(rb(nplot,:)'*fact_relax,nb(nplot,:)'./vs(nplot,:)','--')
loglog(ro(nplot,:)'*fact_relax,no(nplot,:)'./vs(nplot,:)',':')
figure;
slp=[];
for i=nplot
slp=[slp;gradient(log(smooth(rs(i,:),ns(i,:)./vs(i,:),3))',log(rs(i,:)))];
end
% semilogx(rs(nplot,1:end-1)'*fact_relax,diff(log(ns(nplot,:)'./vs(nplot,:)'),1,1)./diff(log(rs(nplot,:)'),1,1),'-')
loglog(rs(nplot,:)'*fact_relax,-1*(slp'),'-')
hold on;
n0=ns(1,:)+no(1,:)+nb(1,:);
slp0=gradient(log(smooth(rs(1,:),n0./vs(1,:),3))',log(rs(1,:)));
loglog(rs(1,:)*fact_relax,-1*slp0,':');
loglog([0.01,100],[2,2]);

loglog(rs(1,:)*fact_relax,n0./vs(1,:),'-');
loglog(rb(nplot,1:end)'*fact_relax,nb(nplot,:)'./vs(nplot,:)','--')
loglog(ro(nplot,1:end)'*fact_relax,no(nplot,:)'./vs(nplot,:)',':')

% figure;
% loglog(rs(nplot,:)'*fact_relax,ns(nplot,:)'./vs(nplot,:)'./repmat(sns(nplot,5)',size(rs,2),1),'-')
% hold on;
% loglog(rb(nplot,:)'*fact_relax,nb(nplot,:)'./vs(nplot,:)'./repmat(sns(nplot,5)',size(rs,2),1),'--')
% loglog(ro(nplot,:)'*fact_relax,no(nplot,:)'./vs(nplot,:)'./repmat(sns(nplot,5)',size(rs,2),1),':')