clear;
addpath('../post');

% runnum='8213';NsnapSF=59;NsnapBT=59;
runnumA='6702DM';NsnapA=99;
runnumB='6702DMNSM';NsnapB=99;

file=['/mnt/A4700/data/',runnumA,'/subcat/anal/MatchTo/',runnumB,'/',num2str(NsnapA),'_',num2str(NsnapB)];
BT2NSM=load(file);
BT2NSM=BT2NSM(1:6689,:);%select first halo
% BT2SF(BT2SF(:,4)==0,:)=[];%exclude fake subhalos
% BT2SF(BT2SF(:,4)<100,:)=[];



file=['/mnt/A4700/data/',runnumB,'/subcat/anal/MatchTo/',runnumA,'/',num2str(NsnapB),'_',num2str(NsnapA)];
NSM2BT=load(file);
NSM2BT=NSM2BT(1:6687,:);
% NSM2BT(NSM2BT(:,4)<100,:)=[];
%%
%select 1-1 match
f=zeros(size(BT2NSM,1),1);
for i=1:size(BT2NSM,1)
    j=BT2NSM(i,2)+1;
    if j>0&&NSM2BT(j,2)+1==i
        f(i)=1;
    end
end
f=logical(f);

% f=BT2NSM(:,2)>0;
figure;
loglog(BT2NSM(f,4),BT2NSM(f,5),'o');

figure;
semilogx(BT2NSM(f,4),BT2NSM(f,5)./BT2NSM(f,4),'o');
%%
%select 1-1 match
f=zeros(size(NSM2BT,1),1);
for i=1:size(NSM2BT,1)
    j=NSM2BT(i,2)+1;
    if j>0&&BT2NSM(j,2)+1==i
        f(i)=1;
    end
end
f=logical(f);

% f=NSM2BT(:,2)>0;
figure;
loglog(NSM2BT(f,4),NSM2BT(f,5),'o');

figure;
semilogx(NSM2BT(f,4),NSM2BT(f,5)./NSM2BT(f,4),'o');