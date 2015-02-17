clear;
addpath('../post');

runnum='8213';NsnapSF=59;NsnapBT=59;
% runnum='6702DM';Nsnap=99;

file=['/mnt/A4700/data/',runnum,'/subcat/anal/MatchTo/SUBFIND/',num2str(NsnapBT),'_',num2str(NsnapSF)];
BT2SF=load(file);
% BT2SF=BT2SF
BT2SF=BT2SF(1:19845,:);%cut first 100 group for 8213
BT2SF(BT2SF(:,4)==0,:)=[];%exclude fake subhalos
% BT2SF(BT2SF(:,4)<100,:)=[];



file=['/mnt/A4700/data/',runnum,'/subcatS/anal/MatchTo/HBT/',num2str(NsnapSF),'_',num2str(NsnapBT)];
SF2BT=load(file);
% SF2BT(SF2BT(:,4)<100,:)=[];
sfmatchfrac=zeros(size(SF2BT,1),1);
sfmatchfrac(SF2BT(:,2)>=0)=SF2BT(SF2BT(:,2)>=0,3)./SF2BT(SF2BT(:,2)>=0,4);
%%
datadir=['/mnt/A4700/data/',runnum,'/subcat/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['main_size_',num2str(NsnapSF,'%03d')]);
btmsize=load_subsize(sizefile);

datadir=['/mnt/A4700/data/',runnum,'/subcatS/profile/'];
basedir='logbin/aqua';
sizefile=fullfile(datadir,basedir,['main_size_',num2str(NsnapSF,'%03d')]);
sfmsize=load_subsize(sizefile);
%% absorption statistics
bm=sortrows(BT2SF,2);
absorb=[];
for i=1:100
    absorb=[absorb;bm(bm(:,2)==sfmsize.id(i)&bm(:,1)~=btmsize.id(i),:)]; %satellites matched to a central
end

bms=sortrows(SF2BT,2);
absorb2=[];
for i=1:100
    absorb2=[absorb2;bms(bms(:,2)==btmsize.id(i)&bms(:,1)~=sfmsize.id(i),:)]; %exclude central
end

for i=1:100
    satid=bms(bms(:,1)==sfmsize.id(i),2);% the subhalo matched by a central
    absorb(absorb(:,1)==satid,:)=[];%kick this out
end
for i=1:100
    satid=bm(bm(:,1)==btmsize.id(i),2);% the subhalo matched by a central
    absorb2(absorb2(:,1)==satid,:)=[];%kick this out
end

x=logspace(1,5,10);
y=histc(absorb(:,4),x);dn=histc(BT2SF(:,4),x);
y2=histc(absorb2(:,4),x);dn2=histc(SF2BT(:,4),x);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
axes('position',[0.25,0.55,0.7,0.4]);
stairs(x,y,'r-','displayname','HBT');hold on;stairs(x,y2,'k-.','displayname','SUBFIND');legend('location','northeast');
set(gca,'ylim',[1,5e3],'ytick',[1;10;100;1000]);
set(gca,'xscale','log','yscale','log','xticklabel','');
ylabel('Counts');
axes('position',[0.25,0.15,0.7,0.4]);
stairs(x,y./dn,'r-','displayname','HBT');hold on;stairs(x,y2./dn2,'k-.','displayname','SUBFIND');%legend('location','northeast');
set(gca,'ytick',[1e-3;1e-2;1e-1],'ylim',[1e-3,1]);
set(gca,'xscale','log','yscale','log');
xlabel('Mass (particles)');ylabel('Relative Counts');
% set(gcf,'paperposition',[0.6,6,20,20]);
outputdir='/home/kam/Projects/HBT/code/data/show';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/absorption_',runnum,'_',num2str(NsnapSF),'.eps']);
%%
% absorb2=sortrows(absorb2,-4);
% absorb2(1:10,:)
% file=['/mnt/A4700/data/',runnum,'/subcatS/anal/MatchTo/HBT/',num2str(NsnapSF),'_',num2str(39)];
% SF2BT2=load(file);
% absorb3=zeros(size(absorb2,1),5);
% for i=1:size(absorb2,1)
%     absorb3(i,:)=SF2BT2(SF2BT2(:,1)==absorb2(i,1),:);
% end
% absorb3=sortrows(absorb3,-4);
% absorb3(1:10,:)
%%  multimatch statistics
bm=sortrows(BT2SF,2);
bm(bm(:,2)<0,:)=[];%exclude background
for i=1:100
    bm(bm(:,2)==sfmsize.id(i),:)=[]; %exclude central absorption
end
id=bm(1,2);
flag=0;
nto=0;
nfr=0;
mdbt=[];
for i=2:size(bm,1)
    if bm(i,2)~=id
       if flag
           flag=0;
           nto=nto+1;
       end
       id=bm(i,2);
    else
        flag=1;
        nfr=nfr+1;
        mdbt=[mdbt;bm(i,4)];%multi-match mass of the minor subs
    end
end
nfr=nfr+nto;
ndup_btfr=nfr;
ndup_btto=nto;
%% multi match statistics
bms=sortrows(SF2BT,2);
bms(bms(:,2)<0,:)=[];
mdsf=[];
for i=1:100
    bms(bms(:,2)==btmsize.id(i),:)=[]; %exclude central
end
id=bms(1,2);
flag=0;
nto=0;
nfr=0;
for i=2:size(bms,1)
    if bms(i,2)~=id
       if flag
           flag=0;
           nto=nto+1;
       end
       id=bms(i,2);
    else
        flag=1;
        nfr=nfr+1;
        mdsf=[mdsf;bms(i,4)];
    end
end
nfr=nfr+nto;
ndup_sffr=nfr;
ndup_sfto=nto;
%% match quality statistics
bm=sortrows(BT2SF,2);
bms=sortrows(SF2BT,2);

for i=1:100
    satid=bms(bms(:,1)==sfmsize.id(i),2);% satellite matched by SF central
    bm(bm(:,2)==sfmsize.id(i)&bm(:,1)~=btmsize.id(i)&bm(:,1)~=satid,:)=[]; %exclude central absorption of satellites
end
btmatchfrac=bm(:,3)./bm(:,4);
btmatchfrac(bm(:,2)<0)=0;%un-matched


for i=1:100
    satid=bm(bm(:,1)==btmsize.id(i),2);% the subhalo matched by a BT central
    bms(bms(:,2)==btmsize.id(i)&bms(:,1)~=sfmsize.id(i)&bms(:,1)~=satid,:)=[]; %exclude central absorption of sat
end
sfmatchfrac=bms(:,3)./bms(:,4);%central absorption has been excluded
sfmatchfrac(bms(:,2)<0)=0;%un-matched

x=0:0.02:1.05;
y=histc(btmatchfrac,x);
ys=cumsum(y(end:-1:1));ys=ys(end:-1:1);
y2=histc(sfmatchfrac,x);
ys2=cumsum(y2(end:-1:1));ys2=ys2(end:-1:1);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
set(gcf,'DefaultAxesXLim',[0,1.05]);
axes('position',[.2,.15,0.32,0.8]);
stairs(x,y,'r-','displayname','HBT');hold on;stairs(x,y2,'k-.','displayname','SUBFIND');legend('location','northwest');
xlabel('$M_{shared}/M_{sub}$');ylabel('Counts');

axes('position',[0.53,0.15,0.32,0.8]);
stairs(x,ys/ys(1),'r-','displayname','HBT');
hold on;
stairs(x,ys2/ys2(1),'k-.','displayname','SUBFIND');
xlabel('$M_{shared}/M_{sub}$');ylabel('$N(>M_{shared}/M_{sub})/N_{all}$');%legend('location','southwest');
set(gca,'yaxislocation','right');
% set(gcf,'paperposition',[0.6,6,20,20]);
outputdir='/home/kam/Projects/HBT/code/data/show';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/submatch_',runnum,'_',num2str(NsnapSF),'.eps']);

x=0:10:100;
y=histc(BT2SF(BT2SF(:,2)<0,4),x);dn=histc(BT2SF(:,4),x);
y2=histc(SF2BT(SF2BT(:,2)<0,4),x);dn2=histc(SF2BT(:,4),x);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
axes('position',[0.25,0.55,0.7,0.4]);
stairs(x,y,'r-','displayname','HBT');hold on;stairs(x,y2,'k-.','displayname','SUBFIND');legend('location','northeast');
ylabel('Counts');
set(gca,'xticklabel','','ytick',[100;200;300;400]);
axes('position',[0.25,0.15,0.7,0.4]);
stairs(x,y./dn,'r-','displayname','HBT');hold on;stairs(x,y2./dn,'k-.','displayname','SUBFIND');
xlabel('Mass (particles)');ylabel('Relative Counts');
% set(gcf,'paperposition',[0.6,6,20,20]);
outputdir='/home/kam/Projects/HBT/code/data/show';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
% print('-depsc',[outputdir,'/subunmatch_',runnum,'_',num2str(NsnapSF),'.eps']);

%% after filtering
BT2SF(BT2SF(:,4)<100,:)=[];
SF2BT(SF2BT(:,4)<100,:)=[];
bm=sortrows(BT2SF,2);
bms=sortrows(SF2BT,2);

for i=1:100
    satid=bms(bms(:,1)==sfmsize.id(i),2);% satellite matched by SF central
    bm(bm(:,2)==sfmsize.id(i)&bm(:,1)~=btmsize.id(i)&bm(:,1)~=satid,:)=[]; %exclude central absorption of satellites
end
btmatchfrac=bm(:,3)./bm(:,4);
btmatchfrac(bm(:,2)<0)=0;%un-matched


for i=1:100
    satid=bm(bm(:,1)==btmsize.id(i),2);% the subhalo matched by a BT central
    bms(bms(:,2)==btmsize.id(i)&bms(:,1)~=sfmsize.id(i)&bms(:,1)~=satid,:)=[]; %exclude central absorption of sat
end
sfmatchfrac=bms(:,3)./bms(:,4);%central absorption has been excluded
sfmatchfrac(bms(:,2)<0)=0;%un-matched

x=0:0.05:1.02;
y=histc(btmatchfrac,x);
ys=cumsum(y(end:-1:1));ys=ys(end:-1:1);
y2=histc(sfmatchfrac,x);
ys2=cumsum(y2(end:-1:1));ys2=ys2(end:-1:1);
figure;
set(gcf,...
    'DefaultLineLineWidth',1,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on',...
    'DefaultTextInterpreter','latex');
set(gcf,'DefaultAxesXLim',[0,1.05]);
axes('position',[.2,.15,0.32,0.8]);
stairs(x,y,'r--','displayname','HBT');hold on;stairs(x,y2,'k-','displayname','SUBFIND');legend('location','northwest');
xlabel('$M_{shared}/M_{sub}$');ylabel('Counts');
set(gca,'yscale','log');
axes('position',[0.53,0.15,0.32,0.8]);
stairs(x,ys/ys(1),'r--','displayname','HBT');
hold on;
stairs(x,ys2/ys2(1),'k-','displayname','SUBFIND');
xlabel('$M_{shared}/M_{sub}$');ylabel('$N(>M_{shared}/M_{sub})/N_{all}$');%legend('location','southwest');
set(gca,'yaxislocation','right');
% set(gcf,'paperposition',[0.6,6,20,20]);
outputdir='/home/kam/Projects/HBT/code/data/show';
% outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
print('-depsc',[outputdir,'/submatch100_',runnum,'_',num2str(NsnapSF),'.eps']);
%%
RunNum='8213';Nsnap=59;SofteningHalo=3.3;
datadir=['/mnt/A4700/data/',RunNum,'/subcat/anal/'];
subpos=load([datadir,'subcen_',num2str(Nsnap,'%03d')]);  
datadir=['/mnt/A4700/data/',RunNum,'/subcatS/anal/'];
sfpos=load([datadir,'subcen_',num2str(Nsnap,'%03d')]);  
%%
bm=sortrows(BT2SF(BT2SF(:,4)>=100,:),2);
bm(bm(:,2)<0,:)=[];%exclude background
for i=1:100
    bm(bm(:,2)==sfmsize.id(i),:)=[]; %exclude central
end
id=bm(1,2);
flag=0;
nto=0;
nfr=0;
mdbt=[];
pair=[];
ifrom=1;
idsf=bms(bms(:,2)==bm(1,1),1);
if numel(idsf)>1
    idsf=idsf(1);
end
for i=2:size(bm,1)
    if bm(i,2)~=id  %new
       if flag
           flag=0;
           nto=nto+1;
       end
       id=bm(i,2);
       ifrom=i;
       idsf=bms(bms(:,2)==bm(i,1),1);
       if numel(idsf)>1
       idsf=idsf(1);
       end
    else
        flag=1;
        nfr=nfr+1;
        mdbt=[mdbt;bm(i,4)];%multi-match mass of the minor subs
        d=sqrt(sum((subpos(bm(ifrom,1)+1,:)-subpos(bm(i,1)+1,:)).^2,2));
        idsf2=bms(bms(:,2)==bm(i,1),1);
        if numel(idsf2)>1
            idsf2=idsf2(1);
        end
        if isempty(idsf)
            d1=-100;
        else
            d1=sqrt(sum((subpos(bm(ifrom,1)+1,:)-sfpos(idsf+1,:)).^2,2))/SofteningHalo;
        end
        if isempty(idsf2)
            idsf2=bms(bms(:,2)==bm(ifrom,1),1);
            if numel(idsf2)>0
                idbt2=repmat(bm(i,1),size(idsf2));
                d2=sqrt(sum((subpos(idbt2+1,:)-sfpos(idsf2+1,:)).^2,2))/SofteningHalo;
                d2=-1*min(d2);
            else    
            d2=-100;
            end
        else
            d2=sqrt(sum((subpos(bm(i,1)+1,:)-sfpos(idsf2+1,:)).^2,2))/SofteningHalo;
        end
%         if d<2*SofteningHalo
          pair=[pair;bm(ifrom,1),bm(i,1),bm(ifrom,4),bm(i,4),d/SofteningHalo,d1,d2];
%         end
    end
end
nfr=nfr+nto;
ndup_btfr=nfr;
ndup_btto=nto;