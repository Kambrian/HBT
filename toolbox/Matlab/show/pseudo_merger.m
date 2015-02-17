%% 假并和，而非分裂
clear;
outputdir='/home/kam/Documents/research/Galaxy/code/BoundTracing/data/show';
% SnapMin=25;SnapMax=45;
SnapMin=40;SnapMax=50;
SnapLoad=46;GrpLoad=218;
% SnapLoad=30;GrpLoad=6280;
% SnapLoad=35;GrpLoad=6508;
% SnapLoad2=56;GrpLoad2=2013;
n=SnapMax-SnapMin+1;
datadir='/mnt/A4700/data/8213/subcat/anal/image/';
Xmap=load([datadir,'posmap_S',num2str(SnapLoad),'G',num2str(GrpLoad,'%d'),'.',num2str(SnapMin)]);
Xmap=zeros(size(Xmap,1),3,n);
% Xmap2=load([datadir,'posmap_S',num2str(SnapLoad2),'G',num2str(GrpLoad2,'%d'),'.',num2str(SnapMin)]);
% Xmap2=zeros(size(Xmap2,1),3,n);
i=1;
for Nsnap=SnapMin:SnapMax
Xmap(:,:,i)=load([datadir,'posmap_S',num2str(SnapLoad),'G',num2str(GrpLoad,'%d'),'.',num2str(Nsnap)]);
% Xmap2(:,:,i)=load([datadir,'posmap_S',num2str(SnapLoad2),'G',num2str(GrpLoad2,'%d'),'.',num2str(Nsnap)]);
i=i+1;
end

xmin=min(min(Xmap,[],1),[],3);
xmax=max(max(Xmap,[],1),[],3);
figure();
set(gca,'nextplot','replacechildren');
Nsnap=SnapMin;
for i=1:n
%     figure;
% plot3(Xmap2(1:1:end,1,i),Xmap2(1:1:end,2,i),Xmap2(1:1:end,3,i),'b.','markersize',2);hold on;
plot3(Xmap(1:1:end,1,i),Xmap(1:1:end,2,i),Xmap(1:1:end,3,i),'r.','markersize',2);hold on;
% plot3(Xmap(1,1,i),Xmap(1,2,i),Xmap(1,3,i),'b+','markersize',6);hold on;
% plot3(Xmap(10,1,i),Xmap(10,2,i),Xmap(10,3,i),'g^','markersize',6);hold on;
% plot3(Xmap(50,1,i),Xmap(50,2,i),Xmap(50,3,i),'c*','markersize',6);hold on;
% plot3(Xmap(100,1,i),Xmap(100,2,i),Xmap(100,3,i),'k>','markersize',6);hold on;
% plot3(Xmap(200,1,i),Xmap(200,2,i),Xmap(200,3,i),'mo','markersize',6);hold on;
hold off;
axis([xmin(1),xmax(1),xmin(2),xmax(2),xmin(3),xmax(3)]);
title(num2str(Nsnap));
Nsnap=Nsnap+1;
fr(i)=getframe(gca,[-10,-10,600,800]);
end
figure();
movie(fr,10,1);