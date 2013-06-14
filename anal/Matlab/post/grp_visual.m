clear;
runnum='6702DM';
grpname='S51G86';snaprange=51:99;
datadir=['/mnt/A4700/data/',runnum,'/subcat/anal/follow_group/',grpname];
scaleF_file=['/mnt/A4700/data/',runnum,'/subcat/Redshift.dat'];
tmp=load(scaleF_file);a=tmp(:,2);  

xmin=[];
xmax=[];
n=numel(snaprange);
m=cell(n,1);
com=cell(n,1);
cen=cell(n,1);
for i=1:n
    snap=snaprange(i);
    data=load([datadir,'/',num2str(snap,'%d')]);
    m{i}=data(:,1:2);
    com{i}=data(:,3:5);
    cen{i}=data(:,6:8);
xmin=[xmin;min(cen{i});];
xmax=[xmax;max(cen{i});];
end
xmin=min(xmin);
xmax=max(xmax);
%
figure('position',[1,1,1200,900]);
h1=subplot(1,2,1);
set(gca,'nextplot','replacechildren');
h2=subplot(1,2,2);
set(gca,'nextplot','replacechildren');
data=load([datadir,'/host']);
Mh=data(:,2:3);
xh=data(:,4:6);

nn=size(m{1},1);
for i=1:n
    snap=snaprange(i);
    axes(h1);
    for j=1:nn
        if m{i}(j,1)>0
            circle(com{i}(j,1),com{i}(j,2),m{i}(j,2),'b-');
%         plot(com{i}(j,1),com{i}(j,2),'o','Markersize',m{i}(j,2));
        hold on;
        end
    end
    plot(xh(i,1),xh(i,2),'rp');
    circle(xh(i,1),xh(i,2),Mh(i,2),'r-');
    hold off;
    axis([xmin(1),xmax(1),xmin(2),xmax(2)]);
    text(xmax(1)-(xmax(1)-xmin(1))/4,xmax(2)-(xmax(2)-xmin(2))/10,['a=',num2str(a(snap+1),'%1.2g'),',    z=',num2str(1/a(snap+1)-1,'%1.2g')]);
    title(['(',grpname,',Sim6702DM),orbits of members']);
    
    axes(h2);
    for j=1:nn
    circle(cen{i}(j,1),cen{i}(j,2),m{1}(j,2),'b-');
%     plot(cen{i}(j,1),cen{i}(j,2),'o','Markersize',m{1}(j,2));
    hold on;
    end
    plot(xh(i,1),xh(i,2),'rp');
    circle(xh(i,1),xh(i,2),Mh(i,2),'r-');
    hold off;
    axis([xmin(1),xmax(1),xmin(2),xmax(2)]);
    set(gca,'xtick',[],'ytick',[]);
    title('Most Bound Position');
       
    fr(i)=getframe(h1,[-40,-20,1000,800]);
%     pause(10);
end
figure('position',[1,1,1200,900]);
movie(fr,10,1);
outputdir='/home/kam/BT/data/trace';
movie2avi(fr,[outputdir,'/grp_orb_',grpname,'.avi'],'fps',10);