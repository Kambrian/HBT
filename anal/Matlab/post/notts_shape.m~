clear;

outdir='/home/kam/notts/plots';
codes={'AdaptaHOP','AHF','H3D','H6D','HBT','HSF','mendieta','Rockstar','STF','subfind','VOBOZ'};
cmap=colormap(hsv(numel(codes)));
lines={'-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.','-','--','-.'};
data=cell(numel(codes),1);
for i=1:numel(codes)
    file=['/home/kam/notts/shape/shape_raw_AqA4_',codes{i}];
    disp(['loading ',codes{i}]);
    tmp=importdata(file,',',1);
    data{i}=tmp.data;
end
%%
myfigure;
for i=1:numel(codes)
    r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100;
    Q=Q(filter,:);
    Q=Q(~logical(sum(isnan(Q),2)),:); %exclude nans;
    Q=Q(logical(sum(abs(Q),2)),:);%exclude Q=0;
    
    y=Q(:,3)./Q(:,1);
    [x,n]=linhist(y,0:0.06:1);
    plot(x,n,lines{i}, 'color',cmap(i,:),'displayname',codes{i}); hold on;
end
legend('location','northeastoutside');
text(0.2,190,'aspherical');
text(0.9,190,'spherical');
xlabel('$r_{min}/r_{max}$');
ylabel('N');
set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,4,22,15]);
print('-depsc',[outdir,'/axis_ratio.eps']);
%%
myfigure;
for i=1:numel(codes)
    r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100;
    Q=Q(filter,:);
    Q=Q(~logical(sum(isnan(Q),2)),:); %exclude nans;
    Q=Q(logical(sum(abs(Q),2)),:);%exclude Q=0;
    
    y=Q(:,3)./Q(:,2);
    [x,n]=linhist(y,0:0.06:1);
    plot(x,n,lines{i}, 'color',cmap(i,:),'displayname',codes{i}); hold on;
end
legend('location','northeastoutside');
xlabel('$r_{min}/r_{mid}$');
ylabel('N');
set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,4,22,15]);
print('-depsc',[outdir,'/axis_ratio2.eps']);
%%
myfigure;
for i=1:numel(codes)
    r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100;
    Q=Q(filter,:);
    Q=Q(~logical(sum(isnan(Q),2)),:); %exclude nans;
    Q=Q(logical(sum(abs(Q),2)),:);%exclude Q=0;
    
    y=Q(:,2)./Q(:,1);
    [x,n]=linhist(y,0:0.06:1);
    plot(x,n,lines{i}, 'color',cmap(i,:),'displayname',codes{i}); hold on;
end
legend('location','northeastoutside');
xlabel('$r_{mid}/r_{max}$');
ylabel('N');
set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,4,22,15]);
print('-depsc',[outdir,'/axis_ratio3.eps']);
%%
myfigure;
for i=1:numel(codes)
    r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100&~logical(sum(isnan(Q),2))&logical(sum(Q~=0,2));
    Q=Q(filter,:);
   
%     y=Q(:,3)./Q(:,1);
    y=std(Q,1,2)./mean(Q,2);
    [x,n]=linhist(y,0:0.04:5);
    plot(x,n,lines{i},'color',cmap(i,:),'displayname',codes{i}); hold on;
end
legend('location','northeastoutside');
text(0.02,180,'spherical');
text(1,180,'aspherical');
xlabel('$\sigma(r_i)/<r_i>$');
ylabel('N');
set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,4,22,15]);
print('-depsc',[outdir,'/axis_dispersion.eps']);
%%
myfigure;
for i=1:numel(codes)
    r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100&~logical(sum(isnan(Q),2))&logical(sum(Q~=0,2));
    Q=Q(filter,:);
    r=r(filter,:);
    m=m(filter,:);
    
    y=Q(:,3)./Q(:,1);
%     y=std(Q,1,2)./mean(Q,2);
    [xmed,ymed,ylim]=skeleton(log10(r),y,8,0.683);
    plot(xmed,ymed,lines{i},'color',cmap(i,:),'displayname',codes{i}); hold on;
end
legend('location','northeastoutside');
text(xmed(1),0.4,'aspherical');
text(xmed(1),0.75,'spherical');
xlabel('log($r[kpc/h]$)');
ylabel('Median($r_{min}/r_{max}$)');
set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,4,22,15]);
print('-depsc',[outdir,'/axis_ratio_profile.eps']);

%%
myfigure;
for i=1:numel(codes)
    r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100&~logical(sum(isnan(Q),2))&logical(sum(Q~=0,2));
    Q=Q(filter,:);
    r=r(filter,:);
    m=m(filter,:);
    
%     y=Q(:,3)./Q(:,1);
    y=std(Q,1,2)./mean(Q,2);
    [xmed,ymed,ylim]=skeleton(log10(r),y,8,0.683);
    plot(xmed,ymed,lines{i},'color',cmap(i,:),'displayname',codes{i}); hold on;
end
legend('location','northeastoutside');
text(xmed(1),0.1,'spherical');
text(xmed(1),0.4,'aspherical');
xlabel('log($r[kpc/h]$)');
ylabel('Median($\sigma(r_i)/<r_i>$)');
set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,4,22,15]);
print('-depsc',[outdir,'/axis_dispersion_profile.eps']);
%%
myfigure;
i=5;
 r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100&~logical(sum(isnan(Q),2))&logical(sum(Q~=0,2));
    Q=Q(filter,:);
    r=r(filter,:);
    m=m(filter,:);
    
    y=Q(:,3)./Q(:,1);
%     y=std(Q,1,2)./mean(Q,2);
    [xmed,ymed,ylim]=skeleton(log10(m),y,8,0.683);
    plot(xmed,ymed,lines{i},'color',cmap(i,:),'displayname',codes{i}); hold on;
fx=~isnan(xb);
h=area(xb(fx),[yb(fx,1),yb(fx,2)-yb(fx,1)]);
set(h(2),'facecolor',[0.9,0.9,0.9],'edgecolor','w');
set(h(1),'edgecolor','w','facecolor','w');
set(gca,'layer','top');    
for i=1:numel(codes)
    r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100&~logical(sum(isnan(Q),2))&logical(sum(Q~=0,2));
    Q=Q(filter,:);
    r=r(filter,:);
    m=m(filter,:);
    
    y=Q(:,3)./Q(:,1);
%     y=std(Q,1,2)./mean(Q,2);
    [xmed,ymed,ylim]=skeleton(log10(m),y,8,0.683);
    plot(xmed,ymed,lines{i},'color',cmap(i,:),'displayname',codes{i}); hold on;
    if(strcmp(codes{i},'HBT'))
        xb=xmed;yb=ylim;
    end
end
legend('location','northeastoutside');
% plot(xb,yb(:,1),'color',[0.95,0.95,0.95],'linewidth',4);plot(xb,yb(:,2),'color',[0.95,0.95,0.95],'linewidth',4);
text(xmed(1),0.45,'aspherical');
text(xmed(1),0.8,'spherical');
xlabel('log($m[particles]$)');
ylabel('Median($r_{min}/r_{max}$)');
set(gca,'ylim',[0.3,0.9]);
set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,4,22,15]);
print('-depsc',[outdir,'/axis_ratio_mass_profile.eps']);
%%
myfigure;
for i=1:numel(codes)
    r=data{i}(:,4);
    m=data{i}(:,5);
    Q=sqrt(data{i}(:,1:3));
    filter=r<500&m>100&~logical(sum(isnan(Q),2))&logical(sum(Q~=0,2));
    Q=Q(filter,:);
    r=r(filter,:);
    m=m(filter,:);
    
 %   y=Q(:,3)./Q(:,1);
    y=std(Q,1,2)./mean(Q,2);
    [xmed,ymed,ylim]=skeleton(log10(m),y,8,0.683);
    plot(xmed,ymed,lines{i},'color',cmap(i,:),'displayname',codes{i}); hold on;
end
legend('location','northeastoutside');
text(xmed(1),0.1,'spherical');
text(xmed(1),0.35,'aspherical');
xlabel('log($m[particles]$)');
ylabel('Median($\sigma(r_i)/<r_i>$)');
set(gcf,'papertype','A4');
set(gcf,'paperunits','centimeters','paperposition',[0.6,4,22,15]);
print('-depsc',[outdir,'/axis_dispersion_mass_profile.eps']);