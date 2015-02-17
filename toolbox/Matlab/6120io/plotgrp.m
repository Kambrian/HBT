function h=plotgrp(cat,grpid,dims,vscale,autoscale)
%function h=plotgrp(cat,grpid,dims)
%  [cat] is of the type cat/srccat/subcat, plotting (sub)halo from cat/srccat/subcat, 
%  [dim] with dimensions by dims: 'xyz'(default if not set),'yz'...,'xyzv','xyv'...etc.
%  <grpid>      grpid starting from 0
%  [vscale] if vscale is present, plot velocity field with vel*vscale,0 to disable
%           vscale
%  [autoscale] if autoscale is present,plot velocity field with automatic scale to
%               prevent overlapping, then scale the vel vectors by autoscale, 
%                i.e.,quiver(x,y,.,vx,vy,.,autoscale), 0 to disable auto scale
%       if neither vscale nor autoscale is present, do not scale.
global Pdat

if nargin<3,dims='xyz';end
if ~strcmp(cat.property.particles,'matlab-index'),error('fresh_id2index before plotting');end

dims=lower(dims);
xdim=~isempty(strfind(dims,'x'));
ydim=~isempty(strfind(dims,'y'));
zdim=~isempty(strfind(dims,'z'));
vdim=~isempty(strfind(dims,'v'));

if nargin>3,vdim=1;end
if nargin<4,vscale=0;end
if nargin<5,autoscale=0;end

if strcmp(cat.property.particles,'id'), error('fresh particle id to index first');end
switch cat.property.type
    case 'cat'
        arr=cat.PIDorIndex(cat.Offset(grpid+1)+(1:cat.Len(grpid+1)));
    case {'subcat','srccat'}
        arr=cat.PSubArr{grpid+1};
    otherwise
        error('unknown catalogue type');
end
        
pos=[];vel=[];i=0;lab=cell(3,1);
if xdim
    pos=Pdat.Pos(arr,1);
    i=i+1;lab(i)={'x'};
    if vdim, vel=Pdat.Vel(arr,1);end
end
if ydim
    pos=[pos,Pdat.Pos(arr,2)];
    i=i+1;lab(i)={'y'};
    if vdim, vel=[vel,Pdat.Vel(arr,2)];end
end
if zdim
    pos=[pos,Pdat.Pos(arr,3)];
    i=i+1;lab(i)={'z'};
    if vdim, vel=[vel,Pdat.Vel(arr,2)];end
end

switch i
    case 1
        if vdim
            h=plot(pos,vel,'o');
            xlabel(['pos',upper(lab{1})]);
            ylabel(['vel',upper(lab{1})]);
        else
            h=plot(pos,'o');
            ylabel(['pos',upper(lab{1})]);
        end
    case 2
        if vdim
            if autoscale
            h=quiver(pos(:,1),pos(:,2),vel(:,1),vel(:,2),autoscale);
            else if vscale
            h=quiver(pos(:,1),pos(:,2),vel(:,1)*vscale,vel(:,2)*vscale);
                else
            h=quiver(pos(:,1),pos(:,2),vel(:,1),vel(:,2));        
                end
            end
        else
            h=plot(pos(:,1),pos(:,2),'o');
        end
        xlabel(lab{1});ylabel(lab{2});
    case 3
        if vdim
            if autoscale
            h=quiver3(pos(:,1),pos(:,2),pos(:,3),vel(:,1),vel(:,2),vel(:,3),autoscale);
            else if vscale
            h=quiver3(pos(:,1),pos(:,2),pos(:,3),vel(:,1)*vscale,vel(:,2)*vscale,vel(:,3)*vscale);
                else
            h=quiver3(pos(:,1),pos(:,2),pos(:,3),vel(:,1),vel(:,2),vel(:,3));
                end
            end
        else
        h=plot3(pos(:,1),pos(:,2),pos(:,3),'o');
        end
        xlabel('x');ylabel('y');zlabel('z');
    otherwise
        error('specify dimensions to plot: ''xyz'',"xy"... etc')
end
        title(['(sub)halo #',num2str(grpid),' from ',cat.property.type,num2str(cat.property.snap),'plotted at snap',num2str(Pdat.snap)]);
        
