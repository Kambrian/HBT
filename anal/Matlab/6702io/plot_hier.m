function h=plot_hier(subcat_in,subid,dims)
global pos subcat rtidal
subcat=subcat_in;
if isempty(rtidal)
   rtidal=comoving_virial_radius(subcat.SubLen);
end
if nargin<3,dims='xy';end
if numel(dims)~=2,error('dims must be two letters from xyz');end

dims=lower(dims);
xdim=~isempty(strfind(dims,'x'));
ydim=~isempty(strfind(dims,'y'));
zdim=~isempty(strfind(dims,'z'));

pos=[];i=0;lab=cell(3,1);
if xdim
    pos=subcat.SubProp.CoM(:,1);
    i=i+1;lab(i)={'x'};
end
if ydim
    pos=[pos,subcat.SubProp.CoM(:,2)];
    i=i+1;lab(i)={'y'};
end
if zdim
    pos=[pos,subcat.SubProp.CoM(:,3)];
    i=i+1;lab(i)={'z'};
end

h=gca();
plot(pos(subid+1,1),pos(subid+1,2),'s');
hold on;
plotsub(subid+1);
xlabel(lab{1});ylabel(lab{2});title(['sub_hierarchy for subid=',num2str(subid)]);
hold off;


function plotsub(subind)
global rtidal pos subcat 
circleplt(pos(subind,1:2),rtidal(subind));
hold on;
sonind=subcat.sub_hierarchy.sub(subind)+1;
nextind=subcat.sub_hierarchy.next(subind)+1;
if sonind>0
    plot([pos(subind,1),pos(sonind,1)],[pos(subind,2),pos(sonind,2)],'r-o');
    plotsub(sonind);
end
if nextind>0
    plot([pos(subind,1),pos(nextind,1)],[pos(subind,2),pos(nextind,2)],'-');
    plotsub(nextind);
end



function circleplt(x,r,linetype)
if nargin<3, linetype='-';end
if numel(x)~=2, error('usage: circleplt([x,y],r,linetype)');end
t=0:0.01:2*pi;
plot(x(1)+r*cos(t),x(2)+r*sin(t),linetype);

