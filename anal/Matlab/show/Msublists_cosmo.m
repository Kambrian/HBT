function Mlist=Msublists_cosmo(massdata,cid,rmin,rmax)
% to list subhalo mass in radial range rmin<r<rmax
% rmin,rmax: comoving radius
% cid: central subhalo's index (subid+1)
% massdata=[Msub,CoM];

pos=massdata(:,2:4);
cpos=pos(cid,:);
np=size(pos,1);
lpos=repmat(cpos-rmax,np,1);
upos=repmat(cpos+rmax,np,1);
idlist=find(all(pos>lpos,2)&all(pos<upos,2));
idlist(idlist==cid)=[];%to make sure central is excluded
if isempty(idlist)
    Mlist=[];
    return;
end
r=sqrt(sum((pos(idlist,:)-repmat(cpos,size(idlist))).^2,2));
Mlist=massdata(idlist(r>rmin&r<rmax),1);
