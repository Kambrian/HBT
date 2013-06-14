function rmin=radial_range(hist_bin,pmass_bin,Snappar)
Ninfall_Min=10;Minfall_Min=1000;

Nhist=size(hist_bin,1);
rmin=zeros(Nhist,1)-1;

for h=1:Nhist
    [Ninfall,t,r,v,vt2,kt,mdm,mhost,chost,Hz,Rvir,virialF,virsubs]=get_DMstrp_history(hist_bin{h},pmass_bin(h),Snappar(h,2),1,1,1);
    if Ninfall>Ninfall_Min&&mdm(1)>Minfall_Min&&all(chost>0)
        rmin(h)=min(r./Rvir);
    end
end
rmin=rmin(rmin>0);