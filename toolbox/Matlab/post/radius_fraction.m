function rt=radius_fraction(ft,csat)

vc2=@(x) (log(1+x)-x./(1+x))./x;

if ft<=0
    rt=0;
    return
end

rt0=ft;
rt=ft*vc2(csat)/vc2(rt0*csat);
while abs(rt0-rt)/rt>0.01
    rt0=rt;
    rt=ft*vc2(csat)/vc2(rt0*csat);
end