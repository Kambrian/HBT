function rb=dm_binding_radius(rb_init,csat)

vc2=@(x) (log(1+x)-x./(1+x))./x;

if rb_init<=0
    rb=0;
    return
end

rb0=rb_init;
% rb0=par_dm(2)*r(i)/Rvir(i)/sqrt(vc2(r(i)/Rvir(i)*chost(i))/vc2(chost(i)))/sqrt(3-mgrad(r(i)/Rvir(i)*chost(i)));
rb=rb_init*sqrt(vc2(rb0*csat)/vc2(csat));

iter=0;          
while abs(rb0-rb)/rb>0.01
    rb0=rb;
    rb=rb_init*sqrt(vc2(rb0*csat)/vc2(csat));
    iter=iter+1;
    if iter>100
        warning('100 iterations exceeded in gas_binding_radius()');
        break;
    end
    if rb<0.001  %almost the center, probably no solution
        break;
    end
end