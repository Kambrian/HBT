function rb=gas_binding_radius(rb_init,csat,rho_host)

mass=@(x) log(1+x)-x./(1+x);
vc2=@(x) (log(1+x)-x./(1+x))./x;
rhox=@(x) x./(1+x).^2;
rhoc=@(c) c./mass(c);

if rb_init<=0
    rb=0;
    return
end

rb0=rb_init;
% rb0=(mdm(1)/mhost(i))^(1/3)*r(i)/Rvir(i)*sqrt(G*mhost(i)*pmass/Rvir(i))/v(i);
% rho_host=rhox(r(i)/Rvir(i)*chost(i))*rhoc(chost(i));
rb=rb_init*sqrt(vc2(rb0*csat)/vc2(csat))*sqrt(rhox(rb0*csat)*rhoc(csat)/rho_host);

iter=0;          
while abs(rb0-rb)/rb>0.01
    rb0=rb;
    rb=rb_init*sqrt(vc2(rb0*csat)/vc2(csat))*sqrt(rhox(rb0*csat)*rhoc(csat)/rho_host);
    iter=iter+1;
    if iter>100
        warning('100 iterations exceeded in gas_binding_radius()');
        break;
    end
    if rb<0.001  %almost the center, probably no solution
        break;
    end
end