function frac=gas_fraction(M,c,z,req)
% req=Req/Rs, where Rho_gas=Rho_DM. in the range [0,c].
Rvir=virial_radius(M,z)
Rs=Rvir/c
Rhos=M/(4*pi*Rs^3*NFWmass(c))
Rho_gas0=Rhos/req/(1+req)^2/gasdens(req,c)
Mgas=4*pi*Rs^3*Rho_gas0*quadl(@(x) (gasdens(x,c).*x.^2),0,c);
frac=Mgas/M;



function rho=gasdens(x,c)
rho=exp(3*c/(NFWmass(c))*(log(1+x)./x-1));


function y=NFWmass(x)
y=log(1+x)-x./(1+x);

function R=virial_radius(M,z)
H0=0.1;
G=43007.1;
Omega0=0.3;OmegaL=0.7;
pmass=0.008848;
a=1/(1+z);
Hratio=sqrt(Omega0 /a^3+ (1 - Omega0-OmegaL) /a^2+OmegaL);
Hz=H0*Hratio;
OmegaZ=Omega0/a^3/Hratio^2;
x=OmegaZ-1;
virialF=18.0*pi^2+82.0*x-39.0*x^2; %<Rho_vir>/Rho_cri
R=(2.0*G*M*pmass/virialF/Hz^2)^(1.0/3)/a;

