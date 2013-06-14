function [w,delta,D,D2,D3,w3]=collapse_barrier(a,OmegaM0)
% output:====
% w: collapse barrier for use with z=0 field, w=delta(z)/(D(z)/D(z=0))
% w2: same as w, but different fitting formula
% delta: collapse barrier
% D,D2,D3: growth factor from various solution/fitting formula
% input:====
% a: scale factor to evaluate at
% OmegaM0


x0=(1/OmegaM0-1)^(1/3);
x=a*x0;

delta=0.15*(12*pi)^(2/3)*(1-0.0123*log10(1+x.^3));  %Nakamura & Suto 97, eq.[C-28]

fun=@(xi,y) (1+xi^3*y.^(6/5)).^(-3/2);
D=zeros(size(x));
for i=1:numel(x)
D(i)=a(i)*sqrt(1+x(i)^3)*quad(@(y)fun(x(i),y),0,1);  %Peebles 84, eq.[14]   
end
D0=sqrt(1+x0^3)*quad(@(y)fun(x0,y),0,1); 
D=D/D0;  %normalize at z=0

w=delta./D;


OmegaMZ=1./(1+x.^3);
OmegaLZ=1-OmegaMZ;
D2=2.5*OmegaMZ.*a./(OmegaMZ.^(4/7)-OmegaLZ+(1-OmegaMZ/2).*(1+OmegaLZ/70));  %Cooray&Sheth, Halo Model Rev.,eq.[25]
D0=2.5*OmegaM0/(OmegaM0.^(4/7)-(1-OmegaM0)+(1-OmegaM0/2).*(1+(1-OmegaM0)/70));
D2=D2/D0;


D3=D;
for i=1:numel(x)
D3(i)=(1./OmegaMZ(i)-1).^(1/3)*hypergeom([1/3,1],11/6,1-1./OmegaMZ(i));  % Nakamura & Suto 1997, eq.[C-25]
end
D0=(1./OmegaM0-1).^(1/3)*hypergeom([1/3,1],11/6,1-1./OmegaM0);
D3=D3/D0;

w3=0.15*(12*pi)^(2/3)./a*hypergeom([1/3,1],11/6,1-1/OmegaM0).*OmegaMZ.^-0.215;  %Nakamura&Suto 97, eq.[C-31]


% D=D3,w=dlt/D=dlt/D3=w3;
% D2 is slightly different, and hence dlt/D2;