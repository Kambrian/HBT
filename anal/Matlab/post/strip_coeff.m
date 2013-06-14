function par=strip_coeff(Chost,Csat)

mass=@(x) log(1+x)-x./(1+x);
mgrad=@(x) x.^2./(1+x).^2./mass(x); %dlnM/dlnx, x=r/rs
vc=@(x) sqrt(mass(x)./x);
rhox=@(x) x./(1+x).^2;
rhoc=@(c) c./mass(c); %F=rhoc*rhox; F=rho*r^2/rv^2;
fdmh=@(x) vc(x).*sqrt(3-mgrad(x)); %dependence of tidal radius on host profile
% fgass=@(x) 1./vc(x)./sqrt(rhox(x));  %thermal timescale
fgass=@(x) 1./sqrt(rhox(x));    %dynamical timescale
fgash=@(x) sqrt(rhox(x));


Adm=quad(vc,0.01*Csat,Csat)/vc(Csat)/Csat/sqrt(2);
Bdms=quad(mgrad,0.01*Csat,Csat)/Csat;
Bdmh=quad(fdmh,0.01*Chost,Chost)/vc(Chost)/Chost;
% Bgass=Bdms*quad(fgass,0.01*Csat,Csat)*vc(Csat)/sqrt(rhoc(Csat))/Csat;     %thermal timescale
Bgass=Bdms*quad(fgass,0.01*Csat,Csat)/sqrt(rhoc(Csat))/Csat;
Bgash=quad(fgash,0.01*Chost,Chost)*sqrt(rhoc(Chost))/Chost;
% Agas=sqrt(5/12);   %thermal timescale
Agas=Adm;

Bdm=Bdms*Bdmh/sqrt(2);
% Bgas=Bgass*Bgash/sqrt(6/5); %thermal timescale
Bgas=Bgass*Bgash; %dynamical timescale

% betag=2;alphag=1;
% betad=1;alphad=0.45;
betag=1;alphag=1;
betad=1;alphad=1;

par=[Adm/betad,Bdm/alphad/betad,Agas/betag,Bgas/alphag/betag];

