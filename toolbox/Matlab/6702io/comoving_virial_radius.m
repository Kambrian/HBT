function rvir=comoving_virial_radius(mass)
global header HUBBLE0 G
    Hratio=header.Hz/HUBBLE0;
    scaleF=header.time;
    OmegaZ=header.Omega0/scaleF^3/Hratio^2;
  	virialF=18.0*3.1416*3.1416+82.0*(1-OmegaZ)-39.0*(1-OmegaZ)^2;
    rvir=(mass*2*G*header.mass(2)/virialF/header.Hz^2).^(1.0/3)/scaleF;
