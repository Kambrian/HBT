function t=cosmo_time(a_start,a_end,OmegaM0,OmegaL0,H0)
% if the input H0 is in units of 100km/s/Mpc, then t in units of
% 1/(100km/s/Mpc)=3.09e17s=9.8Gyr
OmegaK0=1-OmegaM0-OmegaL0;
t=quad(@(a) 1./sqrt(OmegaM0./a+OmegaL0*a.^2+OmegaK0)/H0, a_start,a_end);
