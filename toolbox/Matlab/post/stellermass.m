%%
% partmass=0.000058378;
partmass=0.008848;
h=0.7;
m0=31.6;
a=-0.39;b=-1.96;
k=10^10.35;
Minf=SnapInfall(:,2)*partmass;
Msat=2*k./((Minf/m0).^a+(Minf/m0).^b); %in Units of M_sun/h
Ksat=10^0.42/h^1.033*Msat.^0.967;
