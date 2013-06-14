function [xmass,mfunspec,mfunspecln,mfuncum]=mass_sumcount(Mlist,nbin,xmin,xmax)
% divide Mlist into nbin and count dN
% xmass:size(nbin,1), center of mass bins
% mfunspec: size(nbin,2), [dN/dM,error]
% mfunspecln: same as above for dN/dlnM
% mfuncum: size(nbin,2), cumulative mass counts [N(>M),error]

nbin=nbin+1;
if nargin<3
xmin=min(Mlist);xmax=max(Mlist)*1.001;
end
x=logspace(log10(xmin),log10(xmax),nbin)';dlnx=log(x(2))-log(x(1));
[massfun,bin]=histc(Mlist,x);

for i=1:nbin
    massfun(i)=sum(Mlist(bin==i));  % mass (fraction) weighted massfunction
end

summass=cumsum(massfun(nbin:-1:1));

massfun=massfun(1:nbin-1);
errmass=sqrt(massfun);

xmass=zeros(nbin-1,2);
xmass(:,1)=x(1:nbin-1);
for i=1:nbin-1
    xmass(i,2)=mean(Mlist(bin==i));
    if isnan(xmass(i,2))
        xmass(i,2)=(x(i)+x(i+1))/2;
    end
end

massfun=[massfun,errmass];
mfunspec=massfun./repmat(diff(x),1,2);
mfunspecln=massfun/dlnx;
mfuncum=summass(nbin:-1:2);
mfuncum=[mfuncum,sqrt(mfuncum)];