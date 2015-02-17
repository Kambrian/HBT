function [n,en,mfun,err,gof]=fit_massfun_ratio(data,data0,Mhost,flag_plot)
% to fit for the ratio between mfun and mfun0
% return value:
%   n: ratio
%   ne: error of n
%   mfun: mfun./mfun0
%   err: error for mfun
%   gof: chi2/dof, goodness of fit


% mfun0=data0.mfunspecln/Mhost;
% mfun=data.mfunspecln/Mhost;

mfun0=data0.mfunspec(:,1).*data0.xmass(:,2)/Mhost; %subhalo mass per massbin
mfun0=[mfun0,data0.mfunspecln(:,2)/Mhost]; % adopt poisson error from mfunspecln
mfun=data.mfunspec(:,1).*data.xmass(:,2)/Mhost;
mfun=[mfun,data.mfunspecln(:,2)/Mhost];

err0=(mfun0(:,2)./mfun0(:,1)).^2;  % square of relative error in mass function
err=(mfun(:,2)./mfun(:,1)).^2;

mfun=mfun(:,1)./mfun0(:,1);

err=sqrt(err+err0).*mfun;   %absolute error for the normalized massfun
% err=ones(size(err));
flt=mfun>0;  %exclude null points
% flt(1)=0;  %exclude the first point
n=sum(mfun(flt)./err(flt).^2)./sum(err(flt).^-2);
en=sqrt(1./sum(err(flt).^-2));  % covariance matrix
chi2=sum((mfun(flt)-n).^2./err(flt).^2);
dof=numel(flt(flt))-1;  %degree of freedom
gof=chi2/dof; 
disp('chi2/dof=');disp(gof);
if gof>1, en=en*gof; end 

if nargin<4, flag_plot=0; end

if flag_plot

    xmass=data.xmass(flt,2);
%     errorbar(log10(xmass)+10,mfun(flt),err(flt),'.-');
    plot(log10(xmass)+10,mfun(flt),'.:');
    hold on;
%     errorbar(log10([xmass(1),xmass(end)])+10,[n,n],[en,en],'r--');
    plot(log10([xmass(1),xmass(end)])+10,[n,n],'r-');
end