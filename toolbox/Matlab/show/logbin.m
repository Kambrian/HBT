function [xmean,ymean,xmed,ymed,xl]=logbin(x,y,nbin)

xi=logspace(log10(min(x)),log10(max(x)),nbin+1);
[n,bin]=histc(x,xi);
ymean=zeros(nbin,1);
ymed=ymean;
xmean=ymean;
xmed=ymean;
for i=1:nbin
    ymean(i)=mean(y(bin==i));
    ymed(i)=median(y(bin==i));
    xmean(i)=mean(x(bin==i));
    xmed(i)=median(x(bin==i));
end
xl=xi(1:end-1);