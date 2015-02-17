gamma=-1;
beta=1;
Gamma=@(c) 3*c/beta./(log(1+c)-c./(1+c));
figure;plot(2:10,Gamma(2:10),'.-');
vc2=@(x) (log(1+x)-x./(1+x))./x;
% rho=@(x) 1./x./(1+x).^2;
% rho=@(x) x.^gamma;
rho=@(x,c) (1+x).^(Gamma(c)./x);
% df=@(x) (5*x^2-2*log(1+x)-6*log(1+x)*x-4*log(1+x)*x^2+2*x)/x^3/(1+x)^4;
f=@(x,c) vc2(x).*rho(x,c);
x=logspace(-3,1,100);
colors=['r.-';'g.-';'k.-';'c.-'];colors=repmat(colors,20,1);
figure;
xd=[];yd=[];
for c=5:3:20
    y=f(x,c);[ym,i]=max(y);xm=x(i);xd=[xd,xm];yd=[yd,ym];
    loglog(x,y,colors(c,:),'displayname',num2str(c));hold on;
end
plot(xd,yd,'o-','displayname','peaks');
legend('show');
%%
syms x
f=(log(1+x)-x/(1+x))/x^2/(1+x)^2;
g=diff(f);