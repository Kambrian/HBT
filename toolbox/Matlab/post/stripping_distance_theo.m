x=logspace(-5,log10(0.2),15);
y=9/(2*(3^(1/3))^2)+x.^(1/3)/3+x.^(2/3)/(3*3^(1/3))-(7*x)/(9*(3^(1/3))^2);
y2=3^(1/3)+x.^(1/3)/3+(2*x.^(2/3))/(9*3^(1/3))+x/(9*(3^(1/3))^2);
figure;subplot(2,1,1);semilogx(x,y,'r-',x,y2,'b--');
set(gca,'yminortick','on');
xlabel('m/M');ylabel('D/Rvir(M)');legend('compression distance','stripping distance');
subplot(2,1,2);semilogx(x,y./y2);xlabel('m/M');ylabel('R_{strip}/R_{compress}');
