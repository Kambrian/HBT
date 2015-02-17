m=load('/mnt/A4700/data/6702DM/subcat/anal/Masses.dat');
ms=load('/mnt/A4700/data/6702DM/subcatS/anal/Masses.dat');


%%
f=logical(ones(size(m,1),1));
f=f&(m(:,end)==0);
f=f&(m(:,end-1)==0&m(:,end-2)==0&m(:,end-3)==0);
m1=m(f,3);
m2=ms(f,3);
figure;
x=logspace(-2,5,10);
% x=15;
loghist(m1,x,'1','r-',m1);
hold on;
loghist(m2,x,'1','-',m2);
legend('HBT','SUBFIND');