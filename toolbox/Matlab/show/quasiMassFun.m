clear;
runnum='6702DM';
datadir=['/mnt/A4700/data/',runnum,'/subcat/anal/'];
data=load([datadir,'spmass']);
sdat=data(data(:,1)>=data(:,2),:);
x=logspace(1,2.7,10);
[y,bin]=histc(sdat(:,2),x);
xm=zeros(9,1);
for i=1:9
xm(i)=mean(sdat(bin==i,2));
end
figure;loglog(xm,y(1:9),'. -');
hold on;plot(xm,4.8e5*xm.^-2.5,'r-');

