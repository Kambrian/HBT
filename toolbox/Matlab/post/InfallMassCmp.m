RunName='6702';
Nsnap=99;

basedir=['/mnt/A4700/data/',RunName,'/subcat/anal/'];
Minff=importdata([basedir,'steller/SnapInfall_first_',num2str(Nsnap,'%03d')],',',1);
Minfl=importdata([basedir,'steller/SnapInfall_',num2str(Nsnap,'%03d')],',',1);
Lf=Minff.data(:,1:2);
Ll=Minfl.data(:,1:2);
r=Lf./Ll;
rs=r(:,1);rs=rs(~isnan(rs)&rs<inf&rs>0&rs~=1);
rg=r(:,2);rg=rg(~isnan(rg)&rg<inf&rg>0&rg~=1);
figure;hist(log(rs),-5.5:5.5)
figure;hist(log(rg),-5.5:5.5)