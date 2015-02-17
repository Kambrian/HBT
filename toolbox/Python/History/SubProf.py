import numpy as np
import matplotlib.pyplot as plt
import h5py

datadir='/gpfs/data/jvbq85/HBT/data/AqA4/subcat/anal/'

mbdfile=h5py.File(datadir+'MbdInfall.hdf5','r')
Mbd.m=mbdfile['/mass'][...]
Mbd.x=mbdfile['/x'][...]
Mbd.flag=mbdfile['/DirectInfall'][...]
Mbd.snapTVV=mbdfile['/snapTVV'][...]
Mbd.massTVV=mbdfile['/massTVV'][...]
Mbd.r=np.sqrt(np.sum(Mbd.x**2,1))
mbdfile.close()

f=(Mbd.m>100)&(Mbd.flag>0)
#f=Mbd.mTVV[:,0]>10000
R1=
rbin=np.logspace(-2,0.5,20)*


[xr3,n3,rhon3]=loghist(Mbd.r(f),rbin,[],[],1,3);
plot(xr3/R1, rhon3/(sum(n3)/R1^3),'r--','displayname','Ninfall>1e4');hold on;
% f=Mbd.m>100&Mbd.flag==0;
% [xr3,n3,rhon3]=loghist(Mbd.r(f),rbin,[],[],1,3);
% plot(xr3/R1, rhon3/(sum(n3)/R1^3),'r:','displayname','Ninfall>1e4');hold on;
f=Mbd.mTV(:,1)>1000&Mbd.flag>0;
[xr3,n3,rhon3]=loghist(Mbd.r(f),rbin,[],[],1,3);
plot(xr3/R1, rhon3/(sum(n3)/R1^3),'b:','displayname','Ninfall>1e3');
legend('show')
xscale('log');
yscale('log');
% print('-depsc','/work/Projects/SubProf/plots/A4subprof_unstripped.eps')
%%
f=Mbd.flag>0;
figure();
plot(Mbd.x(f,1), Mbd.x(f,2),'r.')
hold on;
f=Mbd.flag==0;
plot(Mbd.x(f,1), Mbd.x(f,2),'g.')
f=Mbd.flag<0&Mbd.m>0&Mbd.m<4.8e6;
plot(Mbd.x(f,1), Mbd.x(f,2),'k.')