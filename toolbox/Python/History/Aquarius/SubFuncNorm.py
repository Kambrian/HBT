import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
h0=0.73

H=sys.argv[1]
datadir='/gpfs/data/jvbq85/HBT/data/Aq'+H+'2/subcat/anal/'
outdir='/gpfs/data/jvbq85/SubProf/data/'
nbin=50
xbin=np.logspace(np.log10(0.1), np.log10(500*h0), nbin)
xcen=xbin/np.sqrt(xbin[1]/xbin[0])
vol=np.diff(np.hstack([0., xbin])**3)*np.pi*4/3

f=h5py.File(datadir+'allpart.subfind.hdf5','r')
mP=f['/PartMass'][0]*h0
x=f['/x'][...]
r=np.sqrt(np.sum(x**2,1))*h0
countM,tmp=np.histogram(r, np.hstack([0., xbin]))#dM
cumM=countM.cumsum()*mP
data=np.array([xbin, cumM]).T
np.savetxt(outdir+'Aq'+H+'2HaloCum.dat', data, header='r[kpc/h], M(<r)[1e10Msun/h]')
#data=np.loadtxt(outdir+'Ph'+H+'2HaloCum.dat')
#densityHalo=data[:,1].T

f2=h5py.File(datadir+'sublist.subfind.hdf5','r')
m=f2['/PartMass'][...]
xs=f2['/x'][...]
rs=np.sqrt(np.sum(xs**2,1))*h0

countM,tmp=np.histogram(rs[(m>100)&(m<1e6)], np.hstack([0., xbin]))#dM
cumN=countM.cumsum()
data=np.array([xbin, cumN/cumM*0.95/((100*mP)**-0.95-(1e6*mP)**-0.95)]).T
np.savetxt(outdir+'Aq'+H+'2MFNorm.dat', data)

plt.loglog(data.T[0],data.T[1], 'r')