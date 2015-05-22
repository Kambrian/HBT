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
x=f['/x'][...]
r=np.sqrt(np.sum(x**2,1))*h0
countM,tmp=np.histogram(r, np.hstack([0., xbin]))#dM
densityHalo=countM/vol*f['/PartMass'][0]*h0
data=np.array([xcen, densityHalo]).T
np.savetxt(outdir+'Aq'+H+'2Halo.dat', data, header='r[kpc/h], rho[1e10Msun/h/(kpc/h)^3]')
#data=np.loadtxt(outdir+'Ph'+H+'2Halo.dat')
#densityHalo=data[:,1].T

f2=h5py.File(datadir+'sublist.subfind.hdf5','r')
m=f2['/PartMass'][...]
xs=f2['/x'][...]
rs=np.sqrt(np.sum(xs**2,1))*h0

countM,tmp=np.histogram(rs[m>100], np.hstack([0., xbin]))#dM
densitySub=countM/vol
data=np.array([xcen, densitySub/densityHalo]).T
np.savetxt(outdir+'Aq'+H+'2rat100.dat', data)

countM,tmp=np.histogram(rs[m>1000], np.hstack([0., xbin]))#dM
densitySub=countM/vol
data=np.array([xcen, densitySub/densityHalo]).T
np.savetxt(outdir+'Aq'+H+'2rat1000.dat', data)

plt.loglog(xcen, densitySub/densityHalo, 'r')
fl=(densitySub>0)&(xcen>1)&(xcen<200)
p=np.polyfit(np.log(xcen)[fl], np.log(densitySub/densityHalo)[fl],1)
plt.plot(xcen, np.exp(np.polyval(p,np.log(xcen))), 'k')
print p