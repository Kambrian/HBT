import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'
  
iInfall=0

AqA1=HaloData('AqA1')
PhA2=HaloData('PhA2')

zAqA1=np.loadtxt(datadir+'/Aquarius/AqA1Redshift.txt')
zPhA2=np.loadtxt(datadir+'/Phoenix/PhA2Redshift.txt')

AqA1.z=zAqA1[AqA1.snapTVV[:,iInfall], 1]
PhA2.z=zPhA2[PhA2.snapTVV[:,iInfall], 1]


def z_dist(H):
  plt.figure()
  m=H.massTVV[:,iInfall]
  a=[]
  for fmin in [1e-6,1e-5,1e-4,1e-3]:
	f=(m>fmin*H.Mvir/H.mP)&(m<fmin*10*H.Mvir/H.mP)
	plt.hist(np.log10(1.+H.z[f]), 50, normed=True,histtype='step',label='%e'%fmin)
	a.append(np.log(1+H.z[f]).mean())
  plt.legend()
  return a

a=z_dist(AqA1)
b=z_dist(PhA2)
plt.figure()
plt.plot(a)
plt.plot(b)

print 0.44*np.log(1+AqA1.z).mean()
print 0.44*np.log(1+PhA2.z).mean()