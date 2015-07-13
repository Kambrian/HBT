''' generate Monte-Carlo samples of subhaloes '''
import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from nfw import NFWHalo,mean_concentration
from myutils import myhist
plt.ion()
#rootdir='/gpfs/data/jvbq85/SubProf/'
rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'/plots/'

def BoostFactor(Mhalo, model='Gao', Mlim=1e-16):
  '''Mhalo and Mlim in 1e10 Msun/h'''
  if model=='Gao':
	return 1.6e-3*(Mhalo*1e10/0.73)**0.39*(Mlim/1e-16)**-0.226
  elif model=='Pinzke':
	return 0.76*0.023*(Mhalo/Mlim)**0.226
  else:
	raise "unknown model"
  
def getBoost(ifile):
  f=h5py.File(datadir+'MockHalo-%d.hdf5'%ifile)
  count,m,mAcc,R,Rp,phi,mu=f['Common']
  MHost=f['MHost'][...]
  cAcc,rt,Lt,Lv=f['Ludlow']
  b1=np.sum((count*Lt)[(m>1e-16)&(R<1.)])/NFWHalo(MHost).Ls
  cAcc,rt,Lt,Lv=f['Maccio']
  b2=np.sum((count*Lt)[(m>1e-16)&(R<1.)])/NFWHalo(MHost).Ls
  return MHost,b1,b2

#data=[getBoost(i) for i in xrange(20)]
#data=np.array(data).T
#np.savetxt(datadir+'/HaloBoost.dat', data)
data=np.loadtxt(datadir+'/HaloBoost.dat')

plt.figure()
plt.loglog(data[0]*1e10,data[1],'r-', label=r'Ludlow14 $M(c)$')
plt.plot(data[0]*1e10,data[2],'g--', label=r'Maccio08 $M(c)$')

plt.plot(data[0]*1e10, BoostFactor(data[0], 'Gao'), 'k', alpha=0.3, lw=3, label='Gao12')
plt.plot(data[0]*1e10, BoostFactor(data[0], 'Pinzke'), 'k:', label='Pinzke11')

plt.xlabel(r'$M_{200}[M_\odot/h]$')
plt.ylabel('Boost Factor')
x=np.logspace(-4,3)
plt.legend(loc=2,fontsize=15)
plt.savefig(outdir+'/HaloBoost.pdf')

#TODO: solve the inconsistency?????? state that the model likelihood overestimate the boost factor by a factor of ...