import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'


##plot profiles
iInfall=0
#H=HaloData(sys.argv[1])


def get_frac(H, rmin=0.5, rmax=0.8):
  yfrac=[]
  for fMin in x:
	f=(H.massTVV[:,iInfall]>fMin*H.Mvir/H.mP)&(H.r>rmin*H.Rvir)&(H.r<rmax*H.Rvir) #&(H.massTVV[:,iInfall]<fMin*x[2]/x[1]*H.Mvir/H.mP)
	#f=(H.flag>fMin*H.Mvir/H.mP)&(H.r>0.5*H.Rvir)&(H.r<0.8*H.Rvir)&(H.flag<fMin*x[2]/x[1]*H.Mvir/H.mP)
	print np.sum(f)
	yfrac.append(np.sum(f&(H.m>0))*1./np.sum(f))
  return yfrac

colors='krgbc'
plt.figure()
x=np.logspace(-7,-3,30)
for i in range(1,2):
  frac=get_frac(HaloData('AqA%d'%i))
  plt.semilogx(x,frac, label='AqA%d'%i, color=colors[i-1])
  #frac=get_frac(HaloData('AqA%d'%i), rmin=0, rmax=1.)
  #plt.semilogx(x,frac, '--', label='AqA%d'%i, color=colors[i-1])
plt.xlabel(r'$m_{\rm acc}/M_{200}$')
plt.ylabel(r'$f_{\rm s}(>m_{\rm acc}/M_{200})$')
#plt.plot(plt.xlim(), [0.56,0.56], 'k:')
plt.ylim([0.3,0.62])
#H=HaloData('AqA1')
#plt.plot(1000.*H.mP/H.Mvir*np.array([1.,1.]), plt.ylim(), ':')
#plt.legend(loc=3, fontsize=10)
#plt.savefig(outdir+'/SurvivalRate.A1.eps')

#frac=get_frac(HaloData('AqA2',None))
#plt.semilogx(x,frac, label='HBTA2')

#frac=get_frac(HaloData('AqA3',None))
#plt.semilogx(x,frac, label='HBTA3')
plt.figure()
x=np.logspace(-6,-3,20)
y=[]
for h in 'ABCDEF':
  frac=get_frac(HaloData('Aq'+h+'2'))
  y.append(frac)
  #plt.semilogx(x,frac,':')
y=np.array(y)
plt.semilogx(x, np.nanmean(y,0), 'r-', label='Aquarius')
plt.fill_between(x, np.nanmean(y,0)+np.nanstd(y,0), np.nanmean(y,0)-np.nanstd(y,0), color='r', alpha=0.3)

y=[]
for h in 'ABCDEF':
  frac=get_frac(HaloData('Ph'+h+'2'))
  y.append(frac)
  #plt.semilogx(x,frac,':')
y=np.array(y)
plt.semilogx(x, np.nanmean(y,0), 'g-', label='Phoenix')
plt.fill_between(x, np.nanmean(y,0)+np.nanstd(y,0), np.nanmean(y,0)-np.nanstd(y,0), color='g', alpha=0.3)

plt.xlabel(r'$m_{\rm acc}/M_{200}$')
plt.ylabel('Survival Rate')
  
plt.legend(loc=4)
#plt.plot(plt.xlim(), [0.55,0.55], 'k:')
#plt.savefig(outdir+'/SurvivalRate.pdf')