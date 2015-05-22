'''compare HBT and subfind mass functions, and convergence with resolution'''

import sys,os
sys.path.append(os.path.abspath('..'))
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
Rref=1
nbin=20
iInfall=0
xbin=np.logspace(-7,-1,nbin)
def plotMF(H, xbin, iInfall=None, **kwarg):
  '''H: halodata
  iInfall: None for current mass, otherwise specify infall mass'''
  if iInfall is None:
	m=H.m
  else:
	m=H.massTVV[:,iInfall]
  print 'np:', H.Mvir/H.mP
  data=m[(H.r<Rref*H.Rvir)&(m>0)]*H.mP/H.Mvir
  n,_=np.histogram(data, xbin, weights=None)
  y=n/np.log(xbin[1]/xbin[0])
  x=xbin[:-1]*np.sqrt(xbin[1]/xbin[0])
  y=y*x
  #x*=H.Mvir;y/=H.Mvir;
  #plt.step(xbin[:-1]*H.Mvir, y/H.Mvir, where='post', label=H.name)
  l,=plt.loglog(x,y, **kwarg)
  f=n>4
  #f=f&(x<1e-4)
  #w=np.sqrt(n)[f]
  w=None
  polypar,polyV=np.polyfit(np.log(x)[f], np.log(y)[f],1, w=w , cov=True)
  #l2,=plt.plot(x, np.exp(np.polyval(polypar, np.log(x))), '--', color=l.get_color(), label=r'$\alpha=%.2f$'%(-polypar[0]))
  print polypar, np.exp(-polypar[1]/polypar[0])
  return l

plt.figure()
plotMF(HaloData('AqA1'), xbin, label='AqA1', color='k')
for i in '2345':
  H=HaloData('AqA'+i)
  l=plotMF(H,xbin, label=H.name)
  plotMF(HaloData('AqA'+i, None), xbin, color=l.get_color(), linestyle='--')

legend=plt.legend(loc=1,fontsize=12)

plt.xlabel(r'$m/M_{200}$')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}\ln m$')
plt.savefig(outdir+'CurrentMF.Aq.eps')

plt.figure()
plotMF(HaloData('AqA1'), xbin, iInfall=0, label=H.name)
for i in '2345':
  H=HaloData('AqA'+i)
  l=plotMF(H,xbin, iInfall=0, label=H.name)
  plotMF(HaloData('AqA'+i, None), xbin, iInfall=0, color=l.get_color(), linestyle='--')

legend=plt.legend(loc=1,fontsize=12)

plt.xlabel(r'$m/M_{200}$')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}\ln m$')
#plt.savefig(outdir+'CurrentMF.Aq.eps')