#absolute mass on x axis
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

def stackMF(halolist,xbin,nMin=100.,Rref=1,iInfall=0):
  '''process a list of haloes'''
  Nall=[]
  Mall=[]
  for h in halolist:
	H=HaloData(h)
	n=np.histogram(H.massTVV[:,iInfall][H.r<Rref*H.Rvir]*H.mP, xbin)[0]*1.
	n[xbin[:-1]<nMin*H.mP]=np.nan
	Nall.append(n)
	Mall.append(H.Mvir)
  Nall=np.array(Nall)
  Mall=np.array(Mall)
  y=Nall.sum(0)/Mall.sum(0)/np.log(xbin[1]/xbin[0])
  yall=Nall/Mall[:,None]/np.log(xbin[1]/xbin[0])
  x=xbin[:-1]*np.sqrt(xbin[1]/xbin[0])
  f=(Nall.sum(0)>0)
  par=powerlaw_fit(x[f], y[f])
  return {'p':par,'x':x,'y':y*x,'ya':yall*x}


def plotBand(x,data, label=None, color='r', ls='-', alpha=0.3, filt=None):
  y=np.mean(data,0)
  ystd=np.std(data,0)
  f=(y>0)&(ystd>0)
  if filt!=None:
	f=f&filt
  l,=plt.loglog(x[f], y[f], linestyle=ls, color=color, label=label)
  plt.fill_between(x[f], (y+ystd)[f], (y-ystd)[f], color=color, alpha=0.3)
  return l

xbin=np.logspace(-5,0,30)
Aq=stackMF(['Aq'+h+'2' for h in 'ABCDEF'], xbin)
xbin=np.logspace(-2,3,30)
Ph=stackMF(['Ph'+h+'2' for h in 'ABCDEF'], xbin)

print 'AqPars:', Aq['p']
print 'PhPars:', Ph['p']

## Mass Functions
plt.figure()
##Aq
x=Aq['x']
plotBand(x, Aq['ya'], 'Aquarius', 'r', '-')#, filt=(x>1e-6))
plt.plot(x, Aq['y'], 'r--')
plt.plot(x, x*Aq['p'][0][2]*x**Aq['p'][0][1], 'k-')
x=Ph['x']
plotBand(x, Ph['ya'], 'Phoenix', 'g', '-', 0.2)#, filt=(x>1e-6))  
plt.plot(x, Ph['y'], 'g--')
plt.plot(x, x*Ph['p'][0][2]*x**Ph['p'][0][1], 'k-')
plt.yscale('log', nonposy='clip')

legend=plt.legend(loc=1,fontsize=12)

plt.xlabel(r'$m[10^{10}M_{\odot}/h]$')
plt.ylabel(r'$(m/M_{200})\mathrm{d}N/\mathrm{d}\ln m$')
#plt.savefig(outdir+'MF.All.pdf')
