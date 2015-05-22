import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

A1=HaloData('A1')
A2=HaloData('A2')
PhA2=HaloData('PhA2')  
PhB2=HaloData('PhB2')  
PhC2=HaloData('PhC2')
PhD2=HaloData('PhD2')
PhE2=HaloData('PhE2')
PhF2=HaloData('PhF2')
PhG2=HaloData('PhG2')

##plot profiles
Rref=1
iInfall=0
fMin=1e-4

## resolved
nbin=20
xbin=np.logspace(-2, np.log10(2), nbin)
colors=['r','g','b']

def scaled_density(H, fMin, xbin, Rref):
  mMin=fMin*H.Mvir
  rHalo,denHalo,denHaloRef,denHaloErr=H.get_host_density(xbin, Rref)
  rSub,denSubAll,denSubRef0,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP), xbin, Rref)
  f=(denSubAll/denHalo>0)&(rSub<1.8)&(rSub>0.6)
  C=(denSubAll/denHalo)[f].mean()
  denSubRef0*=C
  denSubAll=denSubAll/denHalo/C
  rSub,denSub,denSubRef,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP)&(H.m>0), xbin, Rref)
  denSub/=denHalo*denSubRef0/denSubRef
  return rSub,denSub,denSubAll

plt.figure()
hall=[]
r,y,yall=scaled_density(A1, fMin, xbin, Rref)
plt.plot(r[yall>0], yall[yall>0], 'r-')
htmp,=plt.plot(r[y>0], y[y>0], 'r-', label='A1')
hall.append(htmp)
r,y,yall=scaled_density(A2, fMin, xbin, Rref)
plt.plot(r[yall>0], yall[yall>0], 'b-')
htmp,=plt.plot(r[y>0], y[y>0], 'b-', label='A2')
hall.append(htmp)
prof=[]
profall=[]
for H in [PhA2,PhB2, PhC2, PhD2, PhE2, PhF2, PhG2]:
  r,y,yall=scaled_density(H, fMin, xbin, Rref)
  prof.append(y)
  profall.append(yall)
prof=np.array(prof)
profall=np.array(profall)
y=prof.mean(0)
ystd=prof.std(0)
ylim=[prof.min(0), prof.max(0)]
yall=profall.mean(0)
yallstd=profall.std(0)
yalllim=[profall.min(0),profall.max(0)]
htmp,=plt.plot(r[y>0],y[y>0],'g-', label='Phoenix')
hall.append(htmp)
#plt.fill_between(r, ylim[0], ylim[1], color='g', alpha=0.3)
plt.fill_between(r[y>0], (y-ystd)[y>0], (y+ystd)[y>0], color='g', alpha=0.3)
plt.plot(r,yall,'g-')
#plt.fill_between(r, yalllim[0], yalllim[1], color='g', alpha=0.3)
plt.fill_between(r, yall+yallstd, yall-yallstd, color='g', alpha=0.3)

plt.ylim([5e-2,3])
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend(hall, [h.get_label() for h in hall], loc=4, fontsize=15)
plt.plot(plt.xlim(), [0.56,0.56], 'k:')
plt.plot(plt.xlim(), [1,1], 'k--')  
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
plt.savefig(outdir+'/SubprofInfall_resolved.rat.Phoenix.pdf')
