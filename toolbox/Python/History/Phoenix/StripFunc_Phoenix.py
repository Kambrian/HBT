import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
from myutils import skeleton
from nfw import NFWHalo
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

PhHalo={}
for H in 'ABCDEFG':
  PhHalo[H]=HaloData('Ph'+H+'2')
PhHalo['a']=HaloData('A1')

Rref=1
nbin=30
Rconv=0.3
nMinInfall=0
iInfall=0
fMin=1e-5

Rstar=2.8
Mstar=0.06
alpha=0.99
beta=0.9
sigma=1.
fs=0.56 #1????

colors='rgbcmyk'

#ratio
xbin=np.logspace(-2, np.log10(2), 30)
denHaloAll=[]
denSubAll=[]
denRatAll=[]
for i,H in enumerate('ABCDEFG'):
  rHalo,denHalo,denHaloRef,denHaloErr=PhHalo[H].get_host_density(xbin, Rref)
  rSub,denSub,denSubRef,denSubErr=PhHalo[H].get_sub_density((PhHalo[H].m>fMin*PhHalo[H].Mvir/PhHalo[H].mP), xbin, Rref)
  denRatAll.append(denSub/denHalo)
denRatAll=np.array(denRatAll)
denRat=denRatAll.mean(0)
plt.figure()
plt.loglog(rSub, denRat, 'g')
plt.fill_between(rSub, denRat+denRatAll.std(0), denRat-denRatAll.std(0), color='g', alpha=0.3)
plt.yscale('log', nonposy='clip')
f=(rSub<1)&(denRat>0)
pars=powerlaw_fit(rSub[f], denRat[f])
plt.plot(rSub, (rSub/pars[0][0])**pars[0][1], 'k')
print pars
beta=pars[0][1]
x0=pars[0][0]

##plot profiles

## mass fraction
nbin=30
xbin=np.logspace(-1.5, np.log10(2), nbin)
nMinInfall=1000
x=[]
y=[]
for H in 'ABCDEFG':
  f=PhHalo[H].massTVV[:,iInfall]>nMinInfall
  x.extend(PhHalo[H].r[f]/PhHalo[H].Rvir)
  y.extend((1.*PhHalo[H].m/PhHalo[H].massTVV[:,iInfall])[f])
x=np.array(x)
y=np.array(y)
countAll,_=np.histogram(x, xbin)
plt.figure()
plt.plot(x,y, 'c.', alpha=0.1)
p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
plt.plot(xmid, p[1], 'r-', lw=2, label='Resolved')
plt.plot(xmid, p[[0,2]].T, 'r--',lw=2)
p,xmid=cross_section(x,y,xbin,44+0.56*np.array([(100-68.3)/2,50,(100+68.3)/2]))
plt.plot(xmid, p[1], 'k-', lw=2, label='Disruption-corrected')
plt.plot(xmid, p[[0,2]].T, 'k--', lw=2)
plt.yscale('log', nonposy='clip')
plt.xscale('log')
f0=(xmid<1.)&(xmid>0.01)&(p[1]>100./nMinInfall)
#suppression function
yscale=(p[1][f0]/(xmid[f0]/x0)**beta).mean()
plt.plot(rSub, denRat*yscale, 'bo',label=r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')

ypred=(xbin/x0)**beta*yscale
print 'R*:', x0*yscale**(-1./beta), 'beta:', beta
ysig=1.
plt.plot(xbin, ypred, 'g-', lw=2, label=r'$\beta=%.1f$'%(beta))
plt.plot(xbin, ypred*np.exp(ysig), 'g--', lw=2)
plt.plot(xbin, ypred/np.exp(ysig), 'g--', lw=2)

#tidal limits
Host=NFWHalo(6.57e4,5.96)
R=np.logspace(-2,0.3, 50)*Host.Rv
plt.loglog(R/Host.Rv, Host.strip_func(NFWHalo(1e-4*Host.M), R, k=1), 'm-', label='Tidal limit')
#Host=NFWHalo(183,15)
#R=np.logspace(-2,0.3, 50)*Host.Rv
#plt.loglog(R/Host.Rv, Host.strip_func(NFWHalo(1e-4*Host.M), R, k=1), 'm--', label='Aq Tidal limit')
plt.plot(xbin, (xbin/2.07)**1.43, 'k-', lw=4, alpha=0.2, label=r'$\beta=1.43$')
#plt.plot(xbin, (xbin/3.55)**0.9, 'k:')

plt.legend(loc=4,fontsize=15)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$m/m_0$')
plt.xlim([0.04,3])
plt.ylim([1e-3,1])
plt.plot(plt.xlim(), 20./nMinInfall*np.array([1,1]), 'r:')
plt.plot(plt.xlim(), 100./nMinInfall*np.array([1,1]), 'r:')
plt.text(plt.xlim()[0]*1.1, 20./nMinInfall, '20', color='r')
plt.text(plt.xlim()[0]*1.1, 100./nMinInfall, '100', color='r')
#plt.savefig(outdir+'/StrippingResolved.Phoenix.pdf')
