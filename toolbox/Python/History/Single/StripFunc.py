import sys,os
sys.path.append(os.path.abspath('..'))
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

simu=sys.argv[1] #input 'Aq' or 'Ph'
#simu='Aq'

if simu[:2]=='Aq':
  alpha=0.95
  sigma=1.1
  fs=0.54 #1????
else:
  alpha=0.94
  sigma=1.1
  fs=0.56

hlist=[simu]
Halo={}
for H in hlist:
  Halo[H]=HaloData(H)
Mhalo=np.mean([Halo[h].Mvir for h in hlist])
chalo=np.mean([Halo[h].Cvir for h in hlist])
print 'M=', Mhalo, 'c=', chalo

Rref=1
iInfall=0
fMin=1e-5

colors='rgbcmyk'

#ratio
nbin=30
xbin=np.logspace(-2, np.log10(2), nbin)
denHaloAll=[]
denSubAll=[]
denRatAll=[]
for i,H in enumerate(hlist):
  rHalo,denHalo,denHaloRef,denHaloErr=Halo[H].get_host_density(xbin, Rref)
  rSub,denSub,denSubRef,denSubErr=Halo[H].get_sub_density((Halo[H].m>fMin*Halo[H].Mvir/Halo[H].mP), xbin, Rref)
  denRatAll.append(denSub/denHalo)
denRatAll=np.array(denRatAll)
denRat=denRatAll.mean(0)
plt.figure()
plt.loglog(rSub, denRat, 'g')
plt.fill_between(rSub, denRat+denRatAll.std(0), denRat-denRatAll.std(0), color='g', alpha=0.3)
plt.yscale('log', nonposy='clip')
f=(rSub<1)&(denRat>0)&(rSub>0.1)
pars=powerlaw_fit(rSub[f], denRat[f])
plt.plot(rSub, (rSub/pars[0][0])**pars[0][1], 'k')
print pars
gamma=pars[0][1]
beta=gamma/alpha
x0=pars[0][0]
mustar=pars[0][2]

##plot profiles

## mass fraction
nbin=30
xbin=np.logspace(-1.5, np.log10(2), nbin)
nMinInfall=1000
x=[]
y=[]
for H in hlist:
  f=Halo[H].massTVV[:,iInfall]>nMinInfall
  x.extend(Halo[H].r[f]/Halo[H].Rvir)
  y.extend((1.*Halo[H].m/Halo[H].massTVV[:,iInfall])[f])
x=np.array(x)
y=np.array(y)
countAll,_=np.histogram(x, xbin)
plt.figure()
plt.plot(x,y, 'c.', alpha=0.1)
p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
plt.plot(xmid, p[1], 'r-', lw=2, label='Resolved')
plt.plot(xmid, p[[0,2]].T, 'r--',lw=2)
p,xmid=cross_section(x,y,xbin,100*(1.-fs)+fs*np.array([(100-68.3)/2,50,(100+68.3)/2]))
plt.plot(xmid, p[1], 'k-', lw=2, label='Disruption-corrected')
plt.plot(xmid, p[[0,2]].T, 'k--', lw=2)
plt.yscale('log', nonposy='clip')
plt.xscale('log')
f0=(xmid<1.)&(xmid>0.01)&(p[1]>100./nMinInfall)
#suppression function
yscale=(p[1][f0]/(xmid[f0]/x0)**beta).mean()
plt.plot(rSub, denRat**(1/alpha)*yscale, 'bo',label=r'$(\rho_{\rm Sub}/\rho_{\rm Halo})^{1/\alpha}$')

ypred=(xmid)**beta*mustar*yscale
print 'R*:', x0*yscale**(-1./beta), 'beta:', beta, 'mu*:', mustar*yscale
plt.plot(xmid, ypred, 'g-', lw=2, label=r'$\beta=%.1f$'%(beta))
plt.plot(xmid, ypred*np.exp(sigma), 'g--', lw=2)
plt.plot(xmid, ypred/np.exp(sigma), 'g--', lw=2)

#tidal limits
Host=NFWHalo(Mhalo,chalo)
R=np.logspace(-2,0.3, 50)*Host.Rv
plt.loglog(R/Host.Rv, Host.strip_func(NFWHalo(1e-4*Host.M), R, k=1), 'm-', label='Tidal limit')
#Host=NFWHalo(183,15)
#R=np.logspace(-2,0.3, 50)*Host.Rv
#plt.loglog(R/Host.Rv, Host.strip_func(NFWHalo(1e-4*Host.M), R, k=1), 'm--', label='Aq Tidal limit')
#plt.plot(xbin, (xbin/2.07)**1.43, 'k-', lw=4, alpha=0.2, label=r'$\beta=1.43$')
#plt.plot(xbin, (xbin/3.55)**0.9, 'k:')

plt.legend(loc=4,fontsize=15)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$m/m_{\rm acc}$')
plt.xlim([0.04,3])
plt.ylim([1e-3,1])
plt.plot(plt.xlim(), 20./nMinInfall*np.array([1,1]), 'r:')
plt.plot(plt.xlim(), 100./nMinInfall*np.array([1,1]), 'r:')
plt.text(plt.xlim()[0]*1.1, 20./nMinInfall, '20', color='r')
plt.text(plt.xlim()[0]*1.1, 100./nMinInfall, '100', color='r')
#plt.savefig(outdir+'/StrippingResolved.%s.pdf'%simu)

p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2, 99])
plt.plot(xmid, p[3], 'b-')
plt.figure()
f=(p[1]>100./nMinInfall)
plt.semilogx(xmid, np.log(p[2]/p[1]), 'r--')
plt.semilogx(xmid[f], np.log(p[2]/p[1])[f], 'r-')
plt.semilogx(xmid, np.log(ypred/p[0]), 'k-')
plt.semilogx(xmid, np.log(p[3]/ypred), 'b-')
p,xmid=cross_section(x,y,xbin,100*(1.-fs)+fs*np.array([(100-68.3)/2,50,(100+68.3)/2]))
f=(p[1]>100./nMinInfall)
plt.semilogx(xmid, np.log(p[2]/p[1]), 'g--')
plt.semilogx(xmid[f], np.log(p[2]/p[1])[f], 'g-')
plt.semilogx(xmid, np.log(p[2]/ypred), 'k-')
plt.xlim([0.1,2])
plt.ylim([0.5,1.5])