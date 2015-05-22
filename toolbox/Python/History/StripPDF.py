import numpy as np
import matplotlib.pyplot as plt
from matplotlib.finance import candlestick
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
import scipy.stats as st
from MbdIO import *
from myutils import skeleton,density_of_points,percentile_contour
from nfw import NFWHalo
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

A1=HaloData('AqA1')  
A2=HaloData('AqA2')
A3=HaloData('AqA3')
A4=HaloData('AqA4')
A5=HaloData('AqA5')
fmt={A1:'d',A2:'x',A3:'^',A4:'s',A5:'o'}
colors={A1:'k',A2:'r',A3:'g',A4:'b',A5:'c'}

iInfall=0

xmin=0.5
xmax=0.8
fMinInfall=1e-4
nMin=100.
alpha=0.96
beta=1.5
mustar=0.39

for H in [A1,A2,A3,A4,A5]:
  nMinInfall=fMinInfall*H.Mvir/H.mP
  f=(H.massTVV[:,iInfall]>nMinInfall)&(H.r/H.Rvir>xmin)&(H.r/H.Rvir<xmax)
  x=(1.*H.m/H.massTVV[:,iInfall])[f]
  print x[x>0].min()
  xbin=np.logspace(np.log10(x[x>0].min()),0, 1000)
  y=np.histogram(x,xbin)
  frac=y[0][::-1].cumsum()*1./len(x)
  frac=np.hstack([0, frac])[::-1]
  plt.plot(xbin, frac, '--', color=colors[H])
  plt.plot(xbin[xbin>nMin/nMinInfall], frac[xbin>nMin/nMinInfall], color=colors[H], label=H.name)
  if H is A1:
	par=st.norm.fit(np.log(x[x>0]))
	frac0=1.*np.sum(x>0)/len(x)
	#plt.plot(xbin, frac0*st.norm.sf(np.log(xbin),loc=par[0], scale=par[1]), '-', color=(.7,.7,.7), alpha=1, lw=5)
	plt.plot(xbin, frac0*st.norm.sf(np.log(xbin),loc=np.median(np.log(x[x>0])), scale=par[1]), '-', color=(.7,.7,.7), alpha=1, lw=5)
	print par,frac0
plt.semilogx()
plt.ylim([0,1])
#plt.xlim([4e-3,1])
plt.legend()
plt.xlabel(r'$m/m_0$')
plt.ylabel(r'Fraction$(>m/m_0)$')
#plt.savefig(outdir+'/StripCDF.eps')

##plot profiles
Rref=1
iInfall=0
fs=0.56
sigma=1.2

## mass fraction
nbin=30
xbin=np.logspace(-1.5, np.log10(2), nbin)
nMinInfall=1000.
H=A1
f=(H.massTVV[:,iInfall]>nMinInfall)#&(H.massTVV[:,iInfall]<nMinInfall*10)
x=H.r[f]/H.Rvir
y=(1.0*H.m/H.massTVV[:,iInfall]**1.)[f]
countAll,_=np.histogram(x, xbin)
plt.figure()
plt.plot(x,y, 'c.', alpha=0.1)
p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
plt.plot(xmid, p[1], 'r-', lw=2, label='Resolved')
plt.plot(xmid, p[[0,2]].T, 'r--',lw=2)
p,xmid=cross_section(x,y,xbin,100*(1-fs)+fs*np.array([(100-68.3)/2,50,(100+68.3)/2]))
plt.plot(xmid, p[1], 'k-', lw=2, label='Disruption-corrected')
plt.plot(xmid, p[[0,2]].T, 'k--', lw=2)
plt.yscale('log', nonposy='clip')
plt.xscale('log')
f0=(xmid<1.)&(xmid>0.01)&(p[1]>100./nMinInfall)

#suppression function
Hhost=A2
Hsub=A1
nMin=1000.
xbin=np.logspace(np.log10(0.5), np.log10(500*0.73), nbin)/Hhost.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=Hhost.get_host_density(xbin, Rref)
rSub,denSub,denSubRef,denSubErr=Hsub.get_sub_density((Hsub.m>nMin), xbin, Rref)
denRat=denSub/denHalo
f=(denRat>0)&(rSub<Rref)#&(rSub>0.1)
denRatErr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo
w=1/(denRatErr)[f]
pars=powerlaw_fit(rSub[f], denRat[f],w=w)
print pars
gamma=pars[0][1]
beta=gamma/alpha
x0=pars[0][0]
mustar=pars[0][2]

yscale=(p[1][f0]/(xmid[f0]/x0)**beta).mean()
yerr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)/alpha*denRat**(1/alpha)*yscale
plt.errorbar(rSub, denRat**(1/alpha)*yscale, yerr=yerr, fmt='bo',label=r'$(\rho_{\rm Sub}/\rho_{\rm Halo})^{1/\alpha}$')

ypred=(xmid)**beta*mustar*yscale
print 'R*:', x0*yscale**(-1./beta), 'beta:', beta, 'mu*:', mustar*yscale

#plt.errorbar(rSub, (denSub/denHalo)**(1/alpha)*yscale, yerr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo*yscale, fmt='bo',label=r'$\rho_{\rm Sub}/\rho_{\rm Halo}$') #label=r'$m/M_{200}>10^{%d}$'%(np.log10(Hsub.mP*mMin/Hsub.Mvir)))

plt.plot(xmid, ypred, 'g-', lw=2, label=r'$\beta=%.1f$'%(beta))
plt.plot(xmid, ypred*np.exp(sigma), 'g--', lw=2)
plt.plot(xmid, ypred/np.exp(sigma), 'g--', lw=2)

#tidal limits
Host=NFWHalo(183,15)
#Host=NFWHalo(83,9)
#Host=NFWHalo(100,7)
R=np.logspace(-2,0.3, 50)*Host.Rv
msat=1e-4*Host.M
h=NFWHalo(msat)
cref=h.C
k=1
plt.loglog(R/Host.Rv, Host.strip_func(h, R, k=k), 'm-', label='Tidal limit')
#plt.loglog(R/Host.Rv, Host.strip_func(NFWHalo(msat, cref*10**0.15), R, k=k), 'g--')
#plt.loglog(R/Host.Rv, Host.strip_func(NFWHalo(msat, cref/10**0.15), R, k=k), 'b--')

plt.legend(loc=4,fontsize=15)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$m/m_0$')
plt.xlim([0.04,3])
plt.ylim([1e-3,1])
plt.plot(plt.xlim(), 20./nMinInfall*np.array([1,1]), 'r:')
plt.plot(plt.xlim(), 100./nMinInfall*np.array([1,1]), 'r:')
plt.text(plt.xlim()[0]*1.1, 20./nMinInfall, '20', color='r')
plt.text(plt.xlim()[0]*1.1, 100./nMinInfall, '100', color='r')
#plt.savefig(outdir+'/StrippingResolved.'+H.name+'.nMin1e4.pdf')
#plt.savefig(outdir+'/StrippingResolved.'+H.name+'.pdf')

'''
## infall mass dependence
H=A1
hall=[]
plt.figure()
plt.errorbar(rSub, denSub/denHalo*yscale, yerr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo*yscale, fmt='bo',label=r'$\rho_{\rm Sub}/\rho_{\rm Halo}$') 
for nbin,c,lw,nMinInfall in [(40,'r',1,[1e3,1e4]), (30,'g',2,[1e4,1e5]), (20,'b',3,[1e5,1e7])]:
  xbin=np.logspace(-2, np.log10(2), nbin)
  f=(H.massTVV[:,iInfall]>nMinInfall[0])&(H.massTVV[:,iInfall]<nMinInfall[1])
  x=H.r[f]/H.Rvir
  y=(1.0*H.m/H.massTVV[:,iInfall]**1.)[f]
  countAll,_=np.histogram(x, xbin)
  #plt.plot(x,y, 'c.', alpha=0.1)
  #p,xmid=cross_section(x,y,xbin,44+0.56*np.array([(100-68.3)/2,50,(100+68.3)/2]))
  #plt.plot(xmid, p[1], 'k-', lw=lw, label='Disruption-corrected')
  #plt.plot(xmid, p[[0,2]].T, 'k--', lw=lw)
  p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
  htmp,=plt.plot(xmid, p[1], '-', color=c,lw=lw, label=r'$%s\sim %s$'%(fmtexp10(nMinInfall[0]*H.mP/H.Mvir*1.1),fmtexp10(nMinInfall[1]*H.mP/H.Mvir*1.1)))
  hall.append(htmp)
  #plt.plot(xmid, p[[0,2]].T, '--', color=c, lw=lw)
  plt.yscale('log', nonposy='clip')
  plt.xscale('log')
  plt.plot(plt.xlim(), 100./nMinInfall[0]*np.array([1,1]), ':', color=c, lw=lw)
  #plt.plot(plt.xlim(), 200./nMinInfall[0]*np.array([1,1]), ':')  

plt.legend(loc=2)#hall, (r'$10^{-6}\sim10^{-5}$',r'$10^{-5}\sim10^{-4}$',r'$10^{-4}\sim10^{-2}$'),
'''