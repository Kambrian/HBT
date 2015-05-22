import sys,os
sys.path.append(os.path.abspath('..'))
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

AqHalo=[HaloData('Aq'+h+'2') for h in 'ABCDEF']
PhHalo=[HaloData('Ph'+h+'2') for h in 'ABCDEF']

iInfall=0

rmin=0.5
rmax=0.8
fMinInfall=1e-4
nMin=100.

xbin=np.logspace(-3,0, 1000)
def get_cum_distr(H, fMinInfall, xbin, rmin, rmax):
  nMinInfall=fMinInfall*H.Mvir/H.mP
  f=(H.massTVV[:,iInfall]>nMinInfall)&(H.r/H.Rvir>rmin)&(H.r/H.Rvir<rmax)
  x=(1.*H.m/H.massTVV[:,iInfall])[f]
  #print x[x>0].min()
  y=np.histogram(x,xbin)
  frac=y[0][::-1].cumsum()*1./len(x)
  frac=np.hstack([0, frac])[::-1]
  return frac,x

def batchCDF(HList, color='r',label=None):
  fracAll=[]
  xmin=[]
  xall=[]
  for H in HList:
	nMinInfall=fMinInfall*H.Mvir/H.mP
	frac,x=get_cum_distr(H, fMinInfall, xbin, rmin, rmax)
	fracAll.append(frac)
	xmin.append(nMin/nMinInfall)
	xall.extend(x)
	frac0=1.*np.sum(x>0)/len(x)
	print frac0, len(x)
  x=np.array(xall)
  xmin=np.max(xmin)
  fracAll=np.array(fracAll)
  frac=fracAll.mean(0)
  fraclim=[fracAll.min(0),fracAll.max(0)]
  fracstd=fracAll.std(0)
  
  plt.plot(xbin,frac, '-', color=color)
  plt.plot(xbin[xbin>xmin], frac[xbin>xmin], '-', color=color, label=label)
  plt.fill_between(xbin,frac+fracstd,frac-fracstd,color=color,alpha=0.3)
  #plt.fill_between(xbin,fraclim[0],fraclim[1],color='g',alpha=0.3)
  par=st.norm.fit(np.log(x[x>0]))
  frac0=1.*np.sum(x>0)/len(x)
  #plt.plot(xbin, frac0*st.norm.sf(np.log(xbin),loc=par[0], scale=par[1]), '-', color='k')
  print par,frac0
  print 'fs=',frac[0],'+-',fracstd[0]
  return x

Aq=batchCDF(AqHalo, 'r', 'Aquarius')
Ph=batchCDF(PhHalo, 'g', 'Phoenix')

print 'All:'
x=np.hstack([Aq,Ph])
par=st.norm.fit(np.log(x[x>0]))
frac0=1.*np.sum(x>0)/len(x)
print par, frac0
#plt.plot(xbin, frac0*st.norm.sf(np.log(xbin),loc=par[0], scale=par[1]), '-', color='k')
plt.plot(xbin, frac0*st.norm.sf(np.log(xbin),loc=np.median(np.log(x[x>0])), scale=1.1), '-', color='k')

#if H is PhA2:
  #par=st.norm.fit(np.log(x[x>0]))
  #frac0=1.*np.sum(x>0)/len(x)
  #plt.plot(xbin, frac0*st.norm.sf(np.log(xbin),loc=par[0], scale=par[1]), '-', color=l.get_color(), alpha=0.2, lw=5)
  #print par,frac0
plt.semilogx()
plt.ylim([0,1])
plt.legend()
plt.xlabel(r'$m/m_0$')
plt.ylabel(r'Fraction$(>m/m_0)$')
#plt.savefig(outdir+'/StripCDF.All.pdf')
'''
##plot profiles
Rref=1
iInfall=0

## mass fraction
nbin=30
xbin=np.logspace(-1.5, np.log10(2), nbin)
nMinInfall=1000.
H=PhB2
f=(H.massTVV[:,iInfall]>nMinInfall)#&(H.massTVV[:,iInfall]<nMinInfall*10)
x=H.r[f]/H.Rvir
y=(1.0*H.m/H.massTVV[:,iInfall]**1.)[f]
countAll,_=np.histogram(x, xbin)
plt.figure()
plt.plot(x,y, 'c.', alpha=0.1)
p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
plt.plot(xmid, p[1], 'r-', lw=2, label='Resolved')
plt.plot(xmid, p[[0,2]].T, 'r--',lw=2)
p,xmid=cross_section(x,y,xbin,40+0.6*np.array([(100-68.3)/2,50,(100+68.3)/2]))
plt.plot(xmid, p[1], 'k-', lw=2, label='Disruption-corrected')
plt.plot(xmid, p[[0,2]].T, 'k--', lw=2)
#s=skeleton(x,y, xbin, alpha=0.9)
#plt.plot(s['x']['median'], s['y']['median'], 'k-')
#plt.plot(s['x']['median'], s['y']['CI'].T, 'r--')
#plt.plot(s['x']['median'], 1-1.*s['x']['hist']/countAll, 'b-')
plt.yscale('log', nonposy='clip')
plt.xscale('log')
sln=skeleton(np.log(x[y>0]),np.log(y[y>0]), np.log(xbin))
f0=(~np.isnan(sln['x']['median']))&(sln['x']['median']<np.log(1.))&(sln['x']['median']>np.log(0.01))&(sln['x']['hist']>1)&(sln['y']['median']>np.log(100./nMinInfall))
#w=1/(sln['y']['std']/np.sqrt(sln['x']['hist']))[f]
#p,Vp=np.polyfit(np.log(xmid[f]), np.log(p[1][f]), 1, w=w, cov=True)
#print 'R*:', np.exp(-p[1]/p[0]), 'beta:', p[0]

#suppression function
Hhost=PhB2
Hsub=PhB2
nMin=1000.
xbin=np.logspace(np.log10(0.5), np.log10(5000), nbin)/Hhost.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=Hhost.get_host_density(xbin, Rref)
rSub,denSub,denSubRef,denSubErr=Hsub.get_sub_density((Hsub.m>nMin), xbin, Rref)
f=(denSub/denHalo>0)&(rSub<Rref)
w=None #1/(denSubErr/denSub)[f]
polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],1, w=w, cov=True)
print 'gamma:', polypar[0]
yscale=(p[1][f0]/np.exp(np.polyval(polypar,np.log(xmid[f0])))).mean()
#xref=np.log(0.8)
#yscale=np.exp(np.polyval(p, xref)-np.polyval(polypar, xref))
plt.errorbar(rSub, denSub/denHalo*yscale, yerr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo*yscale, fmt='bo',label=r'$\rho_{\rm Sub}/\rho_{\rm Halo}$') #label=r'$m/M_{200}>10^{%d}$'%(np.log10(Hsub.mP*mMin/Hsub.Mvir)))

polypar[1]+=np.log(yscale)
ypred=np.exp(np.polyval(polypar, np.log(xbin)))
print 'R*:', np.exp(-polypar[1]/polypar[0]), 'beta:', polypar[0]
ysig=1.
plt.plot(xbin, ypred, 'g-', lw=2, label=r'$\beta=%.1f$'%(polypar[0]))
plt.plot(xbin, ypred*np.exp(ysig), 'g--', lw=2)
plt.plot(xbin, ypred/np.exp(ysig), 'g--', lw=2)

#tidal limits
Host=NFWHalo(6.57e4,5.96)
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
#plt.savefig(outdir+'/StrippingResolved.'+H.name+'.pdf')

## infall mass dependence
H=PhA2
plt.figure()
for nbin,c,lw,nMinInfall in [(50,'r',1,[1e3,1e4]), (40,'g',2,[1e4,1e5]), (30,'b',3,[1e5,1e6]), (15,'c',4,[1e6,1e100])]:
  xbin=np.logspace(-2, np.log10(2), nbin)
  f=(H.massTVV[:,iInfall]>nMinInfall[0])&(H.massTVV[:,iInfall]<nMinInfall[1])
  x=H.r[f]/H.Rvir
  y=(1.0*H.m/H.massTVV[:,iInfall]**1.)[f]
  countAll,_=np.histogram(x, xbin)
  #plt.plot(x,y, 'c.', alpha=0.1)
  #p,xmid=cross_section(x,y,xbin,40+0.6*np.array([(100-68.3)/2,50,(100+68.3)/2]))
  #plt.plot(xmid, p[1], 'k-', lw=lw, label='Disruption-corrected')
  #plt.plot(xmid, p[[0,2]].T, 'k--', lw=lw)
  p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
  plt.plot(xmid, p[1], '-', color=c,lw=lw, label='Resolved')
  plt.plot(xmid, p[[0,2]].T, '--', color=c, lw=lw)
  plt.yscale('log', nonposy='clip')
  plt.xscale('log')
  plt.plot(plt.xlim(), 100./nMinInfall[0]*np.array([1,1]), ':', color=c)
  #plt.plot(plt.xlim(), 200./nMinInfall[0]*np.array([1,1]), ':')
  
plt.errorbar(rSub, denSub/denHalo*yscale, yerr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo*yscale, fmt='bo',label=r'$\rho_{\rm Sub}/\rho_{\rm Halo}$') #label=r'$m/M_{200}>10^{%d}$'%(np.log10(Hsub.mP*mMin/Hsub.Mvir)))  
'''