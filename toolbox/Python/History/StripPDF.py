import numpy as np
import matplotlib.pyplot as plt
from matplotlib.finance import candlestick
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
import scipy.stats as st
from MbdIO import *
from myutils import skeleton,density_of_points,percentile_contour
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

A1=HaloData('A1')  
A2=HaloData('A2')
A3=HaloData('A3')
A4=HaloData('A4')
A5=HaloData('A5')
fmt={A1:'d',A2:'x',A3:'^',A4:'s',A5:'o'}
colors={A1:'k',A2:'r',A3:'g',A4:'b',A5:'c'}

iInfall=0

xmin=0.6
xmax=0.9
fMinInfall=1e-4
nMin=100.

for H in [A1,A2,A3,A4,A5]:
  nMinInfall=fMinInfall*H.Mvir/H.mP
  f=(H.massTVV[:,iInfall]>nMinInfall)&(H.r/H.Rvir>xmin)&(H.r/H.Rvir<xmax)
  x=(1.*H.m/H.massTVV[:,iInfall])[f]
  xbin=np.logspace(np.log10(x[x>0].min()),0, 1000)
  y=np.histogram(x,xbin)
  frac=y[0][::-1].cumsum()*1./len(x)
  frac=np.hstack([0, frac])[::-1]
  plt.plot(xbin, frac, '--', color=colors[H])
  plt.plot(xbin[xbin>nMin/nMinInfall], frac[xbin>nMin/nMinInfall], color=colors[H], label=H.name)
  if H is A1:
	par=st.norm.fit(np.log(x[x>0]))
	frac0=1.*np.sum(x>0)/len(x)
	plt.plot(xbin, frac0*st.norm.sf(np.log(xbin),loc=par[0], scale=par[1]), '-', color=(.7,.7,.7), alpha=1, lw=5)
plt.semilogx()
plt.ylim([0,1])
plt.legend()
plt.xlabel(r'$m/m_0$')
plt.ylabel(r'Fraction$(>m/m_0)$')
#plt.savefig(outdir+'/StripCDF.eps')

##plot profiles
Rref=1
iInfall=0

## mass fraction
nbin=30
xbin=np.logspace(-1, np.log10(2), nbin)
nMinInfall=3000
H=A1
f=(H.massTVV[:,iInfall]>nMinInfall)#&(H.massTVV[:,iInfall]<nMinInfall*10)
x=H.r[f]/H.Rvir
y=(1.0*H.m/H.massTVV[:,iInfall]**1.)[f]
countAll,_=np.histogram(x, xbin)
plt.figure()
plt.plot(x,y, 'c.', alpha=0.1)
p,xmid=cross_section(x,y,xbin,40+0.6*np.array([(100-68.3)/2,50,(100+68.3)/2]))
plt.plot(xmid, p[1], 'k-', lw=2, label='Disruption-corrected')
plt.plot(xmid, p[[0,2]].T, 'k--', lw=2)
p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
plt.plot(xmid, p[1], 'r-', lw=2, label='Resolved')
plt.plot(xmid, p[[0,2]].T, 'r--',lw=2)
#s=skeleton(x,y, xbin, alpha=0.9)
#plt.plot(s['x']['median'], s['y']['median'], 'k-')
#plt.plot(s['x']['median'], s['y']['CI'].T, 'r--')
#plt.plot(s['x']['median'], 1-1.*s['x']['hist']/countAll, 'b-')
plt.yscale('log', nonposy='clip')
plt.xscale('log')
sln=skeleton(np.log(x[y>0]),np.log(y[y>0]), np.log(xbin))
f=(~np.isnan(sln['x']['median']))&(sln['x']['median']<np.log(1.))&(sln['x']['median']>np.log(0.1))
p,Vp=np.polyfit(sln['x']['median'][f], sln['y']['median'][f], 1, w=1/(sln['y']['std']/np.sqrt(sln['x']['hist']))[f], cov=True)
ypred=np.exp(np.polyval(p, np.log(xbin)))
ysig=1.
plt.plot(xbin, ypred, 'g-', lw=2, label=r'$\gamma=%.1f\pm%.2f$'%(p[0], np.sqrt(Vp[0,0])))
plt.plot(xbin, ypred*np.exp(ysig), 'g--', lw=2)
plt.plot(xbin, ypred/np.exp(ysig), 'g--', lw=2)

Hhost=A2
Hsub=A1
nMin=1000.
xbin=np.logspace(np.log10(0.5), np.log10(500*0.73), nbin)/Hhost.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=Hhost.get_host_density(xbin, Rref)
rSub,denSub,denSubRef,denSubErr=Hsub.get_sub_density((Hsub.m>nMin), xbin, Rref)
f=(denSub/denHalo>0)&(rSub<Rref)
w=1/(denSubErr/denSub)[f]
polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],1, w=w, cov=True)
xref=np.log(0.8)
yscale=np.exp(np.polyval(p, xref)-np.polyval(polypar, xref))
plt.errorbar(rSub, denSub/denHalo*yscale, yerr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo*yscale, fmt='bo',label=r'$\rho_{\rm Sub}/\rho_{\rm Halo}$') #label=r'$m/M_{200}>10^{%d}$'%(np.log10(Hsub.mP*mMin/Hsub.Mvir)))
plt.legend(loc=4,fontsize=15)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$m/m_0$')
plt.xlim([0.05,3])
plt.ylim([1e-3,1])
plt.plot(plt.xlim(), 20./nMinInfall*np.array([1,1]), 'r:')
plt.plot(plt.xlim(), 200./nMinInfall*np.array([1,1]), 'r:')
plt.text(plt.xlim()[0]*1.1, 20./nMinInfall, '20', color='r')
plt.text(plt.xlim()[0]*1.1, 200./nMinInfall, '200', color='r')
plt.savefig(outdir+'/StrippingResolved.'+H.name+'.pdf')
'''
sln=skeleton(np.log10(x),np.log10(y), np.log10(xbin), alpha=0.68)
f=(~np.isnan(sln['x']['median']))&(sln['x']['median']>np.log10(0.1))&(sln['x']['median']<np.log10(1))
p,Vp=np.polyfit(sln['x']['median'][f], sln['y']['median'][f], 1, w=1/(sln['y']['std']/np.sqrt(sln['x']['hist']))[f], cov=True)
sln2=skeleton(np.log10(x),np.log10(y), np.log10(xbin), alpha=0.9)
quotes=np.hstack((sln['x']['median'][:,None], sln['y']['CI'].T, sln2['y']['CI'][::-1].T))
fig,ax=plt.subplots()
candlestick(ax, quotes, width=0.05, colorup='r')
plt.plot(sln['x']['median'], sln['y']['median'], 'k-')
plt.plot(np.log10(xbin), (np.polyval(p, np.log10(xbin))), 'g-', label=r'$\gamma=%.1f\pm%.2f$'%(p[0], np.sqrt(Vp[0,0])))
plt.plot(plt.xlim(), np.log10(nMin0/nMinInfall*np.array([1,1])), 'r:')
plt.legend(loc=4)
plt.xlabel(r'$\log(R/R_{\rm 200})$')
plt.ylabel(r'$\log(m/m_0)$')
#plt.savefig(outdir+'/StrippingResolved_candle.eps')

## infall mass bin
plt.figure()
for i,fMin in enumerate([1e-6,1e-5,1e-4,1e-3]):
  nMin=fMin*A1.Mvir/A1.mP
  f=(A1.massTVV[:,iInfall]>nMin)&(A1.massTVV[:,iInfall]<nMin*10)&(A1.m>0)
  print f.sum()
  sln=skeleton(np.log10(A1.r/A1.Mvir)[f],np.log10(1.0*A1.m/A1.massTVV[:,iInfall])[f], np.linspace(-1, 0.3, nbin), alpha=0.68)
  plt.plot(sln['x']['median'], sln['y']['mean'], '-o') #average log ratio.
  ff=(~np.isnan(sln['x']['median']))&(sln['x']['median']>-1)&(sln['x']['median']<0)
  p=np.polyfit(sln['x']['median'][ff], sln['y']['mean'][ff],1)
  print p
plt.xlabel(r'$\log(R/R_{\rm 200})$')
plt.ylabel(r'$\log(m/m_0)$')  

## resolution-free median
nbin=50
iInfall=0
xbin=np.logspace(-1, np.log10(3), nbin)
nMin0=100.
nMinInfall=2000
f=(A1.massTVV[:,iInfall]>nMinInfall)&(A1.massTVV[:,iInfall]<nMinInfall*1e10)#&(A1.flag>0)
x=A1.r[f]/A1.Rvir
y=(1.0*A1.m/A1.massTVV[:,iInfall])[f]
p,_=cross_section(x,y,xbin)
plt.figure()
plt.plot(x,y, 'c.', alpha=1)
s=skeleton(x,y, xbin, alpha=0.68)
plt.plot(s['x']['median'], s['y']['median'], 'k-')
plt.plot(s['x']['median'], p[:,[0,-1]], 'r--')
plt.plot(s['x']['median'], p[:,[1,-2]], 'b--')
plt.ylim([1e-3,10])
plt.yscale('log', nonposy='mask')
plt.xscale('log')
plt.plot(plt.xlim(), nMin0/nMinInfall*np.array([1,1]), 'r:')
#plt.plot(s['x']['median'], 10**-0.3*s['x']['median']**1.2, 'r-')
f=~np.isnan(s['x']['median'])&(s['x']['median']<3)&(s['y']['median']>1.*nMin0/nMinInfall)
if f.sum()>2:
  w=1/(np.log(s['y']['CI'][1]/s['y']['median'])/np.sqrt(s['x']['hist']))[f]
  #w=None
  par,Vpar=np.polyfit(np.log(s['x']['median'][f]), np.log(s['y']['median'][f]), 1, w=w, cov=True)
  plt.plot(xbin, np.exp(np.polyval(par, np.log(xbin))), 'g-', label=r'$\gamma=%.1f\pm %.2f$'%(par[0], np.sqrt(Vpar[0][0])))
plt.legend(loc=4)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$m/m_0$')
plt.axis([1e-1,3,1e-3,1])
#plt.savefig(outdir+'/StrippingPercentiles_FirstInfall.pdf')

## joint distribution at a given radius
#A1.m[A1.m==0]=1
for i,rmin in enumerate(np.logspace(-2, 0, 6)):
  f=(A1.r/A1.Rvir>rmin)&(A1.r/A1.Rvir<rmin*10**(2/5.))
  x=np.log10(A1.massTVV[f,iInfall])
  y=np.log10(A1.m[f]*1./A1.massTVV[f, iInfall])
  ff=A1.m[f]>0
  if ff.sum()>2:
	plt.figure()
	plt.plot(x[ff], y[ff], '.', color='b', alpha=0.1)
	X,Y,Z=density_of_points(np.array([x[ff],y[ff]]), bins=30, method='kde')
	percentile_contour(X,Y,Z, [0.1,0.3,0.5,0.7,0.9], colors='k')
	p,xmid=cross_section(x,y,20)
	plt.plot(xmid, p[:,2], 'r-')
	plt.plot(xmid, p[:,[0,-1]], 'g-')
	plt.plot(xmid, p[:,[1,-2]], 'b-')
	try:
	  xmin=xmid[p[:,2]+xmid>np.log10(20)].min()
	  ymed=np.median(y[x>xmin])
	  print np.log10(rmin)+0.2, ymed
	  plt.plot(plt.xlim(), ymed*np.array([1,1.]), 'c--')
	except:
	  pass
  plt.ylabel(r'$\log(m/m_0)$')
  plt.xlabel(r'$\log(m_0)$') 
 ''' 
