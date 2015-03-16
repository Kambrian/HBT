import numpy as np
import matplotlib.pyplot as plt
from matplotlib.finance import candlestick
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
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

##plot profiles
Rref=1
iInfall=0

## mass fraction
nbin=20
xbin=np.logspace(-1, np.log10(2), nbin)
nMin0=20.
nMinInfall=2000
countAll,_=np.histogram(A1.r[A1.massTVV[:,iInfall]>nMinInfall]/A1.Rvir, xbin)
f=(A1.m>nMin0)&(A1.massTVV[:,iInfall]>nMinInfall)
x=A1.r[f]/A1.Rvir
y=(1.0*A1.m/A1.massTVV[:,iInfall])[f]
plt.figure()
plt.plot(x,y, 'c.', alpha=0.2)
s=skeleton(x,y, xbin, alpha=0.9)
plt.plot(s['x']['median'], s['y']['median'], 'k-')
plt.plot(s['x']['median'], s['y']['CI'].T, 'r--')
plt.plot(s['x']['median'], 1-1.*s['x']['hist']/countAll, 'b-')
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.plot(plt.xlim(), nMin0/nMinInfall*np.array([1,1]), 'r:')
#plt.plot(s['x']['median'], 10**-0.3*s['x']['median']**1.2, 'r-')
sln=skeleton(np.log(x),np.log(y), np.log(xbin[(xbin<1.)&(xbin>0.1)]))
f=~np.isnan(sln['x']['median'])
p,Vp=np.polyfit(sln['x']['median'][f], sln['y']['median'][f], 1, w=1/(sln['y']['std']/np.sqrt(sln['x']['hist']))[f], cov=True)
plt.plot(xbin, np.exp(np.polyval(p, np.log(xbin))), 'g-', label=r'$\gamma=%.1f\pm%.2f$'%(p[0], np.sqrt(Vp[0,0])))
plt.legend(loc=4)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$m/m_0$')
plt.xlim([0.07,3])
#plt.savefig(outdir+'/StrippingResolved.pdf')

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
  
