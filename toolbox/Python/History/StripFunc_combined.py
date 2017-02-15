import numpy as np
import matplotlib.pyplot as plt
from matplotlib.finance import candlestick
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
import scipy.stats as st
from scipy.optimize import curve_fit
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

##plot profiles
Rref=1
iInfall=0
fs=0.55
sigma=1.1

## mass fraction
def plot_strip(plot_model=True):
  nbin=30
  xbin=np.logspace(-1.5, np.log10(2), nbin)
  nMinInfall=1000.
  H=A1
  f=(H.massTVV[:,iInfall]>nMinInfall)
  x=H.r[f]/H.Rvir
  y=(1.0*H.m/H.massTVV[:,iInfall]**1.)[f]
  countAll,_=np.histogram(x, xbin)
  if not plot_model:
	plt.plot(x,y, 'c.', alpha=0.1)
  p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
  h,=plt.plot(xmid, p[1], 'r-', lw=2)#, label='Resolved')
  if not plot_model:
	h.set_label('Resolved')
  plt.plot(xmid, p[[0,2]].T, 'r--',lw=2)
  p,xmid=cross_section(x,y,xbin,100*(1-fs)+fs*np.array([(100-68.3)/2,50,(100+68.3)/2]))
  h,=plt.plot(xmid, p[1], 'k-', lw=2)
  if not plot_model:
	h.set_label('Disruption-corrected')
  plt.plot(xmid, p[[0,2]].T, 'k--', lw=2)
  plt.yscale('log', nonposy='clip')
  plt.xscale('log')

  if plot_model:
	f0=(xmid<1.)&(xmid>0.01)&(p[1]>100./nMinInfall)
	#suppression function
	nMin=1000.
	xbin=np.logspace(np.log10(0.5), np.log10(500*0.73), nbin)/H.Rvir
	rHalo,denHalo,denHaloRef,denHaloErr=H.get_host_density(xbin, Rref)
	rSub,denSub,denSubRef,denSubErr=H.get_sub_density((H.m>nMin), xbin, Rref)
	denRat=denSub/denHalo
	f=(denRat>0)&(rSub<Rref)&(rSub>0.1)
	denRatErr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo
	w=None
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

	plt.plot(xmid, ypred, 'g-', lw=2, label=r'$\beta=%.1f$'%(beta))
	plt.plot(xmid, ypred*np.exp(sigma), 'g--', lw=2)
	plt.plot(xmid, ypred/np.exp(sigma), 'g--', lw=2)
  else:
	Host=NFWHalo(183,15)
	R=np.logspace(-2,0.3, 50)*Host.Rv
	msat=1e-4*Host.M
	h=NFWHalo(msat)
	cref=h.C
	k=1
	plt.loglog(R/Host.Rv, Host.strip_func(h, R, k=k), 'm-', label='Tidal limit')

  plt.legend(loc=4,fontsize=15)
  plt.xlabel(r'$R/R_{\rm 200}$')
  plt.xlim([0.04,3])
  plt.ylim([1e-3,1])
  plt.plot(plt.xlim(), 20./nMinInfall*np.array([1,1]), 'r:')
  plt.plot(plt.xlim(), 100./nMinInfall*np.array([1,1]), 'r:')
  if not plot_model:
	plt.text(plt.xlim()[0]*1.1, 20./nMinInfall, '20', color='r')
	plt.text(plt.xlim()[0]*1.1, 100./nMinInfall, '100', color='r')
fig=plt.figure(figsize=(16,8))
ax1=plt.subplot(121)
plot_strip(False)
plt.ylabel(r'$m/m_{\rm acc}$')
ax2=plt.subplot(122)
plot_strip(True)
fig.subplots_adjust(wspace=0)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.savefig(outdir+'/StrippingResolved.combined.pdf')
