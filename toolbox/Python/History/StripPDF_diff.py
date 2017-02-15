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

def strip_pdf(lnx,mu,sig,mumax,f0):
  lnmu=np.log(mu)
  lnmumax=np.log(mumax)
  df=st.norm(loc=lnmu, scale=sig)
  lnx=np.asarray(lnx)
  y=np.zeros_like(lnx)
  y[lnx<lnmumax]=f0*df.pdf(lnx[lnx<lnmumax])/df.cdf(lnmumax)
  return y

def strip_cdf(lnx,mu,sig,mumax,f0):
  lnmu=np.log(mu)
  lnmumax=np.log(mumax)
  df=st.norm(loc=lnmu, scale=sig)
  lnx=np.asarray(lnx)
  y=np.zeros_like(lnx)
  y[lnx<lnmumax]=f0*(df.sf(lnx[lnx<lnmumax])-df.sf(lnmumax))/df.cdf(lnmumax)
  return y

fig,ax = plt.subplots(2, sharex=True, figsize=(8,12))  
for H in [A1,A2,A3,A4,A5]:
  nMinInfall=fMinInfall*H.Mvir/H.mP
  f=(H.massTVV[:,iInfall]>nMinInfall)&(H.r/H.Rvir>xmin)&(H.r/H.Rvir<xmax)
  x=(1.*H.m/H.massTVV[:,iInfall])[f]
  print x[x>0].min()
  xbin=np.logspace(np.log10(x[x>0].min()),0, 1000)
  y=np.histogram(x,xbin)
  frac=y[0][::-1].cumsum()*1./len(x)
  frac=np.hstack([0, frac])[::-1]
  ax[0].plot(xbin, frac, '--', color=colors[H])
  ax[0].plot(xbin[xbin>nMin/nMinInfall], frac[xbin>nMin/nMinInfall], color=colors[H], label=H.name)
  if H is A1:
	par=list(st.norm.fit(np.log(x[x>0])))
	frac0=1.*np.sum(x>0)/len(x)
	ax[0].plot(xbin, frac0*st.norm.sf(np.log(xbin),loc=par[0], scale=par[1]), '-', color=(.7,.7,.7), alpha=1, lw=5, label='fit')
	loc=np.median(np.log(x[x>0]))
	mumax=min(np.log(4.2)+loc, 0)
	df=st.norm(loc=loc, scale=par[1])
	df=st.norm(loc=loc+0., scale=1.25)
	strip_cdffit=lambda lnx,mu,sig,mumax: strip_cdf(lnx,mu,sig,mumax,frac0)
	par2=curve_fit(strip_cdffit, np.log(xbin), frac, p0=[np.exp(loc), 1.2, np.exp(mumax)])#, frac0])
hall=[];
for H in [A1,A2,A3,A4,A5]:
  nMinInfall=fMinInfall*H.Mvir/H.mP
  f=(H.massTVV[:,iInfall]>nMinInfall)&(H.r/H.Rvir>xmin)&(H.r/H.Rvir<xmax)
  x=(1.*H.m/H.massTVV[:,iInfall])[f]
  xbin=np.logspace(np.log10(x[x>0].min()),0, 10)
  print H.m[H.m>0].min(), x[x>0].min()
  y=np.histogram(x,xbin)
  frac=y[0]*1./len(x)/np.log(xbin[1]/xbin[0])
  xcen=xbin[:-1]*np.sqrt(xbin[1]/xbin[0])
  ax[1].plot(xcen, frac, '--', color=colors[H])
  ax[1].plot(xcen[xbin[:-1]>nMin/nMinInfall], frac[xbin[:-1]>nMin/nMinInfall], color=colors[H], label=H.name)
  if H is A1:
	par=list(st.norm.fit(np.log(x[x>0])))
	frac0=1.*np.sum(x>0)/len(x)
	xfit=np.logspace(np.log10(xbin[0]), 0,50)
	ax[1].plot(xfit, frac0*st.norm.pdf(np.log(xfit),loc=par[0], scale=par[1]), '-', color=(.7,.7,.7), alpha=1, lw=5, label='fit')
ax[0].semilogx()
ax[1].semilogx()
ax[0].set_ylabel(r'Fraction$(>\mu)$')	
ax[1].set_xlabel(r'$\mu=m/m_{\rm acc}$')
ax[1].set_ylabel(r'$\rm{d}P/\rm{d}\ln\mu$')
ax[1].legend(loc=2)
fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
ax[1].set_ylim([0,0.28])
ax[0].set_ylim([0,0.8])
#nbins = 7 #len(ax[0].get_yticklabels())
#plt.setp([a.yaxis for a in ax], major_locator=MaxNLocator(nbins=nbins, prune='lower',symmetric=True))
#plt.savefig(outdir+'/StripPDF_combined.eps')
