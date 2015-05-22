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

A1=HaloData('A1')
PhA2=HaloData('PhA2')
PhB2=HaloData('PhB2')
PhC2=HaloData('PhC2')
PhD2=HaloData('PhD2')
PhE2=HaloData('PhE2')
PhF2=HaloData('PhF2')
PhG2=HaloData('PhG2')
#fmt={A1:'d',A2:'x',A3:'^',A4:'s',A5:'o'}
#colors={A1:'k',A2:'r',A3:'g',A4:'b',A5:'c'}

iInfall=0

xmin=0.5
xmax=0.8
fMinInfall=1e-4
nMin=100.

for H in [A1,PhA2]:#,PhB2,PhC2,PhD2,PhE2,PhF2,PhG2]:
  nMinInfall=fMinInfall*H.Mvir/H.mP
  f=(H.massTVV[:,iInfall]>nMinInfall)&(H.r/H.Rvir>xmin)&(H.r/H.Rvir<xmax)
  x=(1.*H.m/H.massTVV[:,iInfall])[f]
  print x[x>0].min()
  xbin=np.logspace(np.log10(x[x>0].min()),0, 30)
  y=np.histogram(np.log(x),np.log(xbin),density=True)
  y=np.hstack([y[0],0])
  l,=plt.step(xbin, y, '--')
  plt.step(xbin[xbin>nMin/nMinInfall], y[xbin>nMin/nMinInfall], color=l.get_color(), label=H.name)
  if H is PhA2:
	par=st.norm.fit(np.log(x[x>0]))
	frac0=1.*np.sum(x>0)/len(x)
	plt.plot(xbin, st.norm.pdf(np.log(xbin),loc=par[0], scale=par[1]), '-', color=l.get_color(), alpha=0.2, lw=5)
	print par,frac0
plt.semilogx()
plt.ylim([0,1])
plt.legend(loc=2)
plt.xlabel(r'$m/m_0$')
plt.ylabel(r'Fraction$(>m/m_0)$')
#plt.savefig(outdir+'/StripCDF.eps')
