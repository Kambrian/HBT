import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

H=HaloData(sys.argv[1])

##plot profiles
Rref=1
iInfall=0
nbin=80
Rconv=0.3
fMin=1e-6

xbin=np.logspace(-2,0.4, nbin)

## density profile
rHalo,denHalo,denHaloRef,denHaloErr=H.get_host_density(xbin, Rref)
rSub,denSub,denSubRef,denSubErr=H.get_sub_density(H.massTVV[:,iInfall]>fMin*H.Mvir/H.mP, xbin, Rref)
plt.errorbar(rSub[denSub>0],denSub[denSub>0],denSubErr[denSub>0], fmt='o', label=H.name+r',$N_0>%s$'%fmtexp10(fMin*H.Mvir/H.mP,'%.1f'))
plt.yscale('log', nonposy='clip')
#plt.yscale('log', nonposy='mask')
plt.xscale('log')
htmp,=plt.plot(rHalo,denHalo, '--', linewidth=3, alpha=0.5, label=H.name+'Halo')
plt.legend(loc=1,fontsize=15)

plt.plot([0.6/H.Rvir,0.6/H.Rvir], [1,2])

plt.ylim([0.05,1e7])
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho(R)/\rho(R_{200})$')
#plt.savefig(outdir+'/SubprofInfall.pdf')


##ratios
plt.figure()
f=(denSub/denHalo>0)&(rSub<1)&(rSub>1e-1)
polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],0, w=1/(denSubErr/denSub)[f], cov=True)
C=np.exp(polypar[0])
plt.errorbar(rSub, denSub/denHalo/C, yerr=denSubErr/denHalo/C, label=H.name+',$N_0>%s$'%fmtexp10(fMin*H.Mvir/H.mP,'%.1f'))
plt.ylim([0.1,100])
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend()
plt.plot(plt.xlim(), [1,1], 'k--')  
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
#plt.savefig(outdir+'/SubprofInfall_rat.eps')

## mass bins
plt.figure()
syms='^dso'
for i,fMin in enumerate([1e-6,1e-5,1e-4,1e-3]):
  nMin=fMin*H.Mvir/H.mP
  rSub,denSub,denSubRef,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>nMin)&(H.massTVV[:,iInfall]<nMin*10), xbin, Rref)
  f=(denSub/denHalo>0)&(rSub<1)&(rSub>1e-1)
  polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],0, w=1/(denSubErr/denSub)[f], cov=True)
  C=np.exp(polypar[0])
  f=denSub/denHalo>0
  if fMin<1e-3:
	plt.plot(rSub[f], (denSub/denHalo/C)[f], '-'+syms[i], label=r'$10^{%d}<m_0/M_{200}<10^{%d}$'%(np.log10(fMin),np.log10(fMin)+1))
  else:
	plt.errorbar(rSub[f], (denSub/denHalo/C)[f], yerr=(denSubErr/denHalo/C)[f], fmt='-'+syms[i], alpha=0.4, label=r'$10^{%d}<m_0/M_{200}<10^{%d}$'%(np.log10(fMin),np.log10(fMin)+1))
#plt.plot(rSub, np.exp(np.polyval(polypar, np.log(rSub))),'r--', label=r'$\gamma=%.1f$'%polypar[0])  
plt.ylim([0.1,100])
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend()
plt.plot(plt.xlim(), [1,1], 'k--')  
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
#plt.savefig(outdir+'/SubprofInfall_bin.pdf')