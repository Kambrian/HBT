import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

A1=HaloData('A1')  
A2=HaloData('A2')
A3=HaloData('A3')
A4=HaloData('A4')
A5=HaloData('A5')
B4=HaloData('B4')

##plot profiles
Rref=1
iInfall=0
nbin=80
Rconv=0.3
mMin=1e-4*A2.Mvir

xbin=np.logspace(np.log10(Rconv), np.log10(500*0.73), nbin)/A2.Rvir
#xbin=np.linspace(Rconv, 500*0.73, nbin)/A2.Rvir
fmt={A1:'d',A2:'x',A3:'^',A4:'s',A5:'o'}
colors={A1:'k',A2:'r',A3:'g',A4:'b',A5:'c'}

## density profile
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, Rref)
for H in [A1,A2,A3,A4,A5]:
#for H in [A4]:
  rSub,denSub,denSubRef,denSubErr=H.get_sub_density(H.massTVV[:,iInfall]>mMin/H.mP, xbin, Rref)
  if H==A1:
	denSubRef0=denSubRef
	f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
	C=(denSub/denHalo)[f].mean()
	C0=C
  C=denSubRef0*C0/denSubRef #same normalization
  plt.plot(rSub[denSub>0],denSub[denSub>0]/C,'-',color=colors[H], alpha=0.5, label=H.name+r',$N_0>%s$'%fmtexp10(mMin/H.mP,'%.1f'))
  #plt.errorbar(rSub,denSub, denSubErr, fmt=fmt[H], color=colors[H], label=H.name)
#plt.yscale('log', nonposy='clip')
plt.yscale('log', nonposy='mask')
plt.xscale('log')
legend1=plt.legend(loc=1,fontsize=15)

halolines=[]
for H in [A2,A4]:
  rHalo,denHalo,denHaloRef,denHaloErr=H.get_host_density(xbin, Rref)
  htmp,=plt.plot(rHalo,denHalo, '--', linewidth=3, color=colors[H], alpha=1, label=H.name+'Halo')
  halolines.append(htmp)
  #plt.errorbar(rHalo,denHalo,denHaloErr, fmt='--', linewidth=3, color=colors[H], alpha=1, label=H.name+'Halo')
plt.legend(halolines, [a.get_label() for a in halolines], loc=3,fontsize=15)
plt.gca().add_artist(legend1)

##add B4
#mMinB=1e-4*B4.Mvir
#rHalo,denHalo,denHaloRef,denHaloErr=B4.get_host_density(xbin, Rref)
#rSub,denSub,denSubRef,denSubErr=B4.get_sub_density(B4.massTVV[:,iInfall]>mMinB/B4.mP, xbin, Rref)
#f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
#C=(denSub/denHalo)[f].mean()
#plt.plot(rSub[denSub>0],denSub[denSub>0]/C,'-',color='m', alpha=0.5, label=B4.name+r',$N_0>%s$'%fmtexp10(mMinB/B4.mP,'%.1f'))
#plt.plot(rHalo,denHalo, '--', linewidth=3, color='m', alpha=1, label=B4.name+'Halo')

plt.ylim([0.05,1e7])
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho/\rho(R_{200})$')
#plt.savefig(outdir+'/SubprofInfall.pdf')


##ratios
plt.figure()
H=A1
rSub,denSub,denSubRef,denSubErr=H.get_sub_density(H.massTVV[:,iInfall]>mMin/H.mP, xbin, Rref)
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, Rref)
f=(denSub/denHalo>0)&(rSub<1)&(rSub>1e-1)
polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],0, w=1/(denSubErr/denSub)[f], cov=True)
C=np.exp(polypar[0])
#print polypar, np.sqrt(polyV.diagonal())
plt.errorbar(rSub, denSub/denHalo/C, yerr=denSubErr/denHalo/C, fmt=fmt[H], color=colors[H], label=H.name+',$N_0>%s$'%fmtexp10(mMin/H.mP,'%.1f'))
#plt.plot(rSub, np.exp(np.polyval(polypar, np.log(rSub))),'r--', label=r'$\gamma=%.1f$'%polypar[0])  
plt.ylim([0.1,100])
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend()
plt.plot(plt.xlim(), [1,1], 'k--')  
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
#plt.savefig(outdir+'/SubprofInfall_rat.eps')

## mass bins
nbin=30
xbin=np.logspace(np.log10(Rconv), np.log10(500*0.73), nbin)/A2.Rvir
plt.figure()
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, Rref)
syms='^dso'
H=A1
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