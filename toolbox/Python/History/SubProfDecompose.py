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

##plot profiles
Rref=1
iInfall=0
nbin=80
Rconv=0.3
mMin=1e-4*A2.Mvir

xbin=np.logspace(np.log10(Rconv), np.log10(500*0.73), nbin)/A2.Rvir
#xbin=np.linspace(Rconv, 500*0.73, nbin)/A2.Rvir
fmt={A1:'d',A2:'x',A3:'^',A4:'s',A5:'o'}
colors={A1:'m',A2:'r',A3:'g',A4:'b',A5:'c'}

## resolved
nbin=20
xbin=np.logspace(np.log10(Rconv), np.log10(500*0.73), nbin)/A2.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, Rref)
rSub,denSub,denSubRef0,denSubErr=A1.get_sub_density((A1.massTVV[:,iInfall]>mMin/A1.mP), xbin, Rref)
f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
C=(denSub/denHalo)[f].mean()
denSubRef0*=C
plt.figure()
hall=[]
f=denSub/denHalo>0
htmp,=plt.plot(rSub[f], (denSub/denHalo/C)[f], 'k-', lw=5, alpha=0.4, label='A1 all')
y0=denSub/denHalo/C
hall.append(htmp)
for H in [A1,A2,A3,A4,A5]:
  rSub,denSub,denSubRef,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP/1)&(H.m>0), xbin, Rref)
  #f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
  #C=(denSub/denHalo)[f].mean() #normalize individually.
  C=denSubRef0/denSubRef #normalize by A1all
  f=denSub/denHalo>0
  y=(denSub/denHalo/C/y0)[(rSub>0.6)&(rSub<1.)].mean()
  print y
  if H!=A1:
	htmp,=plt.plot(rSub[f], (denSub/denHalo/C)[f], fmt[H]+'-', color=colors[H], label=H.name)
  else:
	htmp=plt.errorbar(rSub[f], (denSub/denHalo/C)[f], (denSubErr/denHalo/C)[f], fmt=fmt[H]+'-', color=colors[H], alpha=0.7, label=H.name)
	#plt.xscale('log')
	#y=0.6
	#plt.plot(plt.xlim(), [y,y], 'k--')
  hall.append(htmp)
plt.ylim([5e-2,3])
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend(hall, [h.get_label() for h in hall], loc=3, fontsize=15)
plt.plot(plt.xlim(), [1,1], 'k--')  
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
#plt.savefig(outdir+'/SubprofInfall_resolved.rat.samenorm.pdf')

plt.figure()
nbin=80
xbin=np.logspace(np.log10(Rconv), np.log10(500*0.73), nbin)/A2.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, Rref)
plt.plot(rHalo, denHalo, 'k-')
nbin=50
xbin=np.logspace(np.log10(Rconv), np.log10(500*0.73), nbin)/A2.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, Rref)
rSub,denSub,denSubRef0,denSubErr=A1.get_sub_density((A1.massTVV[:,iInfall]>mMin/A1.mP), xbin, Rref)
f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
C=(denSub/denHalo)[f].mean()
denSubRef0*=C
for H in [A1,A2,A3,A4,A5]:
  rSub,denSub,denSubRef,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP/1)&(H.m>0), xbin, Rref)
  #f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
  #C=(denSub/denHalo)[f].mean() #normalize individually
  C=denSubRef0/denSubRef #normalize by A1 all
  plt.plot(rSub[denSub>0],denSub[denSub>0]/C,'-',color=colors[H], alpha=1, label=H.name+r',$N_0>%s$'%fmtexp10(mMin/H.mP,'%.1f'))
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend(loc=3,fontsize=15)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho(R)/\rho(R_{200})$')
#plt.savefig(outdir+'/SubprofInfall_resolved.samenorm.eps')

#decompose A2

plt.figure()
H=A1
nbin=80
xbin=np.logspace(np.log10(Rconv), np.log10(500*0.73), nbin)/A2.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, Rref)
plt.plot(rHalo, denHalo, 'k-')
rSub,denSub,denSubRef0,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP), xbin, Rref)
f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
C=(denSub/denHalo)[f].mean()
plt.plot(rSub[denSub>0],denSub[denSub>0]/C,'b-',alpha=0.3, linewidth=4, label=r'$m_0/M_{200}>%s$'%fmtexp10(mMin/A2.Mvir))
rSub,denSub,denSubRef1,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP)&(H.m>0), xbin, Rref)
C1=C*denSubRef0/denSubRef1
plt.plot(rSub[denSub>0],denSub[denSub>0]/C1,'r--', alpha=1, label='Resolved')
rSub,denSub,denSubRef2,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP)&(H.m==0), xbin, Rref)
C2=C*denSubRef0/denSubRef2
plt.plot(rSub[denSub>0],denSub[denSub>0]/C2,'g:', alpha=1, label='Orphan')

#plt.errorbar(rSub,denSub, denSubErr, fmt=fmt[H], color=colors[H], label=H.name)
#plt.yscale('log', nonposy='clip')
#plt.ylim([5e-2,3])
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend(loc=1)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho(R)/\rho(R_{200})$')
#plt.savefig(outdir+'/SubprofInfall_Decompose.pdf')

## mass dependence (only depend on stripping fraction?) todo.............
nbin=20
xbin=np.logspace(np.log10(Rconv), np.log10(500*0.73), nbin)/A2.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, Rref)
hall=[]
plt.figure()
for H in [A1,A2,A3,A4,A5]:
  rSub,denSub,denSubRef,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP/1)&(H.m>0), xbin, Rref)
  #plt.plot(rSub[denSub>0],denSub[denSub>0],'-',color=colors[H], alpha=0.5, label=H.name+r',$N_0>%s$'%fmtexp10(mMin/H.mP,'%.1f'))
  f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
  C=(denSub/denHalo)[f].mean()
  polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],0, w=1/(denSubErr/denSub)[f], cov=True)
  #C=np.exp(polypar[0])
  f=denSub/denHalo>0
  if H!=A1:
	htmp,=plt.plot(rSub[f], (denSub/denHalo/C)[f], fmt[H]+'-', label=H.name)
  else:
	htmp=plt.errorbar(rSub[f], (denSub/denHalo/C)[f], (denSubErr/denHalo/C)[f], fmt=fmt[H]+'-', alpha=0.7, label=H.name)
  hall.append(htmp)
  #plt.errorbar(rSub,denSub, denSubErr, fmt=fmt[H], color=colors[H], label=H.name)
#plt.yscale('log', nonposy='clip')
plt.ylim([5e-2,3])
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend(hall, [h.get_label() for h in hall], loc=3, fontsize=15)
plt.plot(plt.xlim(), [1,1], 'k--')  
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
#plt.savefig(outdir+'/SubprofInfall_resolved.rat.pdf')
