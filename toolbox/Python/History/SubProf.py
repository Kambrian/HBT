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
nbin=30
Rconv=0.3
mMin=1000.
nMinInfall=0
iInfall=0

Hhost=A2
Hsub=A1

xbin=np.logspace(np.log10(0.5), np.log10(500*0.73), nbin)/Hhost.Rvir

plt.figure()
rHalo,denHalo,denHaloRef,denHaloErr=Hhost.get_host_density(xbin, Rref)
plt.plot(rHalo,denHalo, 'k', label=Hhost.name+'Halo')

rSub,denSub,denSubRef,denSubErr=Hsub.get_sub_density((Hsub.m>mMin)&(Hsub.massTVV[:,iInfall]>nMinInfall), xbin, Rref)
plt.errorbar(rSub,denSub, denSubErr, fmt='bo', label=Hsub.name+'Sub')
plt.yscale('log', nonposy='clip')
plt.xscale('log')

f=(denSub/denHalo>0)&(rSub<Rref)
w=1/(denSubErr/denSub)[f]
#w=None
polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],1, w=w, cov=True)
print polypar, np.sqrt(polyV.diagonal())
plt.plot(rHalo[denHalo>0], (denHalo*(rHalo/Rref)**polypar[0])[denHalo>0],'r--', label='Model')
#plt.plot(rHalo, denHalo*(rHalo/Rref)**1.2,'r--', label=r'$\beta=1.2$')

plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho/\rho(R_{200})$')
plt.legend()
#plt.savefig(outdir+'/'+Hsub.name+'subprof.eps')

plt.figure() 
plt.errorbar(rSub, denSub/denHalo, yerr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo, fmt='bo',label=Hsub.name+'sub') #label=r'$m/M_{200}>10^{%d}$'%(np.log10(Hsub.mP*mMin/Hsub.Mvir)))
plt.yscale('log', nonposy='clip')
plt.xscale('log')

plt.plot(rSub, np.exp(np.polyval(polypar, np.log(rSub))),'r--', label=r'$\gamma=%.2f\pm%.2f$'%(polypar[0],np.sqrt(polyV[0,0])))

REinasto=199*0.73/Hhost.Rvir
alpha=0.678
denEinasto=einasto(rSub/REinasto,alpha)/einasto(Rref/REinasto, alpha)
plt.plot(rSub[denHalo>0], (denEinasto/denHalo)[denHalo>0], 'g-',label='Einasto')
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
legend1=plt.legend(loc=2)
plt.ylim([1e-3,10])
#binlines=[]
#for i,fMin in enumerate([1e-6,1e-5,1e-4,1e-3]):
  #nMin=fMin*Hsub.Mvir/Hsub.mP
  #rSub,denSub,denSubRef,denSubErr=Hsub.get_sub_density((Hsub.m>nMin)&(Hsub.m<nMin*1e10), xbin, Rref)
  #f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
  #C=(denSub/denHalo)[f].mean()
  #f=denSub/denHalo>0
  #htmp,=plt.plot(rSub[f], (denSub/denHalo/C)[f], ':', linewidth=3, alpha=1, label=r'$10^{%d}<m/M_{200}<10^{%d}$'%(np.log10(fMin),np.log10(fMin)+1)) 
  #binlines.append(htmp)
#plt.legend(binlines, [a.get_label() for a in binlines], loc=4,fontsize=15)
#plt.gca().add_artist(legend1)

#plt.savefig(outdir+'/'+Hsub.name+'subprof_rat.eps')