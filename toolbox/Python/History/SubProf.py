import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'
  
A1=HaloData('AqA1')
A2=HaloData('AqA2')
A3=HaloData('AqA3')
A4=HaloData('AqA4')
A5=HaloData('AqA5')
A2H=HaloData('AqA2', None)
A4H=HaloData('AqA4', None)

##plot profiles
Rref=1
nbin=30
Rconv=0.3
fMin=1e-6
nMinInfall=0
iInfall=0

Hsub=HaloData('AqA1')
Hhost=Hsub

xbin=np.logspace(np.log10(0.5), np.log10(500*0.73), nbin)/Hhost.Rvir

plt.figure()
rHalo,denHalo,denHaloRef,denHaloErr=Hhost.get_host_density(xbin, Rref)
plt.plot(rHalo,denHalo, 'k', label=Hhost.name+'Halo')
#for Htmp in [A2,A4]:
  #rSub,denSub,denSubRef,denSubErr=Htmp.get_sub_density((Htmp.m>fMin*Htmp.Mvir/Htmp.mP)&(Htmp.massTVV[:,iInfall]>nMinInfall), xbin, Rref)
  #plt.plot(rSub[denSub>0],denSub[denSub>0], '-', label=Htmp.name+'Sub')

rSub,denSub,denSubRef,denSubErr=Hsub.get_sub_density((Hsub.m>fMin*Hsub.Mvir/Hsub.mP)&(Hsub.massTVV[:,iInfall]>nMinInfall), xbin, Rref)
plt.errorbar(rSub,denSub, denSubErr, fmt='bo', label=Hsub.name+'Sub')
plt.yscale('log', nonposy='clip')
plt.xscale('log')

f=(denSub/denHalo>0)&(rSub<Rref)
w=1/np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)[f]
#w=None
polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],1, w=w, cov=True)
print polypar, np.sqrt(polyV.diagonal())
plt.plot(rHalo[denHalo>0], (denHalo*(rHalo)**polypar[0]*np.exp(polypar[1]))[denHalo>0],'r--', label='Model')
ylnerr=np.sqrt(np.log(rHalo)**2*polyV[0,0]+polyV[1,1]+2*np.log(rSub)*polyV[0,1])
y=(denHalo*(rHalo)**polypar[0]*np.exp(polypar[1]))
#plt.fill_between(rHalo[denHalo>0], y/np.exp(ylnerr), y*np.exp(ylnerr), color='r', alpha=0.2)
#plt.plot(rHalo, denHalo*(rHalo/Rref)**1.2,'r--', label=r'$\beta=1.2$')

plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho(R)/\rho(R_{200})$')
plt.legend()
#plt.savefig(outdir+'/'+Hsub.name+'subprof.eps')

plt.figure()  
plt.errorbar(rSub, denSub/denHalo, yerr=np.sqrt((denSubErr/denSub)**2+(denHaloErr/denHalo)**2)*denSub/denHalo, fmt='bo',label=Hsub.name+'Sub') #label=r'$m/M_{200}>10^{%d}$'%(np.log10(fMin)))  
#for Htmp in [A2,A4]:
  #rSub,denSub,denSubRef,denSubErr=Htmp.get_sub_density((Htmp.m>fMin*Htmp.Mvir/Htmp.mP)&(Htmp.massTVV[:,iInfall]>nMinInfall), xbin, Rref)
  #plt.plot(rSub[denSub>0],(denSub/denHalo)[denSub>0], '-', label=Htmp.name+'Sub')
plt.yscale('log', nonposy='clip')
plt.xscale('log')

#plt.plot(rSub, np.exp(np.polyval(polypar, np.log(rSub))),'r--', label=r'$\gamma=%.2f\pm%.2f$'%(polypar[0],np.sqrt(polyV[0,0])))
plt.plot(rSub, np.exp(np.polyval(polypar, np.log(rSub))),'r--', label=r'$\gamma=%.1f$'%(polypar[0]))
ylnerr=np.sqrt(np.log(rSub)**2*polyV[0,0]+polyV[1,1]+2*np.log(rSub)*polyV[0,1])
y=np.exp(np.polyval(polypar, np.log(rSub)))
#plt.fill_between(rSub, y/np.exp(ylnerr), y*np.exp(ylnerr), color='r', alpha=0.2)
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

#REinasto=0.64
#alpha=1.0
#denEinasto=einasto(rSub/REinasto,alpha)/einasto(Rref/REinasto, alpha)
#plt.plot(rSub[denHalo>0], (denEinasto/denHalo)[denHalo>0], 'b-',label='Einasto2')
#legend1=plt.legend(loc=2)

#REinasto=0.58
#alpha=1.0
#denEinasto=einasto(rSub/REinasto,alpha)/einasto(Rref/REinasto, alpha)
#c=[5.72,4.41,4.27,3.88,3.48,3.81,0.78,1.98,4.18]
#a=[0.215,0.235,0.181,0.205,0.149,0.186,0.097,0.117,0.190]
#denall=[]
#for i in xrange(9):
  #denHalo=einasto(rSub*c[i], a[i])/einasto(Rref*c[i], a[i])
  #denall.append(denHalo)
#denall=np.array(denall)  
#plt.plot(rSub, (denEinasto/denall.mean(0))/5, 'c-',label='Einasto3')
#legend1=plt.legend(loc=2)