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
#A3=HaloData('A3')
#A4=HaloData('A4')
#A5=HaloData('A5')
#A2H=HaloData('A2', None)
#A4H=HaloData('A4', None)

##plot profiles
Rref=1
nbin=30
Rconv=0.3
nMinInfall=0
iInfall=0

Rstar=2.07
Mstar=0.062
alpha=0.99
beta=1.43
sigma=1.
fs=0.56 #1????

Hhost=A2
Hsub=A1

xbin=np.logspace(np.log10(1), np.log10(500*0.73), nbin)/Hhost.Rvir

plt.figure()
rHalo,denHalo,denHaloRef,denHaloErr=Hhost.get_host_density(xbin, Rref)
for fMin,fmt in [(1e-7,'o'),(1e-6,'d'),(1e-5,'s'),(1e-4,'^')]:
  rSub,denSub,denSubRef,denSubErr=Hsub.get_sub_density((Hsub.m>fMin*Hsub.Mvir/Hsub.mP)&(Hsub.massTVV[:,iInfall]>nMinInfall), xbin, Rref)
  #plt.errorbar(rSub, denSub*denSubRef, denSubErr*denSubRef, fmt='bo', label=Hsub.name+'Sub')
  l,=plt.plot(rSub, denSub*denSubRef, fmt, label=r'$m_0/M_{200}>%s$'%fmtexp10(fMin))
  plt.yscale('log', nonposy='mask')
  plt.xscale('log')
  rho0=fs*np.exp(sigma**2*alpha**2/2)*((fMin/Mstar)**(-alpha)-(0.1/Mstar)**(-alpha))/alpha*denHalo*denHaloRef/Hhost.Mvir*Hhost.Rvir**3 # dN/d^3(R/Rv)
  rho=rho0*(rSub/Rstar)**(alpha*beta)
  plt.plot(rSub, rho,'--', color=l.get_color())
  print ((fMin/Mstar)**(-alpha)-(0.1/Mstar)**(-alpha))/alpha
plt.ylim([1,1e6])  
#plt.plot(rHalo, denHalo*denSubRef, 'g-', label='Halo')
plt.legend(loc=3, fontsize=15)
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$dN/d(R/R_{200})^3$')
#plt.savefig(outdir+'AqProfPred.eps')
#nbin=30
#nMinInfall=2000
#H=A1
#f=(H.massTVV[:,iInfall]>nMinInfall)#&(H.massTVV[:,iInfall]<nMinInfall*10)
#x=H.r[f]/H.Rvir
#y=(1.0*H.m/H.massTVV[:,iInfall]**1.)[f]
#countAll,_=np.histogram(x, xbin)
#p,xmid=cross_section(x,y,xbin,40+0.6*np.array([(100-68.3)/2,50,(100+68.3)/2]))
#fall=[xmid, p[1]]
#p,xmid=cross_section(x[y>0],y[y>0],xbin,[(100-68.3)/2,50,(100+68.3)/2])
#fres=[xmid, p[1]]
#plt.plot(rHalo, rho0*fres[1]**alpha, 'b--')
#plt.plot(rHalo, rho0*fall[1]**alpha, 'c--')
#plt.savefig(outdir+'/'+Hsub.name+'subprof_rat.eps')