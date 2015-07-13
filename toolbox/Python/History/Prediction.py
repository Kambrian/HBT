import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'
  
Rref=1
nbin=30
Rconv=0.3
nMinInfall=0
iInfall=0

mustar=0.35
beta=1.23
A=0.073
alpha=0.96
#alpha=0.99; A=0.056 #weighted fit to MF
sigma=1.2
fs=0.56 #this 1.2 fudge factor is not good....
B=0.6

Halo=HaloData('AqA1')

xbin=np.logspace(np.log10(1), np.log10(500*0.73), nbin)/Halo.Rvir

plt.figure()
rHalo,denHalo,denHaloRef,denHaloErr=Halo.get_host_density(xbin, Rref)
for fMin,fmt in [(1e-7,'o'),(1e-6,'d'),(1e-5,'s'),(1e-4,'^')]:
  mMin=fMin*Halo.Mvir
  rSub,denSub,denSubRef,denSubErr=Halo.get_sub_density((Halo.m>mMin/Halo.mP)&(Halo.massTVV[:,iInfall]>nMinInfall), xbin, Rref)
  #plt.errorbar(rSub, denSub*denSubRef, denSubErr*denSubRef, fmt='bo', label=Halo.name+'Sub')
  l,=plt.plot(rSub, denSub*denSubRef, fmt, label=r'$m/M_{200}>%s$'%fmtexp10(fMin))
  plt.yscale('log', nonposy='mask')
  plt.xscale('log')
  rho0=A*B*Halo.Mvir*fs*np.exp(sigma**2*alpha**2/2.)*(mMin**(-alpha)-(0.01*Halo.Mvir)**(-alpha))/alpha*denHalo*denHaloRef/Halo.Mvir*Halo.Rvir**3 # dN/d^3(R/Rv)
  rho=rho0*(rSub)**(alpha*beta)*mustar**alpha
  plt.plot(rSub, rho,'--', color=l.get_color())
  print A*(mMin**(-alpha)-(0.1*Halo.Mvir)**(-alpha))/alpha
  x=rSub
  f=x<1
  print 'Mint/M200=', np.sum((denHalo*denHaloRef*Halo.Rvir**3*x**3)[f])*4*np.pi*np.log(x[1]/x[0])/Halo.Mvir
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
#plt.savefig(outdir+'/'+Halo.name+'subprof_rat.eps')