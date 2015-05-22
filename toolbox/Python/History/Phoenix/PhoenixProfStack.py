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

hlist='ABCDEF'

PhHalo={}
for H in hlist:
  PhHalo[H]=HaloData('Ph'+H+'2')
PhHalo['a']=HaloData('AqA1')

Rref=1
nMinInfall=0
iInfall=0
fMin=1e-5

Rstar=3.55
Mstar=0.09
alpha=0.99
beta=0.9
sigma=1.
fs=0.56 #1????

colors='rgbcmyk'

#ratio
xbin=np.logspace(-2, np.log10(2), 30)
denHaloAll=[]
denSubAll=[]
denRatAll=[]
for i,H in enumerate(hlist):
  rHalo,denHalo,denHaloRef,denHaloErr=PhHalo[H].get_host_density(xbin, Rref)
  rSub,denSub,denSubRef,denSubErr=PhHalo[H].get_sub_density((PhHalo[H].m>fMin*PhHalo[H].Mvir/PhHalo[H].mP), xbin, Rref)
  denRatAll.append(denSub/denHalo)
  denHaloAll.append(denHalo*denHaloRef/PhHalo[H].Mvir*PhHalo[H].Rvir**3)#dP/d^3(R/Rv)
  denSubAll.append(denSub*denSubRef)
denHaloAll=np.array(denHaloAll)
denHalo=denHaloAll.mean(0)
denSubAll=np.array(denSubAll)
denSub=denSubAll.mean(0)
denRatAll=np.array(denRatAll)
denRat=denRatAll.mean(0)
plt.figure()
for i,H in enumerate(hlist):
  plt.subplot(131)
  plt.loglog(rSub, denHaloAll[i])
  plt.subplot(132)
  plt.loglog(rSub, denSubAll[i])
  plt.subplot(133)
  plt.loglog(rSub, denRatAll[i])
plt.figure()
#plt.loglog(rSub, (denSub/denHalo)/(denSub/denHalo)[20]*denRat[20], 'r')
plt.loglog(rSub[denRat>0], denRat[denRat>0], 'gs', label='Phoenix')
plt.fill_between(rSub, denRat+denRatAll.std(0), denRat-denRatAll.std(0), color='g', alpha=0.3)
plt.yscale('log', nonposy='clip')
f=(rSub<Rref)&(denRat>0)&(rSub>0.1)
pars=powerlaw_fit(rSub[f], denRat[f])
plt.plot(rSub, (rSub/pars[0][0])**pars[0][1], 'g-', label=r'$\gamma=%.1f$'%pars[0][1])
print pars
gamma=pars[0][1]
def plot_haloArat():
  H='a'
  fMin=1e-6
  rHalo,denHalo,denHaloRef,denHaloErr=PhHalo[H].get_host_density(xbin, Rref)
  rSub,denSub,denSubRef,denSubErr=PhHalo[H].get_sub_density((PhHalo[H].m>fMin*PhHalo[H].Mvir/PhHalo[H].mP), xbin, Rref)
  plt.plot(rSub[denSub>0], (denSub/denHalo)[denSub>0], 'ro', label=PhHalo[H].name)
  f=(rSub<Rref)&(denSub>0)
  pars=powerlaw_fit(rSub[f], (denSub/denHalo)[f])
  plt.plot(rSub, (rSub/pars[0][0])**pars[0][1], 'r-', label=r'$\gamma=%.1f$'%pars[0][1])
plot_haloArat()  
plt.legend(loc=4)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
#plt.savefig(outdir+'/subprof_rat.Phoenix.pdf')

#unevolve MF
plt.figure()
xbin=np.logspace(-7,-2,100)
x=xbin[:-1]*(xbin[1]/xbin[0])
NAll=[]
for h in hlist:
  H=PhHalo[h]
  m=H.massTVV[:,iInfall]
  data=m[(H.r<Rref*H.Rvir)&(m>100)]*H.mP/H.Mvir
  n,_=np.histogram(data, xbin)
  NAll.append(n)
NAll=np.array(NAll)/np.log(xbin[1]/xbin[0])
plt.loglog(x, NAll.mean(0), 's', mfc='none', mec='g', label='Phoenix')
plt.fill_between(x, np.percentile(NAll, 25, axis=0), np.percentile(NAll, 75, axis=0), color='g', alpha=0.3)
plt.yscale('log', nonposy='clip')
#for i in range(7):
  #plt.plot(x, NAll[i])
f=(NAll.sum(0)>0)&(x>1e-6)
pars=powerlaw_fit(x[f],NAll.mean(0)[f], w=1/np.sqrt(NAll.sum(0))[f])
plt.plot(x, (x/pars[0][0])**pars[0][1], 'g-', label=r'$\alpha=%.2f$'%(-pars[0][1]))
print pars
alpha=-pars[0][1]
Mstar=pars[0][0]
def plot_haloAmf():
  H=PhHalo['a']
  m=H.massTVV[:,iInfall]
  data=m[(H.r<Rref*H.Rvir)&(m>100)]*H.mP/H.Mvir
  n,_=np.histogram(data, xbin)
  plt.plot(x, n/np.log(xbin[1]/xbin[0]), 'o', mfc='none', mec='r', label=H.name)
  f=n>0  
  pars=powerlaw_fit(x[f],n[f]/np.log(xbin[1]/xbin[0]), w=1/np.sqrt(n)[f])
  plt.plot(x, (x/pars[0][0])**pars[0][1], 'r-', label=r'$\alpha=%.2f$'%(-pars[0][1]))
plot_haloAmf()
plt.legend(loc=1)
plt.xlabel(r'$m_{\rm acc}/M_{200}$')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}\ln m_{\rm acc}$')
#plt.savefig(outdir+'InfallMF.Phoenix.pdf')

beta=gamma/alpha
print alpha, beta, Mstar, Rstar
#final model
plt.figure()
plt.loglog(rSub, denSub, 'r')
rho0=fs*np.exp(sigma**2*alpha**2/2)*((fMin/Mstar)**(-alpha)-(0.1/Mstar)**(-alpha))/alpha*denHalo # dN/d^3(R/Rv)
rho=rho0*(rSub/Rstar)**(alpha*beta)
plt.plot(rSub, rho,'--', color=colors[i])
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$dN/d(R/R_{200})^3$')
#plt.savefig(outdir+'PhProfPred.eps')
