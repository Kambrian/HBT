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

Halo={}
HList=['Aq'+h+'2' for h in 'ABCDEF']+['Ph'+h+'2' for h in 'ABCDEF']

for H in HList:
  Halo[H]=HaloData(H)
MAll=[Halo[h].Mvir for h in HList]
cAll=[Halo[h].Cvir for h in HList]

muAll=np.array([0.34, 0.61, 0.56, 0.42, 0.50, 0.31, 0.25, 0.48, 0.41, 0.39, 0.31, 0.34])
betaAll=np.array([1.37, 1.59, 1.95, 1.40, 1.30, 1.07, 0.88, 1.13, 1.31, 0.95, 0.90, 1.09])

#mu=[0.42, 0.34]
#beta=[1.4, 1.0]
#M=[np.mean([HaloData('Aq'+h+'2').Mvir for h in 'ABCDEF']), np.mean([HaloData('Ph'+h+'2').Mvir for h in 'ABCDEF'])]
#R=[np.mean([HaloData('Aq'+h+'2').Rvir for h in 'ABCDEF']), np.mean([HaloData('Ph'+h+'2').Rvir for h in 'ABCDEF'])]
#c=[np.mean([HaloData('Aq'+h+'2').Cvir for h in 'ABCDEF']), np.mean([HaloData('Ph'+h+'2').Cvir for h in 'ABCDEF'])]
#x=np.logspace(-2,0,10)
#plt.loglog(x, mu[0]*x**beta[0], 'r')
#plt.loglog(x, mu[1]*x**beta[1], 'g')

Rref=1
nMinInfall=0
iInfall=0
fMin=1e-5

colors='rgbcmyk'

#ratio
xbin=np.logspace(-2, np.log10(2), 30)
denHaloAll=[]
denSubAll=[]
denRatAll=[]
parAll=[]
for i,H in enumerate(HList):
  rHalo,denHalo,denHaloRef,denHaloErr=Halo[H].get_host_density(xbin, Rref)
  rSub,denSub,denSubRef,denSubErr=Halo[H].get_sub_density((Halo[H].m>fMin*Halo[H].Mvir/Halo[H].mP), xbin, Rref)
  denRatAll.append(denSub/denHalo)
  denHaloAll.append(denHalo*denHaloRef/Halo[H].Mvir*Halo[H].Rvir**3)#dP/d^3(R/Rv)
  denSubAll.append(denSub*denSubRef)
  f=(rSub<Rref)&(denSub>0)&(rSub>0.1)
  pars=powerlaw_fit(rSub[f], (denSub/denHalo)[f])
  parAll.append([pars[0][1], pars[1][1]])
  
denHaloAll=np.array(denHaloAll)
denHalo=denHaloAll.mean(0)
denSubAll=np.array(denSubAll)
denSub=denSubAll.mean(0)
denRatAll=np.array(denRatAll)
parAll=np.array(parAll)
MAll=np.array(MAll)
cAll=np.array(cAll)
plt.figure()
f=MAll<1000
plt.errorbar(MAll[f], parAll[f,0], parAll[f,1], fmt='ro')
f=MAll>1000
plt.errorbar(MAll[f], parAll[f,0], parAll[f,1], fmt='go')
plt.xscale('log')
plt.xlabel(r'$M_{200}[10^{10}M_\odot/h]$')
plt.ylabel(r'$\gamma$')
#plt.savefig(outdir+'/Gamma-Mass.eps')

plt.figure()
f=MAll<1000
plt.errorbar(MAll[f], parAll[f,0], parAll[f,1], fmt='ro')
f=MAll>1000
plt.errorbar(MAll[f], parAll[f,0], parAll[f,1], fmt='go')
plt.xscale('log')
plt.xlabel(r'$M_{200}[10^{10}M_\odot/h]$')
plt.ylabel(r'$\gamma$')
x=np.logspace(-1,6)
ff=lambda x: 1.8*x**-0.05 #fit to average
plt.plot(x, ff(x)**0.95, 'k')
fff=lambda x: 2.32*x**-0.069 #fit to data points
plt.plot(x, fff(x)**0.95, 'b')
a=np.diff([1.0,1.4])/np.diff(np.log([6.7e4,100]))
b=1.4-a*np.log(100.) #fit to average
plt.plot(x, (a*np.log(x)+b)**0.95, 'c')

plt.figure()
f=MAll<1000
plt.plot(MAll[f], muAll[f], 'ro')
f=MAll>1000
plt.plot(MAll[f], muAll[f], 'go')
plt.xscale('log')
plt.xlabel(r'$M$')
plt.ylabel(r'$\mu_\star$')

plt.figure()
f=MAll<1000
plt.plot(cAll[f], muAll[f], 'ro')
f=MAll>1000
plt.plot(cAll[f], muAll[f], 'go')
plt.xscale('log')
plt.xlabel(r'$c$')
plt.ylabel(r'$\mu_\star$')
#plt.savefig(outdir+'/Gamma-c.eps')

#plt.figure()
#f=MAll<1000
#plt.errorbar(MAll[f], betaAll[f], parAll[f,1], fmt='ro')
#f=MAll>1000
#plt.errorbar(MAll[f], betaAll[f], parAll[f,1], fmt='go')
#plt.xscale('log')
#plt.xlabel(r'$M$')
#plt.ylabel(r'$\beta$')

#plt.figure()
#f=MAll<1000
#plt.errorbar(cAll[f], betaAll[f], parAll[f,1], fmt='ro')
#f=MAll>1000
#plt.errorbar(cAll[f], betaAll[f], parAll[f,1], fmt='go')
#plt.xscale('log')
#plt.xlabel(r'$c$')
#plt.ylabel(r'$\beta$')


plt.figure()
f=MAll<1000
plt.errorbar(cAll[f], parAll[f,0], parAll[f,1], fmt='ro')
f=MAll>1000
plt.errorbar(cAll[f], parAll[f,0], parAll[f,1], fmt='go')
plt.xscale('log')
plt.xlabel(r'$c$')
plt.ylabel(r'$\gamma$')
#plt.savefig(outdir+'/Gamma-c.eps')

plt.figure()
plt.loglog(rSub, denRatAll[0], 'r-', label='A')
plt.loglog(rSub, denRatAll[2], 'c-', label='C')
plt.loglog(rSub, denRatAll[1], 'b--', label='B')
plt.loglog(rSub, denRatAll[3], 'g--', label='D')
plt.loglog(rSub, denRatAll[4], 'm--', label='E')
plt.loglog(rSub, denRatAll[5], 'y--', label='F')
plt.legend(loc=2)

plt.figure()
for i,H in enumerate(HList):
  plt.subplot(131)
  plt.loglog(rSub, denHaloAll[i])
  plt.subplot(132)
  plt.loglog(rSub, denSubAll[i])
  plt.subplot(133)
  plt.loglog(rSub, denRatAll[i])
plt.figure()
#plt.loglog(rSub, (denSub/denHalo)/(denSub/denHalo)[20]*denRat[20], 'r')
denRat=denRatAll[:6].mean(0)
plt.loglog(rSub[denRat>0], denRat[denRat>0], 'ro', label='Aquarius')
plt.fill_between(rSub, denRat+denRatAll[:6].std(0), denRat-denRatAll[:6].std(0), color='r', alpha=0.3)
f=(rSub<Rref)&(denRat>0)&(rSub>0.01)
pars=powerlaw_fit(rSub[f], denRat[f])
plt.plot(rSub, (rSub/pars[0][0])**pars[0][1], 'r-', label=r'$\gamma=%.1f$'%pars[0][1])
denRat=denRatAll[6:].mean(0)
plt.loglog(rSub[denRat>0], denRat[denRat>0], 'go', label='Phoenix')
plt.fill_between(rSub, denRat+denRatAll[6:].std(0), denRat-denRatAll[6:].std(0), color='g', alpha=0.3)
plt.yscale('log', nonposy='clip')
f=(rSub<Rref)&(denRat>0)&(rSub>0.01)
pars=powerlaw_fit(rSub[f], denRat[f])
plt.plot(rSub, (rSub/pars[0][0])**pars[0][1], 'g-', label=r'$\gamma=%.1f$'%pars[0][1])
print pars
gamma=pars[0][1]
plt.legend(loc=4)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho_{\rm Sub}/\rho_{\rm Halo}$')
#plt.savefig(outdir+'/subprof_rat.Aq.pdf')
plt.figure()
for i in xrange(len(HList)):
  plt.loglog(MAll[i], cAll[i], 'o', markersize=parAll[i][0]*10)