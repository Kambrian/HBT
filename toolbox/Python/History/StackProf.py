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

HList='ABCDEF'
AqHalo=[HaloData('Aq'+H+'2') for H in HList]
PhHalo=[HaloData('Ph'+H+'2') for H in HList]

Rref=1

sigma=1.1
B=0.65

#ratio
xbin=np.logspace(-2, np.log10(2),100)
def stackHalo(HList, pars, fMin=1e-5, color='r', fmt='-', sym='o'):
  denHaloAll=[]
  denSubAll=[]
  denRatAll=[]
  Mh=np.mean([h.Mvir for h in HList])
  mMin=fMin*Mh
  for i,H in enumerate(HList):
	rHalo,denHalo,denHaloRef,denHaloErr=H.get_host_density(xbin, Rref)
	rSub,denSub,denSubRef,denSubErr=H.get_sub_density((H.m>mMin/H.mP), xbin, Rref)
	denRatAll.append(denSub/denHalo)
	denHaloAll.append(denHalo*denHaloRef*H.Rvir**3)#dM/d^3(R/Rv)
	denSubAll.append(denSub*denSubRef)
  denHaloAll=np.array(denHaloAll)
  denHalo=denHaloAll.mean(0)
  denSubAll=np.array(denSubAll)
  denSub=denSubAll.mean(0)
  denRatAll=np.array(denRatAll)
  denRat=denRatAll.mean(0)
  #return rHalo,denHaloAll,denSubAll,denRatAll,Mh
  x=rSub
  l,=plt.loglog(x[denSub>0], denSub[denSub>0], sym, color=color, label=r'$m/M_{200}>%s$'%fmtexp10(fMin))
  fs,A,alpha,mustar,beta=pars
  rho0=A*B*fs*np.exp(sigma**2*alpha**2/2.)*(mMin**(-alpha)-(0.01*Mh)**(-alpha))/alpha*denHalo# dN/d^3(R/Rv)
  rho=rho0*(mustar*x**beta)**alpha
  plt.plot(x, rho, fmt, color=color)
  print 'factor=', np.sum(denHalo*x**3*(mustar*x**beta)**alpha)/np.sum(denHalo*x**3)
  return l

#Aq=stackHalo(AqHalo)
#Ph=stackHalo(PhHalo)
AqPars=[0.54, 0.089,0.95,0.42,1.37]
PhPars=[0.56, 0.08,0.94,0.34,1.0]
AvPars=[0.55, 0.084, 0.95, 0.38, 1.2]
plt.figure()

for fMin,color in [(1e-6, 'r'), (1e-5, 'g'), (1e-4, 'b')]:
  l=stackHalo(AqHalo, AqPars, fMin, 'r', '--', '-')
  l2=stackHalo(PhHalo, PhPars, fMin, 'g', '--', '-')
#legend1=plt.legend(lall, [h.get_label() for h in lall], loc=1)
plt.legend([l,l2], ['Aquarius','Phoenix'])
#plt.gca().add_artist(legend1)
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$dN/d(R/R_{200})^3$')
#plt.savefig(outdir+'ProfPred.All.eps')

#different stripping functions
plt.figure()
x=np.logspace(-1,0,20)
#AqA1
plt.loglog(x, 0.35*x**1.23)
#AqA to F
plt.loglog(x, 0.34*x**1.37, 'r')
plt.loglog(x, 0.61*x**1.59, 'r')
plt.loglog(x, 0.56*x**1.95, 'r')
plt.loglog(x, 0.42*x**1.40, 'r')
plt.loglog(x, 0.50*x**1.30, 'r')
plt.loglog(x, 0.31*x**1.07, 'r')

#PhF to A
plt.loglog(x, 0.34*x**1.09, 'g')
plt.loglog(x, 0.31*x**0.90, 'g')
plt.loglog(x, 0.39*x**0.95, 'g')
plt.loglog(x, 0.41*x**1.31, 'g')
plt.loglog(x, 0.48*x**1.13, 'g')
plt.loglog(x, 0.25*x**0.88, 'g')

plt.loglog(x, 0.43*x**1.37, 'r--')
plt.loglog(x, 0.34*x**1.0, 'g--')
muAll=[0.34, 0.61, 0.56, 0.42, 0.50, 0.31, 0.25, 0.48, 0.41, 0.39, 0.31, 0.34]
betaAll=[1.37, 1.59, 1.95, 1.40, 1.30, 1.07, 0.88, 1.13, 1.31, 0.95, 0.90, 1.09]
