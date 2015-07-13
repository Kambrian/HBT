import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.optimize import fsolve
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
B=0.6

def virial_density(virtype='c200', scaleF=1., Omega0=0.25):
  '''return virial density'''
  OmegaL=1.-Omega0
  G=43007.1
  HUBBLE0=0.1
  Hz=HUBBLE0 * np.sqrt(Omega0 /scaleF**3+ (1. -Omega0 -OmegaL) / scaleF**2 +OmegaL);
  Hratio=Hz/HUBBLE0
  OmegaZ=Omega0/scaleF**3/Hratio**2

  virialF={'tophat': 18.0*np.pi**2+82.0*(OmegaZ-1)-39.0*(OmegaZ-1)**2,
		  'c200': 200,
		  'b200': 200*OmegaZ}[virtype]

  RhoCrit=3*Hz**2/8./np.pi/G
  #print virialF
  return RhoCrit*virialF

#ratio
xbin=np.logspace(-3, np.log10(3),100)
def AFactor(HList, pars, fMin=1e-5):
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
  x=rSub
  fs,A,alpha,mustar,beta=pars
  f=x<1
  #print 'factor=', np.sum((denHalo*x**3*(mustar*x**beta)**alpha)[f])/np.sum((denHalo*x**3)[f])
  Arat=B*fs*np.exp(sigma**2*alpha**2/2.)*np.sum((denHalo*x**3*(mustar*x**beta)**alpha)[f])/np.sum((denHalo*x**3)[f])
  print 'A/A_acc=', Arat, '; A=', A*Arat
  f=~np.isnan(denHalo)
  FuncHalo=UnivariateSpline(np.log(x)[f],(denHalo*x**3)[f],s=0, ext=0)
  FuncHaloStrip=UnivariateSpline(np.log(x)[f], (denHalo*x**3*(mustar*x**beta)**alpha)[f], s=0, ext=0)
  CumDens=lambda r: FuncHalo.integral(-100., np.log(r))/(r**3)/FuncHalo.integral(-100., 0)*virial_density()
  rB200=fsolve(lambda r: CumDens(r)-virial_density('b200'), 1.)
  rVir=fsolve(lambda r: CumDens(r)-virial_density('tophat'), 1.)
  y=[FuncHaloStrip.integral(-100., a)/FuncHalo.integral(-100.,a) for a in np.log(x)[1:]]
  y=np.array(y)*A*B*fs*np.exp(sigma**2*alpha**2/2.)
  denUp=denHalo*x**3*(mustar*x**beta)**alpha
  denDown=denHalo*x**3
  yy=[np.sum(denUp[:i])/np.sum(denDown[:i]) for i in xrange(1,len(x))]
  yy=np.array(yy)*A*B*fs*np.exp(sigma**2*alpha**2/2.)
  print "Two integrals:", FuncHalo.integral(-10., 1), np.sum(denDown[x<1])*np.log(x[1]/x[0])
  return x[1:],y,yy,rVir,rB200

#Aq=stackHalo(AqHalo)
#Ph=stackHalo(PhHalo)
AqPars=[0.54, 0.089,0.95,0.42,1.37]
PhPars=[0.56, 0.08,0.95,0.34,1.0]
AvPars=[0.55, 0.084, 0.95, 0.38, 1.2]
plt.figure()

x0,y0,y00,RvVir,RvB200=AFactor(AqHalo, AqPars)
x1,y1,y11,_,_=AFactor(PhHalo, PhPars)
#plt.plot(x0,y0,'r-')
#plt.plot(x1,y1,'g-')
plt.plot(x0,y00,'r-')
plt.plot(x1,y11,'g-')
plt.xlim([0.05,3])
plt.ylim([1e-4,1e-1])
plt.loglog()
plt.plot([RvVir, RvVir], plt.ylim(), 'k--')
plt.plot([RvB200, RvB200], plt.ylim(), 'k:')
plt.plot([1,1], plt.ylim(), 'k-')
#plt.plot([Rv1, Rv1], plt.ylim(), 'g--')
#legend1=plt.legend(lall, [h.get_label() for h in lall], loc=1)
plt.legend(['Aquarius','Phoenix'], loc=2)
#plt.gca().add_artist(legend1)
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$A$')
plt.loglog()
plt.savefig(outdir+'MFAmplitude.eps')