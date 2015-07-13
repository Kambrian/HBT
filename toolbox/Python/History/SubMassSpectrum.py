import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.optimize import fsolve,curve_fit
from scipy.integrate import quad
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

lognfw=lambda r,rhos,rs: np.log(nfw(r,rhos,rs))
nfw=lambda r,rhos,rs: rhos/(r/rs)/(1+r/rs)**2
#ratio
xbin=np.logspace(-3, np.log10(3),100)

def FitHost(HList):
  denHaloAll=[]
  for i,H in enumerate(HList):
	rHalo,denHalo,denHaloRef,denHaloErr=H.get_host_density(xbin, Rref)
	denHaloAll.append(denHalo*denHaloRef*H.Rvir**3)#dM/d^3(R/Rv)
  denHaloAll=np.array(denHaloAll)
  denHalo=denHaloAll.mean(0)
  x=rHalo
  plt.figure()
  plt.loglog(x, denHalo,'ro')
  f=~np.isnan(denHalo)
  pars=curve_fit(lognfw, x[f], np.log(denHalo)[f])[0]
  plt.plot(x, np.exp(lognfw(x, *pars)),'r-')
  print pars
  return pars

a=FitHost(AqHalo)
b=FitHost(PhHalo)
def PDFSurv(x, mu, mustar, beta, sigma):
  mubar=mustar*x**beta
  p=0.
  if (mu<4.2*mubar):
	p=np.exp(-0.5*(np.log(mu/mubar)/sigma)**2)/np.sqrt(2*np.pi)/sigma
  #p[mu>4.2*mubar]=0.
  #try:
	#p[mu>4.2*mubar]=0.
  #except:
	#if mu>4.2*mubar:
	  #p=0.
  return p

def mass_spec(mu,nfwpar, pars):
  fs,A,alpha,mustar,beta=pars
  scale=quad(lambda x: nfw(x, *nfwpar)*x**2, 0,1)[0]
  return quad(lambda x: nfw(x, *nfwpar)*PDFSurv(x,mu,mustar,beta,sigma)*x**2, 0,1)[0]/scale*fs

def cum_mass_spec(mu, nfwpar, pars):
  return quad(lambda x: mass_spec(np.exp(x), nfwpar, pars), np.log(mu), 0)[0]

def mass_spec2(mu,HList, pars):
  denHaloAll=[]
  for i,H in enumerate(HList):
	rHalo,denHalo,denHaloRef,denHaloErr=H.get_host_density(xbin, Rref)
	denHaloAll.append(denHalo*denHaloRef*H.Rvir**3)#dM/d^3(R/Rv)
  denHaloAll=np.array(denHaloAll)
  denHalo=denHaloAll.mean(0)
  r=rHalo
  f=~np.isnan(denHalo)
  FuncHalo=UnivariateSpline(r[f],denHalo[f],s=0, ext=0)
  fs,A,alpha,mustar,beta=pars
  scale=quad(lambda x: FuncHalo(x)*x**2, 0,1)[0]
  return quad(lambda x: FuncHalo(x)*PDFSurv(x, mu, mustar, beta, sigma)*x**2, 0,1)[0]/scale*fs


AqPars=[0.54, 0.089,0.95,0.42,1.37]
PhPars=[0.56, 0.08,0.95,0.34,1.0]
AvPars=[0.55, 0.084, 0.95, 0.38, 1.2]
plt.figure()
mu=np.logspace(-3,0,20)
PAq=np.array([mass_spec(x, a, AqPars) for x in mu])
PPh=np.array([mass_spec(x, b, PhPars) for x in mu])
#PAq=[mass_spec2(x, AqHalo, AqPars) for x in mu]
#PPh=[mass_spec2(x, PhHalo, PhPars) for x in mu]
plt.semilogx(mu, PAq*np.log(10.), 'r-')
plt.semilogx(mu, PPh*np.log(10.), 'g-')

plt.legend(['Aquarius','Phoenix'], loc=4)
plt.xlabel(r'$m/m_{\rm acc}$')
plt.ylabel(r'$\mathrm{d}P/\mathrm{d}\log10(m)$')


def SubMassSpecData(HList, nMin=1e4, iInfall=0):
  muAll=[]
  for i,H in enumerate(HList):
	f=(H.massTVV[:,iInfall]>nMin)&(H.r<H.Rvir)
	mu=H.m[f]*1./H.massTVV[f,iInfall]
	muAll.extend(mu)
  muAll=np.array(muAll)
  #y,x=np.histogram(np.log(muAll), np.linspace(np.log(20./nMin),0))
  #y=np.cumsum(y[::-1])[::-1]
  #print np.sum(muAll==0)*1./len(muAll) #has resolution effect!!
  #return x,y*1./len(muAll)
  y,x=np.histogram(np.log(muAll), np.linspace(np.log(20./nMin),0,30), density=True)
  return x,y*np.sum(muAll>20./nMin)*1./len(muAll)

x,y=SubMassSpecData(AqHalo)
#plt.plot(np.exp(x[1:]), y*np.log(10.), 'r:')
x,y=SubMassSpecData(PhHalo)
#plt.plot(np.exp(x[1:]), y*np.log(10.), 'g:')

HaloA=[HaloData('AqA%s'%i) for i in '12345']
def SubMassSpecData2(HList, fMin=1e-4, iInfall=0):
  for i,H in enumerate(HList):
	nMin=H.Mvir*fMin/H.mP
	f=(H.massTVV[:,iInfall]>nMin)&(H.r<H.Rvir)
	mu=H.m[f]*1./H.massTVV[f,iInfall]
	y,x=np.histogram(np.log(mu), np.linspace(np.log(20./nMin),0,20), density=True)
	plt.plot(np.exp(x)[1:], y*np.sum(mu>20./nMin)*1./len(mu),':')
	#y,x=np.histogram(np.log(mu), np.linspace(np.log(20./nMin),0))
	#y=np.cumsum(y[::-1])[::-1]
	#plt.plot(np.exp(x)[1:],y*1./len(mu),':')

#SubMassSpecData2(HaloA)
plt.loglog()
plt.savefig(outdir+'SubFate0.eps')