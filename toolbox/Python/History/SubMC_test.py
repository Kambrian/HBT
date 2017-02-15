''' generate Monte-Carlo samples of subhaloes '''
import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.optimize import fsolve,curve_fit
from scipy.integrate import quad
from nfw import NFWHalo,mean_concentration
from MbdIO import HaloData

plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

HaloA=HaloData('AqA2')
HaloA.spline=interp1d(np.log(HaloA.DM[:,0]), np.log(HaloA.DM[:,1]), bounds_error=False)

MASS_HOST=HaloA.Mvir
#CONC_HOST=5
NUM_SUBS=1e5
C_OFFSET=0.

MASS_MIN_INFALL=1e-6*HaloA.Mvir #in units of 1e10Msun/h
MASS_MAX_INFALL=MASS_HOST/10. #MASS_HOST/10. #
R_MIN=0
R_MAX=2.

sigma=1.1
#B=0.6
fs=0.55
A=0.1*MASS_HOST**-0.02
alpha=0.95
mustar=0.5*MASS_HOST**-0.03
beta=1.7*MASS_HOST**-0.04

Host=NFWHalo(MASS_HOST)#,Chost)
Host.density=lambda R: np.exp(HaloA.spline(np.log(R))) #interp in log-space

def lnPDF(x):
  ''' R in units of Rv'''
  lnm, lnm0, lnc0, lnR=x
  if lnm0>np.log(MASS_MAX_INFALL):
	return -np.inf
  if lnm0<np.log(MASS_MIN_INFALL): #lower bound 1e-3
	return -np.inf
  lnmubar=np.log(mustar)+beta*lnR
  lnmu=lnm-lnm0-lnmubar
  if lnm-lnm0>0:
	return -np.inf
  if lnmu>np.log(4.2):
	return -np.inf
  if lnR>np.log(R_MAX): #log(3)=1.1
	return -np.inf
  
  lnPDFm=-0.5*(lnmu/sigma)**2
  lnPDFm0=-alpha*lnm0
  lncbar=np.log(mean_concentration(np.exp(lnm0),'Duffy'))
  lnPDFc0=-0.5*((lnc0-lncbar)/0.35+C_OFFSET)**2 #log10-scatter=0.15, ln-scatter=0.35 #shift by 1sigma downward
  lnPDFR=3.*lnR+np.log(Host.density(np.exp(lnR)*Host.Rv))
  return lnPDFm+lnPDFm0+lnPDFc0+lnPDFR

#os.environ['OMP_NUM_THREADS']='32'
import emcee
from numpy import *

Host.L=Host.luminosity(Host.Rv)
Host.Msample=Host.mass(R_MAX*Host.Rv)-Host.mass(R_MIN*Host.Rv)
print 'Fraction of mass sampled inside Rvir:', (Host.mass(Host.Rv)-Host.mass(R_MIN*Host.Rv))/Host.M
NUM_SUBS_PRED=A*(MASS_MIN_INFALL**-alpha-MASS_MAX_INFALL**-alpha)/alpha*Host.Msample #infall number scale with host mass (in the given radial range: Mv*M/Mv=M)
NUM_HOSTS=NUM_SUBS/NUM_SUBS_PRED/fs #only sample survived subhaloes.
print 'sample size: %.1g hosts'%NUM_HOSTS

nwalkers=10
nburn=300
nsteps=int(NUM_SUBS/nwalkers+nburn)
print 'running %d steps'%nsteps
ndim=4

x00=array([np.log(MASS_MIN_INFALL*10.),np.log(MASS_MIN_INFALL*20.),np.log(mean_concentration(MASS_MIN_INFALL*20.,'Duffy')),-0.5])
labels=[r"$\ln m$",r"$\ln m_\mathrm{acc}$",r"$\ln c_{\rm acc}$",r"$\ln R/R_{200}$"]
x0=kron(ones([nwalkers,1]),x00)#repmat to nwalkers rows
x0+=random.rand(ndim*nwalkers).reshape(nwalkers,ndim)-0.5 #random offset, [-0.5,0.5]
sampler=emcee.EnsembleSampler(nwalkers,ndim,lnPDF)

sampler.run_mcmc(x0,nsteps)

from matplotlib.pyplot import *
ion()
figure()
for i in range(ndim):
  subplot(ndim,1,i)
  for j in range(nwalkers):
    plot(range(nsteps),sampler.chain[j,:,i],'.')
  ylabel(labels[i])
xlabel('Step')  

sample=sampler.chain[:,nburn:,:]
flatchain=sample.reshape([-1,ndim])
flatchain=np.exp(flatchain) #from log() to linear
nsub=flatchain.shape[0]
nsub_disrupt=nsub/fs*(1-fs) #number of disrupted subhaloes.
print '%d disrupted subhaloes not sampled'%nsub_disrupt

m,mInfall,cInfall,R=flatchain.T
plt.figure()
plt.subplot(121)
y,x=np.histogram(np.log(R),30)
y/=(x[1]-x[0])
x=np.exp(x[1:])
plt.loglog(x,y/4/np.pi/x**3/NUM_HOSTS)
plt.plot(x, Host.density(x*Host.Rv)/Host.Msample*Host.Rv**3*NUM_SUBS_PRED*fs)
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$\mathrm{d}P/\mathrm{d}^3(R/R_{200})$')
y,x=np.histogram(np.log(R[m>MASS_MIN_INFALL]),30)
y/=(x[1]-x[0])
x=np.exp(x[1:])
y=y/4/np.pi/x**3/NUM_HOSTS
plt.loglog(x,y)
rSub,denSub,denSubRef,denSubErr=HaloA.get_sub_density((HaloA.m>1e-6*HaloA.Mvir/HaloA.mP), x, 1)
plt.errorbar(rSub,denSub*denSubRef, denSubErr*denSubRef, fmt='bo', label=HaloA.name+'Sub')

plt.figure()
plt.subplot(121)
y,x=np.histogram(np.log(R),30)
y/=(x[1]-x[0])
x=np.exp(x[1:])
plt.loglog(x,y/4/np.pi/x**3/np.sum(R<1))
plt.plot(x, Host.density(x*Host.Rv)/Host.M*Host.Rv**3)
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$\mathrm{d}P/\mathrm{d}^3(R/R_{200})$')
y,x=np.histogram(np.log(R[m>MASS_MIN_INFALL]),30)
y/=(x[1]-x[0])
x=np.exp(x[1:])
y=y/4/np.pi/x**3/np.sum(R[m>MASS_MIN_INFALL]<1)
plt.loglog(x,y)
rSub,denSub,denSubRef,denSubErr=HaloA.get_sub_density((HaloA.m>1000), x, 1)
plt.errorbar(rSub,denSub/denSub[10]*y[10], denSubErr/denSub[10]*y[10], fmt='bo', label=HaloA.name+'Sub')
plt.subplot(122)
y=y/(Host.density(x*Host.Rv)/Host.M*Host.Rv**3)
plt.loglog(x, y, 'o')
pars=np.polyfit(np.log(x[y>0]), np.log(y[y>0]), 1)
plt.plot(x, np.exp(pars[1])*x**(alpha*beta), '-')
plt.legend(loc=2)
  
plt.figure()
y,x=np.histogram(np.log(m[(m>0)&(R<1)]),100)
y/=(x[1]-x[0])
x=np.exp(x[1:])
FracCorrection=0.9 #scipy.stats.norm.cdf(np.log(4.2)/sigma), to correct for the normalization due to mu_max cut-off
plt.loglog(x,y*x/NUM_HOSTS*FracCorrection, 'g', label='current')
plt.plot(x, x**-alpha*0.008*Host.M*x, 'g')

y,x=np.histogram(np.log(mInfall[R<1]),50)
f=y>10
y/=(x[1]-x[0])
y/=fs #account for disrupted part.
x=np.exp(x[1:])
plt.loglog(x,y*x/NUM_HOSTS, 'r', label='infall')
plt.plot(x, x**-alpha*A*Host.M*x, 'r')
#f=(y>0)&(m<1e-6)
print 'mass func:', np.polyfit(np.log(x)[f],np.log(y)[f],1)
plt.xlabel(r'$m$')
plt.ylabel(r'$m\mathrm{d}N/\mathrm{d}\ln m$')

