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

HaloA=HaloData('PhA2')
HaloA.spline=interp1d(np.log(HaloA.DM[:,0]), np.log(HaloA.DM[:,1]), bounds_error=False)

MASS_HOST=6.7e4
#CONC_HOST=5
NUM_SUBS=float(sys.argv[1]) #1e5
C_OFFSET=0.

MASS_MIN_INFALL=float(sys.argv[2]) #in units of 1e10Msun/h
MASS_MAX_INFALL=MASS_HOST/10. #min(MASS_MIN_INFALL*1e4, MASS_HOST/10.)
R_MIN=1e-3
R_MAX=2.

sigma=1.1
#B=0.6
fs=0.55
A=0.1*MASS_HOST**-0.02
alpha=0.95
mustar=0.5*MASS_HOST**-0.03
beta=1.7*MASS_HOST**-0.04

Host=NFWHalo(MASS_HOST)#,Chost)
#Host.density=lambda R: np.exp(HaloA.spline(np.log(R))) #interp in log-space

def lnPDF(x):
  ''' R in units of Rv'''
  lnm, lnm0, lnR=x
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
  if lnR<np.log(R_MIN): #log(1e-3)
	return -np.inf
  
  lnPDFm=-0.5*(lnmu/sigma)**2
  lnPDFm0=-alpha*lnm0
  lnPDFR=3.*lnR+np.log(Host.density(np.exp(lnR)*Host.Rv))
  return lnPDFm+lnPDFm0+lnPDFR

#os.environ['OMP_NUM_THREADS']='32'
import emcee
from numpy import *

Host.Msample=Host.mass(R_MAX*Host.Rv)-Host.mass(R_MIN*Host.Rv)
print 'Fraction of mass sampled inside Rvir:', (Host.mass(Host.Rv)-Host.mass(R_MIN*Host.Rv))/Host.M
NUM_SUBS_PRED=A*(MASS_MIN_INFALL**-alpha-MASS_MAX_INFALL**-alpha)/alpha*Host.Msample #infall number scale with host mass (in the given radial range: Mv*M/Mv=M)
NUM_HOSTS=NUM_SUBS/NUM_SUBS_PRED/fs #only sample survived subhaloes.
print 'sample size: %.1g hosts'%NUM_HOSTS

nwalkers=10
nburn=300
nsteps=int(NUM_SUBS/nwalkers+nburn)
print 'running %d steps'%nsteps
ndim=3

x00=array([np.log(MASS_MIN_INFALL*10.),np.log(MASS_MIN_INFALL*20.),-0.5])
labels=[r"$\ln m$",r"$\ln m_\mathrm{acc}$",r"$\ln R/R_{200}$"]
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

m,mInfall,R=flatchain.T
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
plt.loglog(x,y/4/np.pi/x**3/np.sum(R[m>MASS_MIN_INFALL]<1))
plt.subplot(122)
y=y/4/np.pi/x**3/np.sum(R[m>MASS_MIN_INFALL]<1)/(Host.density(x*Host.Rv)/Host.M*Host.Rv**3)
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

plt.figure()
plt.loglog(mInfall, cInfall, '.')
x=np.array(plt.xlim())
plt.plot(x, mean_concentration(x,'Duffy'))
plt.plot(x, mean_concentration(x,'Prada'))
plt.xlabel(r'$m_{\rm acc}$')
plt.ylabel(r'$c_{\rm acc}$')

plt.figure()
halos=[NFWHalo(mInfall[i], cInfall[i]) for i in range(nsub)]
rt=np.array([halos[i].radius(m[i]) for i in range(nsub)])
mpred=np.array([halos[i].mass(rt[i]) for i in range(nsub)])
L=np.array([halos[i].luminosity(rt[i]) for i in range(nsub)])
Lv=np.array([halos[i].luminosity(1000.) for i in range(nsub)])
f=mInfall>np.log(1e-8)
plt.hist(np.log10(L[f]),50,log=True)
#plt.figure()
#plt.hist(np.log10(L/Lv),50,log=True)
#x=L/Lv
#x=np.sort(x)[::-1]
#plt.loglog(x, 1.*np.arange(len(x))/len(x))

#theta=np.random.rand(nsub)*np.pi*2.
phi=np.arccos(np.random.rand(nsub)*2-1.) #random angle around the z-axis
Rp=R*np.sin(phi)

b0=1.6e-3*(Host.M*1e10/0.73)**0.39

plt.figure()
y,x=np.histogram(np.log(Rp), 20, weights=L)
y/=L[(Rp<1.)].sum() #normalize to Lsub(<Rv)
y/=(x[1]-x[0]) #dJ/dln(Rp)
x=np.exp(x[1:])
plt.loglog(x, y/x**2/2/np.pi, 'r', label='Sub') 
plt.plot(x, 16/np.pi/np.log(17.)/(1+16.*x**2), 'k', label='Gao')
plt.xlabel(r'$R_p/R_{200}$')
plt.ylabel(r'$j$')
y,x=np.histogram(np.log(Rp), 20, weights=Lv)
y/=Lv[(Rp<1.)].sum() #normalize to Lsub(<Rv)
y/=(x[1]-x[0]) #dJ/dln(Rp)
x=np.exp(x[1:])
plt.loglog(x, y/x**2/2/np.pi, 'g', label='SubAcc') 
y,x=np.histogram(np.log(Rp), 20)
y=y/(1.*(Rp<1.).sum()) #normalize to Nsub(<Rv)
y/=(x[1]-x[0]) #dJ/dln(Rp)
x=np.exp(x[1:])
plt.loglog(x, y/x**2/2/np.pi, 'b', label='SubNum') 
plt.legend(loc=3)

plt.figure()
f=(m>MASS_MIN_INFALL)&(R<1)
y,x=np.histogram(np.log(m)[f], 20, weights=L[f])
y/=(NUM_HOSTS*Host.L) #normalize
y/=(x[1]-x[0]) #db/dln(m)
x=np.exp(x[1:])
plt.loglog(x, y, 'r', label='Sub Cut') 
#plt.plot(x, (x)**(-0.226)*0.1, 'k--')
plt.plot(x, 0.226*b0*(x/1e-16)**-0.226, 'k-', label='Gao')

#f=(m>MASS_MIN_INFALL)&(Rp<1)
#yp,xp=np.histogram(np.log(m)[f], 20, weights=L[f])
#yp/=(NUM_HOSTS*Host.L) #normalize
#yp/=(xp[1]-xp[0]) #db/dln(m)
#xp=np.exp(xp[1:])
#plt.loglog(xp, yp, 'r--')

plt.xlabel(r'$m$')
plt.ylabel(r'$\mathrm{d}b/\mathrm{d}\ln m$')

f=(R<1)
y0,x0=np.histogram(np.log(mInfall)[f], 20, weights=L[f])
y0/=(NUM_HOSTS*Host.L) #normalize
y0/=(x0[1]-x0[0]) #db/dln(m)
x0=np.exp(x0[1:])
plt.loglog(x0, y0, 'r--', label='Infall Cut') 
plt.legend(loc=3)

#np.save('PradaMCDown1sig_min%.0e'%MASS_MIN_INFALL, np.array([x,y,x0,y0]))

a=np.load('min1e+00.npy')
b=np.load('min1e-05.npy')
c=np.load('min1e-10.npy')
d=np.load('min1e-16.npy')
data_old=np.hstack([a,b,c,d])
a=np.load('PradaMC_min1e+00.npy')
b=np.load('PradaMC_min1e-05.npy')
c=np.load('PradaMC_min1e-10.npy')
d=np.load('PradaMC_min1e-16.npy')
data=np.hstack([a,b,c,d])
a=np.load('PradaMCDown1sig_min1e+00.npy')
b=np.load('PradaMCDown1sig_min1e-05.npy')
c=np.load('PradaMCDown1sig_min1e-10.npy')
d=np.load('PradaMCDown1sig_min1e-16.npy')
dataDown=np.hstack([a,b,c,d])
plt.figure()
plt.loglog(data[2], data[3], 'gs')
plt.loglog(dataDown[2], dataDown[3], 'ro')
f=data[1]>0
par=np.polyfit(np.log(data[2][f]), np.log(data[3][f]), 1)
par[1]=np.exp(par[1])
x=np.array(plt.xlim())
plt.plot(x, par[1]*x**par[0], 'k')
plt.plot(x, 0.226*b0*(x/1e-16)**-0.226, 'k--')

plt.figure()
x=np.logspace(np.log10(MASS_MIN_INFALL), np.log10(MASS_MAX_INFALL), 10)
yp=[np.sum(L[(mInfall>xi)&(Rp<1.)]) for xi in x]
y=[np.sum(L[(mInfall>xi)&(R<1.)]) for xi in x]
ybound=[np.sum(L[(m>xi)&(R<1.)]) for xi in x]
y/=(NUM_HOSTS*Host.L)
yp/=(NUM_HOSTS*Host.L)
ybound/=(NUM_HOSTS*Host.L)
plt.loglog(x, y, 'g', label='Infall')
#plt.plot(x, yp, 'r--')
plt.loglog(x, ybound, 'r', label='Sub') 
#plt.plot(x, (x/x[0])**-0.226*yp[0], 'k-')
plt.xlabel(r'$m_{\rm min}$')
plt.ylabel(r'$b$')
plt.loglog(x, par[1]/(-par[0])*x**par[0], 'k-', label='Fit')
plt.plot(x, b0*(x/1e-16)**-0.226, 'k--', label='Gao')
plt.legend(loc=3)

###check truncation radius ##
#plt.figure()
#plt.hist(m/np.exp(flatchain[:,0]),100)