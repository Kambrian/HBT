''' generate Monte-Carlo samples of subhaloes '''
import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from nfw import NFWHalo,mean_concentration
from MbdIO import cross_section
import emcee
from myutils import myhist

plt.ion()

#rootdir='/gpfs/data/jvbq85/SubProf/'
rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

#mlist=np.logspace(8,15, 20)/1e10
ithread=int(sys.argv[1])
MASS_HOST=1e2 #mlist[ithread]
#CONC_HOST=5
NUM_SUBS=1e5 #float(sys.argv[1]) #1e5

MASS_MIN_INFALL=1e-16 #1e-7*MASS_HOST #float(sys.argv[2]) #in units of 1e10Msun/h
MASS_MAX_INFALL=MASS_HOST/10. #min(MASS_MIN_INFALL*1e4, MASS_HOST/10.)
R_MAX=2.
R_MIN=0.

sigma=1.1
fs=0.55
A=0.1*MASS_HOST**-0.02
alpha=0.95
mustar=0.5*MASS_HOST**-0.03
beta=1.7*MASS_HOST**-0.04

Host=NFWHalo(MASS_HOST)#,Chost)
Host.L=Host.luminosity(1000.)

#=================================Generate mu and R=======================================
def lnPDF(x):
  ''' R in units of Rv'''
  lnmu,lnR=x
  lnmubar=np.log(mustar)+beta*lnR
  dlnmu=lnmu-lnmubar
  if lnmu>0:
	return -np.inf
  if dlnmu>np.log(4.2):
	return -np.inf
  if lnR>np.log(R_MAX): #log(3)=1.1
	return -np.inf
  
  lnPDFmu=-0.5*(dlnmu/sigma)**2
  lnPDFR=3.*lnR+np.log(Host.density(np.exp(lnR)*Host.Rv))
  return lnPDFmu+lnPDFR

nwalkers=8
nburn=200
nsteps=int(NUM_SUBS/nwalkers+nburn)
print 'running %d steps'%nsteps
ndim=2

x00=np.array([-0.5,-0.5])
labels=[r"$\ln \mu$",r"$\ln R/R_{200}$"]
x0=np.kron(np.ones([nwalkers,1]),x00)#repmat to nwalkers rows
x0+=(np.random.rand(ndim*nwalkers).reshape(nwalkers,ndim)-0.5)*0.1 #random offset, [-0.5,0.5]*0.1
sampler=emcee.EnsembleSampler(nwalkers,ndim,lnPDF)
sampler.run_mcmc(x0,nsteps)
from matplotlib.pyplot import *
#ion()
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

mu,R=flatchain.T
plt.figure()
plt.subplot(121)
y,x=myhist(np.log(R),30)
y/=(x[1]-x[0])
x=np.exp(x[1:])
plt.loglog(x,y/4/np.pi/x**3/np.sum(R<1))
plt.plot(x, Host.density(x*Host.Rv)/Host.M*Host.Rv**3)
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$\mathrm{d}P/\mathrm{d}^3(R/R_{200})$')
plt.subplot(122)
plt.loglog(R, mu, '.')
xbin=np.logspace(-3, np.log(R_MAX),20)
p,xmid=cross_section(R,mu,xbin,[(100-68.3)/2,50,(100+68.3)/2])
plt.plot(xmid, p[1], 'r-', lw=2, label='Resolved')
plt.plot(xmid, p[[0,2]].T, 'r--',lw=2)
plt.plot(xmid, mustar*xmid**beta, 'k-')
plt.plot(xmid, mustar*xmid**beta*np.exp(sigma), 'k--')
plt.plot(xmid, mustar*xmid**beta/np.exp(sigma), 'k--')
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$\mu$')
#==========projections==========================
phi=np.arccos(np.random.rand(nsub)*2-1.) #random angle around the z-axis
Rp=R*np.sin(phi)
#======================================================Now generate Infall mass==================================
lnmmin,lnmmax=np.log(MASS_MIN_INFALL), np.log(MASS_MAX_INFALL)
lnmAcc=np.random.rand(NUM_SUBS)*(lnmmax-lnmmin)+lnmmin
mAcc=np.exp(lnmAcc)

Host.Msample=Host.mass(R_MAX*Host.Rv)
NUM_SUBS_PRED=fs*A*Host.Msample*(MASS_MIN_INFALL**-alpha-MASS_MAX_INFALL**-alpha)/alpha #infall number scale with host mass (in the given radial range: Mv*M/Mv=M)

#count=mAcc**-alpha
#count=count/count.sum()*NUM_SUBS_PRED #w/sum(w)*Npred, equals to dN/dlnm*[delta(lnm)/N] as N->inf
count=fs*A*Host.Msample*mAcc**-alpha*(lnmmax-lnmmin)/NUM_SUBS #the weight is dN/dlnm*[delta(lnm)/N]
print np.sum(count), NUM_SUBS_PRED

y,x=myhist(np.log(mAcc),50,weights=count)
f=y>10
print 'mass func slope:', np.polyfit(x[f],np.log(y)[f],1)[0]
plt.figure()
y,x=myhist(np.log(mAcc),50,weights=count)
xmax=np.exp(x[-1])
x=np.exp(x[:-1])
y=y[::-1].cumsum()[::-1]
plt.loglog(x,y/fs*x, 'ro', label='infall')
plt.plot(x, (x**-alpha-xmax**-alpha)/alpha*A*Host.Msample*x, 'k')
plt.xlabel(r'$m$')
plt.ylabel(r'$mN(>m)$')	
#========== generate final mass =========================
m=mAcc*mu

plt.figure()
plt.subplot(121)
y,x=myhist(np.log(R),30)
y/=(x[1]-x[0])
x=np.exp(x[1:])
plt.loglog(x,y/4/np.pi/x**3/np.sum(R<1))
plt.plot(x, Host.density(x*Host.Rv)/Host.M*Host.Rv**3)
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$\mathrm{d}P/\mathrm{d}^3(R/R_{200})$')
y,x=myhist(np.log(R[m>MASS_MIN_INFALL]),30, weights=count[m>MASS_MIN_INFALL])
y/=(x[1]-x[0])
x=np.exp(x[1:])
y=y/4/np.pi/x**3/np.sum(count[(m>MASS_MIN_INFALL)&(R<1)])
plt.loglog(x,y)
plt.subplot(122)
y=y/(Host.density(x*Host.Rv)/Host.M*Host.Rv**3)
plt.loglog(x, y, 'o')
pars=np.polyfit(np.log(x[y>0]), np.log(y[y>0]), 1)
plt.plot(x, np.exp(pars[1])*x**(alpha*beta), '-')
plt.legend(loc=2)
#============================gen concentration==========================
sigmaC=0.3#log10-scatter=0.13, ln-scatter=0.3 
deltalnC=np.random.normal(0, sigmaC, NUM_SUBS)
#=========================combine data===========
common=np.array([count,m,mAcc,R,Rp,phi,mu])

lncMaccio=np.log(mean_concentration(mAcc,'MaccioW1'))
#lncDuffy=np.log(mean_concentration(mAcc,'Duffy'))
#lncPrada=np.log(mean_concentration(mAcc,'Prada'))
lncLudlow=np.log(mean_concentration(mAcc,'Ludlow'))
def eval_luminosity(lncbar):
  cAcc=np.exp(lncbar+deltalnC)
  #============================evaluate truncation and luminosity=========
  halos=[NFWHalo(mAcc[i], cAcc[i]) for i in range(nsub)]
  rt=np.array([halos[i].radius(m[i]) for i in range(nsub)])
  mpred=np.array([halos[i].mass(rt[i]) for i in range(nsub)])
  Lt=np.array([halos[i].luminosity(rt[i]) for i in range(nsub)])
  Lv=np.array([halos[i].luminosity(1000.) for i in range(nsub)])
  return cAcc,rt,Lt,Lv
Maccio=eval_luminosity(lncMaccio)
#Duffy=eval_luminosity(lncDuffy)
#Prada=eval_luminosity(lncPrada)
#Prada1=eval_luminosity(lncPrada-sigmaC)
Ludlow=eval_luminosity(lncLudlow)
Ludlow1=eval_luminosity(lncLudlow-sigmaC*2)
#Ludlow1=eval_luminosity(lncLudlow-0.5)

b0=1.6e-3*(Host.M*1e10/0.73)**0.39

def plot_lum(cL,mMin=MASS_MIN_INFALL):
  z=R*np.cos(phi)
  L,Lv=cL[2],cL[3]
  plt.figure()
  Lw=L*count
  f=mAcc>mMin
  y,x=myhist(np.log(Rp)[f], 20, weights=Lw[f])
  y/=Lw[(Rp<1.)&f].sum() #normalize to Lsub(<Rv)
  y/=(x[1]-x[0]) #dJ/dln(Rp)
  x=np.exp(x[1:])
  plt.loglog(x, y/x**2/2/np.pi, 'r:', label='Sub-NoSelect') 
  f=(m>mMin)&(z<1.)
  y,x=myhist(np.log(Rp)[f], 20, weights=Lw[f])
  y/=Lw[(Rp<1.)&f].sum() #normalize to Lsub(<Rv)
  y/=(x[1]-x[0]) #dJ/dln(Rp)
  x=np.exp(x[1:])
  plt.loglog(x, y/x**2/2/np.pi, 'r', label='Sub') 
  plt.plot(x, 16/np.pi/np.log(17.)/(1+16.*x**2), 'k', label='Gao')
  plt.xlabel(r'$R_p/R_{200}$')
  plt.ylabel(r'$j$')
  Lvw=Lv*count
  y,x=myhist(np.log(Rp)[f], 20, weights=Lvw[f])
  y/=Lvw[(Rp<1.)&f].sum() #normalize to Lsub(<Rv)
  y/=(x[1]-x[0]) #dJ/dln(Rp)
  x=np.exp(x[1:])
  plt.loglog(x, y/x**2/2/np.pi, 'g', label='SubAcc') 
  y,x=myhist(np.log(Rp)[f], 20)
  y=y/(1.*f[(Rp<1.)].sum()) #normalize to Nsub(<Rv)
  y/=(x[1]-x[0]) #dJ/dln(Rp)
  x=np.exp(x[1:])
  plt.loglog(x, y/x**2/2/np.pi, 'b', label='SubNum') 
  plt.legend(loc=3)

  plt.figure()
  x=np.logspace(np.log10(MASS_MIN_INFALL), np.log10(MASS_MAX_INFALL), 50)
  yp=[np.sum(Lw[(mAcc>xi)&(Rp<1.)]) for xi in x]
  y=[np.sum(Lw[(mAcc>xi)&(R<1.)]) for xi in x]
  ybound=[np.sum(Lw[(m>xi)&(R<1.)]) for xi in x]
  y/=Host.L
  yp/=Host.L
  ybound/=Host.L
  plt.loglog(x, y, 'g', label='Acc')
  #plt.plot(x, yp, 'r--')
  plt.loglog(x, ybound, 'r', label='Sub') 
  #plt.plot(x, (x/x[0])**-0.226*yp[0], 'k-')
  plt.xlabel(r'$m_{\rm min}$')
  plt.ylabel(r'$b$')
  #plt.loglog(x, par[1]/(-par[0])*x**par[0], 'k-', label='Fit')
  #plt.plot(x, b0*(x/1e-16)**-0.226, 'k--', label='Gao')
  plt.plot(x, 0.76*0.023*(MASS_HOST/x)**0.226, 'k--', label='Pinzke')
  plt.legend(loc=3)
  plt.ylim(1e-4,1e3)

plot_lum(Maccio)
plot_lum(Ludlow)
plot_lum(Ludlow1)
#============================save========================================
outfile=h5py.File(datadir+'MockHalo-%d.hdf5'%ithread,'w')
outfile.create_dataset('MHost',data=MASS_HOST)
dset=outfile.create_dataset('Common', data=common)
dset.attrs['rows']='count,m,mAcc,R,Rp,phi,mu'

dset=outfile.create_dataset('Maccio',data=Maccio)
dset.attrs['rows']='cAcc,rt,Lt,Lv'
dset.attrs['sigma_lnC']=sigmaC

dset=outfile.create_dataset('Ludlow',data=Ludlow)
dset.attrs['rows']='cAcc,rt,Lt,Lv'
dset.attrs['sigma_lnC']=sigmaC

dset=outfile.create_dataset('Ludlow1',data=Ludlow1)
dset.attrs['rows']='cAcc,rt,Lt,Lv'
dset.attrs['sigma_lnC']=sigmaC
outfile.close()
