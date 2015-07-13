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
from myutils import myhist

#plt.ion()

rootdir='/gpfs/data/jvbq85/SubProf/'
#rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

MASS_HOST=6.7e4
#CONC_HOST=5
NUM_SUBS=float(sys.argv[1]) #1e5
C_OFFSET=0.

MASS_MIN_INFALL=float(sys.argv[2]) #in units of 1e10Msun/h
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

lnmmin,lnmmax=np.log(MASS_MIN_INFALL), np.log(MASS_MAX_INFALL)
lnmAcc=np.random.rand(NUM_SUBS)*(lnmmax-lnmmin)+lnmmin
mAcc=np.exp(lnmAcc)

Host.Msample=Host.mass(R_MAX*Host.Rv)-Host.mass(R_MIN*Host.Rv)
print 'Fraction of mass sampled inside Rvir:', (Host.mass(Host.Rv)-Host.mass(R_MIN*Host.Rv))/Host.M
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