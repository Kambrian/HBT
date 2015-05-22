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

H=HaloData(sys.argv[1])

Rref=1
iInfall=0
fMin=1e-5
xbin=np.logspace(-2, np.log10(2), 20)

plt.figure()
mMin=fMin*H.Mvir

rHalo,denHalo,denHaloRef,denHaloErr=H.get_host_density(xbin, Rref)
plt.loglog(rHalo, denHalo, 'k-')

rSub,denSub,denSubRef0,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP), xbin, Rref)
f=(denSub/denHalo>0)&(rSub<1.8)&(rSub>0.6)
C=(denSub/denHalo)[f].mean()
denSub/=C
denSubRef0*=C
plt.plot(rSub[denSub>0],denSub[denSub>0],'b-',alpha=0.3, linewidth=4, label=r'$m_0/M_{200}>%s$'%fmtexp10(mMin/H.Mvir))

rSub,denSub,denSubRef1,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP)&(H.m>0), xbin, Rref)
C1=denSubRef0/denSubRef1
plt.plot(rSub[denSub>0],denSub[denSub>0]/C1,'r--', alpha=1, label='Resolved')

rSub,denSub,denSubRef2,denSubErr=H.get_sub_density((H.massTVV[:,iInfall]>mMin/H.mP)&(H.m==0), xbin, Rref)
C2=denSubRef0/denSubRef2
plt.plot(rSub[denSub>0],denSub[denSub>0]/C2,'g:', alpha=1, label='Orphan')

plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.legend(loc=1)
plt.xlabel(r'$R/R_{\rm 200}$')
plt.ylabel(r'$\rho(R)/\rho(R_{200})$')
#plt.savefig(outdir+'/SubprofInfall_Decompose.pdf')
