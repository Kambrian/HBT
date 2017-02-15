import numpy as np
import matplotlib.pyplot as plt
from nfw import *
plt.ion()
rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

H=NFWHalo(1e4)

r=np.logspace(-3,0.5)
y=H.density(r*H.Rv)*10

plt.figure()
plt.loglog(r, y*r**0.95, 'r--', label=r'$m$ selection')
Aref=np.logspace(-2.5, 0, 3)**0.95
plt.loglog(r, Aref[0]*y, 'k', lw=8, alpha=0.4)
plt.loglog(r, Aref[1]*y, 'k', lw=4, alpha=0.4, label=r'$m_{\rm acc}$ selection')
plt.loglog(r, Aref[2]*y, 'k', lw=2, alpha=0.4)
plt.xlabel(r'$R/R_{200}$')
plt.ylabel('Number Density')
plt.legend(loc=3)
plt.savefig(outdir+'/SubProfDemo.pdf')
  

