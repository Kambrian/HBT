import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'
  
A1=HaloData('A1')
A2=HaloData('A2')
A3=HaloData('A3')
A4=HaloData('A4')
A5=HaloData('A5')
A2H=HaloData('A2', None)
A4H=HaloData('A4', None)

nbin=20
xbin=np.logspace(np.log10(0.1), np.log10(500*0.73), nbin)/A2.Rvir
rHalo,denHalo,denHaloRef,denHaloErr=A2.get_host_density(xbin, 1)
rSub,denSub,denSubRef,denSubErr=A1.get_sub_density(A1.m>1000, xbin, 1)

AqRv={'A':179.49, 'B':137.02, 'C':177.26, 'D':177.28, 'E':154.96, 'F':152.72}
PhRv={'A':1414.,'B':1526., 'C':1332., 'D':1386., 'E':1369., 'F':1509., 'G':1704.}

AqRat={}
PhRat={}
for H in 'ABCDE':
  AqRat[H]=np.loadtxt(datadir+'Aquarius/Aq'+H+'2rat1000.dat').T
  AqRat[H][0]/=AqRv[H]
  sp=interp1d(AqRat[H][0][AqRat[H][1]>0], AqRat[H][1][AqRat[H][1]>0])
  AqRat[H][1]/=sp(1)
  
for H in 'ABCDEFG':
  PhRat[H]=np.loadtxt(datadir+'Phoenix/Ph'+H+'2rat1000.dat').T
  PhRat[H][0]/=PhRv[H]
  sp=interp1d(PhRat[H][0][PhRat[H][1]>0], PhRat[H][1][PhRat[H][1]>0])
  PhRat[H][1]/=sp(1)

plt.figure()
plt.loglog(rSub, denSub/denHalo, 'o')
for H in 'ABCDE':
  plt.loglog(AqRat[H][0], AqRat[H][1], '-')

for H in 'ABCDEFG':
  plt.loglog(PhRat[H][0], PhRat[H][1], ':')

#CoCo
#for color,m in [('r','11'),('g','12'),('b','13')]:
  #halo=np.loadtxt(datadir+'CoCo/halo'+m+'.stats')
  #sub=np.loadtxt(datadir+'CoCo/sat'+m+'.stats')
  #r=sub[:,1]
  #f=(r<1)
  #f2=(r>0.3)
  #x=r[f]
  #ymean=np.exp(np.interp(np.log(r), np.log(halo[:,1]), np.log(halo[:,3])))
  #ymed=np.exp(np.interp(np.log(r), np.log(halo[:,1]), np.log(halo[:,8])))
  ##plt.errorbar(x, sub[f,3]/ymean, sub[f, 14]/ymean, fmt='r')
  #sig=sub[:, 16]
  ##sig=sub[:, 18]-sub[:,16]
  ##sig=sub[:, 8]/np.sqrt(sub[:,19])
  #yref=np.exp(np.interp(np.log(1.), np.log(x), np.log(sub[f,8]/ymed[f])))
  ##plt.errorbar(x, sub[f,8]/ymed[f]/yref, (sig/ymed)[f]/yref, color=color, fmt='o')
  #plt.fill_between(x[f], sub[f,11]/ymed[f]/yref, sub[f,12]/ymed[f]/yref, color=color, alpha=0.2, label=m)
  #polypar,polyV=np.polyfit(np.log(r[f&f2]), np.log(sub[:,8]/ymed/yref)[f&f2],1, w=1/(sig/ymed/yref)[f&f2], cov=True)
  ##plt.plot(x, np.exp(np.polyval(polypar, np.log(x))), '--', color=color, label=m+r':$\gamma=%.2f\pm %.2f$'%(polypar[0],np.sqrt(polyV[0,0])))
  #plt.xscale('log')
  #plt.yscale('log', nonposy='mask')
  
plt.xlabel(r'$R/R_{200}$')
plt.ylabel(r'$\rho_{\rm sub}/\rho_{\rm host}$')
plt.xlim([1e-2,3])
plt.ylim([1e-2,3])

f=(denSub/denHalo>0)&(rSub<1)
w=None
polypar,polyV=np.polyfit(np.log(rSub[f]), np.log(denSub/denHalo)[f],1, w=w, cov=True)
plt.plot(rHalo[denHalo>0], ((rHalo)**polypar[0]*np.exp(polypar[1]))[denHalo>0],'r--', label='Model')
plt.legend()