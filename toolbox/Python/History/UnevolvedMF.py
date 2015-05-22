import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d, UnivariateSpline
from MbdIO import *
plt.ion()

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'
  
##plot profiles
Rref=1
nbin=50
iInfall=0

plt.figure()
H=HaloData('AqA1')
xmin=100*H.mP/H.Mvir
m=H.massTVV[:,iInfall]
data=m[(H.r<Rref*H.Rvir)&(m>0)]*H.mP/H.Mvir
#xbin=np.logspace(np.log10(data.min()), np.log10(data.max()), nbin)
xbin=np.logspace(-7,-2,nbin)
n,_=np.histogram(data, xbin, weights=None)
y=n/np.log(xbin[1]/xbin[0])/H.Mvir
x=xbin[:-1]*np.sqrt(xbin[1]/xbin[0])*H.Mvir

f=(xbin[:-1]>xmin)
#plt.step(xbin[:-1][f], y[f], where='post', label=H.name)
plt.plot(x[f],(y*x)[f],'r')
plt.loglog()
f=f&(n>4)
#w=1./np.sqrt(n)[f]
w=None
pars=powerlaw_fit(x[f], y[f], w=w)
print 'alpha=', pars[0][1], '+-', pars[1][1]
print 'A=', pars[0][2], '+-', pars[1][2]
plt.plot(x, pars[0][2]*x**pars[0][1]*x, 'k-',label=r'$\alpha=%.2f$'%(-pars[0][1]))
		 
#colors={A1:'k',A2:'r',A3:'g',A4:'b',A5:'c'}
#for H in [A3,A5]:
  #m=H.massTVV[:,iInfall]
  #data=m[(H.r<Rref*H.Rvir)&(m>100)]*H.mP/H.Mvir
  #n,xbin=np.histogram(data, np.logspace(np.log10(data.min()), np.log10(data.max()), nbin/3), weights=None)
  #y=n/np.log(xbin[1]/xbin[0])
  #x=xbin[:-1]*(xbin[1]/xbin[0])
  #plt.step(xbin[:-1], y, where='post', color=colors[H], label=H.name)
  #f=y>0
  ##f=f&(x<1e-4)
  #w=np.sqrt(n)[f]
  ##w=None
  #polypar,polyV=np.polyfit(np.log(x)[f], np.log(y)[f],1, w=w , cov=True)
  #plt.plot(x, np.exp(np.polyval(polypar, np.log(x))), '-', color=colors[H], label=r'$\alpha=%.2f$'%(-polypar[0]))

plt.legend(loc=1)
plt.xlabel(r'$m_{\rm acc}/M_{200}$')
plt.ylabel(r'$(m_{\rm acc}/M_{200})\mathrm{d}N/\mathrm{d}\ln m_{\rm acc}$')
#plt.savefig(outdir+'InfallMF.eps')


plt.figure()
H=HaloData('AqA1')
xmin=100*H.mP/H.Mvir
m=H.massTVV[:,iInfall]
data=m[(H.r<Rref*H.Rvir)&(m>0)]*H.mP/H.Mvir
xbin=np.logspace(-7,0,nbin)
n,_=np.histogram(data, xbin, weights=None)
y=n[::-1].cumsum()[::-1]
x=xbin[:-1]*H.Mvir

f=(xbin[:-1]>xmin)
#plt.step(xbin[:-1][f], y[f], where='post', label=H.name)
plt.plot(x[f],(y*x)[f],'r')
plt.loglog()
f=f&(n>4)
#w=np.sqrt(n)[f]
w=None
pars=powerlaw_fit(x[f], y[f], w=w)
print pars
plt.plot(x, pars[0][2]*x**pars[0][1]*x, 'k-',label=r'$\alpha=%.2f$'%(-pars[0][1]))
alpha=-pars[0][1]
A=pars[0][2]/H.Mvir*alpha
print 'alpha=', alpha, 'A=', A

plt.figure()
xbin=np.logspace(-10,0,nbin)
Mv=[134.3, 134.5, 134.1, 134.2, 135.7]
Rv=[179.41,179.49,179.31, 179.36, 180.05]
for i in range(5):
  H=HaloData('AqA%d'%(i+1))
  H.Mvir=Mv[i]
  H.Rvir=Rv[i]
  m=H.massTVV[:,0]
  #m=H.flag
  data=m[(H.r<Rref*H.Rvir)&(m>0)]*H.mP/H.Mvir
  n,_=np.histogram(data, xbin, weights=data)
  y=n[::-1].cumsum()[::-1]
  x=xbin[:-1]
  f=x>20*H.mP/H.Mvir
  plt.plot(x[f],y[f], label=H.name)
plt.legend(loc=3)  
plt.xscale('log')
plt.xlabel(r'$m_{\rm acc}/M_{200}$')
plt.ylabel(r'Mass Fraction $(>m_{\rm acc})$')