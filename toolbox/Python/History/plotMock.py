''' generate Monte-Carlo samples of subhaloes '''
import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from nfw import NFWHalo,mean_concentration
from myutils import myhist
plt.ion()
#rootdir='/gpfs/data/jvbq85/SubProf/'
rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'/plots/'

Galaxy=h5py.File(datadir+'MockGalaxy.hdf5','r')
Cluster=h5py.File(datadir+'MockCluster.hdf5','r')

def getLum(common,cL,mMin=1e-16, Lnorm=None):
  count,m,mAcc,R,Rp,phi,mu=common
  cAcc,rt,Lt,Lv=cL
  f=(m>mMin)
  xbin=np.linspace(np.log(1e-3), np.log(2),50)
  #stripped profile
  Lw=Lt*count
  y,x=myhist(np.log(R)[f], xbin, weights=Lw[f])
  if Lnorm is None:
	Lnorm=Lw[(R<1.)&f].sum() #normalize to Lsub(<Rv)
  y=y/Lnorm 
  #unstripped profile
  Lw=Lv*count
  yv,_=myhist(np.log(R)[f], x, weights=Lw[f])
  if Lnorm is None:
	Lnorm=Lw[(R<1.)&f].sum() #normalize to Lsub(<Rv)
  yv=yv/Lnorm #normalize to Lsub(<Rv) #dJ/dln(Rp)
  
  x=np.exp(x)
  area=np.diff(4./3.*np.pi*x**3)
  y=y/area
  yv=yv/area
  xmid=x[:-1]*np.exp((xbin[1]-xbin[0])/2)
  
  return xmid,y,yv

#2-d profile
def getSB(common,cL,mMin=1e-16, Lnorm=None):
  count,m,mAcc,R,Rp,phi,mu=common
  cAcc,rt,Lt,Lv=cL
  #z=R*np.cos(phi)
  f=(m>mMin)#&(z<1.)
  xbin=np.linspace(np.log(1e-3), np.log(2),50)
  #stripped profile
  Lw=Lt*count
  y,x=myhist(np.log(Rp)[f], xbin, weights=Lw[f])
  if Lnorm is None:
	Lnorm=Lw[(Rp<1.)&f].sum() #normalize to Lsub(<Rv)
  y=y/Lnorm 
  #unstripped profile
  Lw=Lv*count
  yv,_=myhist(np.log(Rp)[f], x, weights=Lw[f])
  if Lnorm is None:
	Lnorm=Lw[(Rp<1.)&f].sum() #normalize to Lsub(<Rv)
  yv=yv/Lnorm 
  
  x=np.exp(x)
  area=np.diff(np.pi*x**2)
  y=y/area
  yv=yv/area
  xmid=x[:-1]*np.exp((xbin[1]-xbin[0])/2)
  
  return xmid,y,yv

plt.figure()
plt.loglog()
halo=Galaxy
x,yl=getLum(halo['Common'], halo['Ludlow'], mMin=1e-6*halo['MHost'][...])[:2]
x,y,yv=getLum(halo['Common'], halo['Ludlow'])[:3]
plt.fill_between(x,yl,y, color='r', alpha=0.5, label='Galaxy')#, label=r'Ludlow14 $M(c)$')
#plt.plot(x,yv,'r--',label='Galaxy')
halo=Cluster
x,yl=getLum(halo['Common'], halo['Ludlow'], mMin=1e-6*halo['MHost'][...])[:2]
x,y,yv=getLum(halo['Common'], halo['Ludlow'])[:3]
plt.fill_between(x,yl,y,color='g', alpha=0.8, label='Cluster')
#plt.plot(x,yv,'g--',label='Cluster')
x,y=getLum(halo['Common'], halo['Maccio'])[:2]
plt.fill_between(x,yl,y,hatch='+',color="none", edgecolor='b', label=r'Maccio08 $M(c)$')#, label='Ludlow14 $M(c)$')
plt.plot(x, 4.53/(1+16.*x**2)**(3./2), 'k-', alpha=0.3, lw=3, label='Gao12') #this differs from the 2-D normalization by 11.3/4.53=2.5
a,b=0.95,-0.27
plt.plot(x, a*x**(a*x**b+b)*(1.+ b*np.log(x))/4/np.pi/x**3, 'k:', label='Pinzke11')

plt.xlabel(r'$x=R/R_{200}$')
plt.ylabel(r'$\mathrm{d}\tilde{L}_{\rm sub}/\mathrm{d}^3 x$')
plt.legend(loc=3)
plt.xlim([1e-3,2])
plt.ylim([1e-3,1e3])
plt.savefig(outdir+'LumDensity.pdf')

def BoostFactor(Mhalo, model='Gao', Mlim=1e-16):
  '''Mhalo and Mlim in 1e10 Msun/h'''
  if model=='Gao':
	return 1.6e-3*(Mhalo*1e10/0.73)**0.39*(Mlim/1e-16)**-0.226
  elif model=='Pinzke':
	return 0.76*0.023*(Mhalo/Mlim)**0.226
  else:
	raise "unknown model"

def SubLumProf(x, model='Gao'):
  '''x=R/R200; return dL/d^3x/L(x<1) of subhaloes'''
  if model=='Gao':
	return 4.53/(1+16.*x**2)**(3./2)
  elif model=='Pinzke':
	a,b=0.95,-0.27
	return a*x**(a*x**b+b)*(1.+ b*np.log(x))/4/np.pi/x**3
  else:
	raise "unknown model"

plt.figure()
plt.loglog()
#halo,name=Cluster,'Cluster' 
halo,name=Galaxy,'Galaxy'
host=NFWHalo(halo['MHost'][...])
x,y,yv=getLum(halo['Common'], halo['Ludlow'], Lnorm=host.Ls)[:3]
plt.plot(x,y, 'r-', label='Sub')
plt.plot(x,yv,'g--', label='Sub:NoStrip')
#x,y,yv=getLum(halo['Common'], halo['Ludlow'], mMin=1e-6*host.M, Lnorm=host.Ls)[:3]
#plt.plot(x,y,'g:')
#plt.plot(x,yv,'g--')
plt.plot(x, host.density(x*host.Rv)**2*host.Rv**3/host.Ls, 'b-', label='Halo')
plt.plot(x, BoostFactor(host.M,'Gao')*SubLumProf(x,'Gao'), 'k-', alpha=0.3, lw=3, label='Gao12') #this differs from the 2-D normalization by 11.3/4.53=2.5
plt.plot(x, BoostFactor(host.M,'Pinzke')*SubLumProf(x, 'Pinzke'), 'k:', label='Pinzke11')
#plt.plot(x, BoostFactor(host.M,'Pinzke', 1e-6*host.M)*SubLumProf(x, 'Pinzke'), 'k:')
plt.xlabel(r'$x=R/R_{200}$')
plt.ylabel(r'$\mathrm{d}(L/L_{\rm halo})/\mathrm{d}^3 x$')
plt.legend(loc=3)
plt.xlim([1e-3,2])
plt.savefig(outdir+'LumDensity'+name+'.pdf')

#plt.figure()
#halo=Galaxy
#x,y,yv,n=getSB(halo['Common'], halo['Ludlow'], mMin=1e-6*halo['MHost'][...])
#plt.loglog(x,y,'r', label='Galaxy')
#x,y,yv,n=getSB(halo['Common'], halo['Ludlow'])
#plt.loglog(x,y,'r--')#, label=r'Ludlow14 $M(c)$')
#halo=Cluster
#x,y,yv,n=getSB(halo['Common'], halo['Ludlow'], mMin=1e-6*halo['MHost'][...])
#plt.loglog(x,y,'g', label='Cluster')
#x,y,yv,n=getSB(halo['Common'], halo['Ludlow'])
#plt.loglog(x,y,'g--')#, label=r'Ludlow14 $M(c)$')
#x,y,yv,n=getSB(halo['Common'], halo['Maccio'])
#plt.loglog(x,y,'b--', label=r'Maccio08 $M(c)$')#, label=r'Ludlow14 $M(c)$')

#plt.plot(x, 16/np.pi/np.log(17.)/(1+16.*x**2), 'k-', alpha=0.7, lw=3, label='Gao12')
#plt.xlabel(r'$x_{\rm p}=R_{\rm p}/R_{200}$')
#plt.ylabel(r'$j$')
#plt.ylabel(r'$\mathrm{d}\tilde{L}_{\rm sub}/\mathrm{d}^2 x_{\rm p}$')
#plt.legend(loc=3)
#plt.xlim([1e-3,2])
#plt.savefig(outdir+'SB.pdf')

def getBoost(common, cL, MHost):
  count,m,mAcc,R,Rp,phi,mu=common
  cAcc,rt,Lt,Lv=cL
  Lw=count*Lt
  x=np.logspace(-16, np.log10(0.1*MHost), 50)
  y=[np.sum(Lw[(m>xi)&(R<1.)]) for xi in x]
  y/=NFWHalo(MHost).Ls*10. #stacked 10 haloes.
  return x,y
  
  
plt.figure()
halo=Galaxy
x,y=getBoost(halo['Common'], halo['Ludlow'], halo['MHost'][...])
h1,=plt.loglog(x*1e10, y, 'r', label=r'Ludlow14 $M(c)$')
x,y=getBoost(halo['Common'], halo['Maccio'], halo['MHost'][...])
h2,=plt.loglog(x*1e10, y, 'r--', label=r'Maccio08 $M(c)$')
#plt.plot(1e-6, BoostFactor(halo['MHost'][...],'Gao'), 'k*', mfc="none", markersize=12)
plt.errorbar(3e-6, BoostFactor(halo['MHost'][...],'Gao'), xerr=2e-6, xuplims=True, capsize=3, color='r')
plt.plot(x*1e10, BoostFactor(halo['MHost'][...],'Pinzke', x), 'r:')
halo=Cluster
x,y=getBoost(halo['Common'], halo['Ludlow'], halo['MHost'][...])
plt.loglog(x*1e10, y, 'g')
x,y=getBoost(halo['Common'], halo['Maccio'], halo['MHost'][...])
plt.loglog(x*1e10, y, 'g--')
plt.errorbar(3e-6, BoostFactor(halo['MHost'][...],'Gao'), xerr=2e-6, xuplims=True, capsize=3, color='g')
plt.plot(x*1e10, BoostFactor(halo['MHost'][...],'Pinzke', x), 'g:')
plt.xlabel(r'$m_{\rm lim}[M_\odot/h]$')
plt.ylabel('Boost Factor')
x=np.logspace(-4,3)
plt.plot(x*1e10, 8000*(x/1e-16)**-0.226, color=[0.7,0.7,0.7], label=r'slope $-0.226$')
plt.legend(loc=3)
plt.savefig(outdir+'/Boost.eps')

#TODO: solve the inconsistency?????? state that the model likelihood overestimate the boost factor by a factor of ...