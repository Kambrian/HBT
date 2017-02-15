''' generate Monte-Carlo samples of subhaloes '''
import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from nfw import NFWHalo,mean_concentration
from matplotlib.patches import Ellipse
from myutils import myhist,contour_handle,weighted_crosssection
plt.ion()
#rootdir='/gpfs/data/jvbq85/SubProf/'
rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'/plots/'

Galaxy=h5py.File(datadir+'MockGalaxy.hdf5','r')
Cluster=h5py.File(datadir+'MockCluster.hdf5','r')
fs=0.55
#===========================gen stellar mass====================
MvInMvb=0.893
MvcInMvb=0.733
Pwang06=[2*10.**10.27/3.33e11*0.73**2*MvcInMvb,3.33e11/MvcInMvb,1-2.59,1-0.276,1]#2006 original paper
Pwang2=[2*10**10.48/6.31e11*0.73**2*MvcInMvb,6.31e11/MvcInMvb,1-2.42,1-0.29,1]#vvds
Pwang3=[2*10**10.17/3.21e11*0.73**2*MvcInMvb,3.21e11/MvcInMvb,1-2.42,1-0.29,1]#unified, close to ling
Pwang4=[2*10**10.23/3.43e11*0.73**2*MvcInMvb,3.43e11/MvcInMvb,1-2.56,1-0.34,1]#DR7, 2013 central
PwangSat=[2*10**10.30/5.23e11*0.73**2*MvcInMvb,5.23e11/MvcInMvb,1-1.99,1-0.298,1] #DR7, 2013 sat
Pguo=[0.129*0.73*MvcInMvb,10**11.4*0.73/MvcInMvb,-0.926,0.261,2.440]#close to Yang at >M0 ; M200c
Pling=[2*10**-1.73,10**11.70,-1.16,0.71,1]#M200b
halo2starP=lambda M,A,M0,alpha,beta,gamma: A/((M/M0)**alpha+(M/M0)**beta)**gamma*M/0.73 #Msun/h
halo2star=lambda M: halo2starP(M*1e10, *PwangSat) #input M in units of 1e10Msun/h, output in Msun/h
sigmaMstar=0.192/np.log10(np.exp(1))

def lens_Mstar(H, MstarMin=1e10):
  count,m,mAcc,Rp=H['Common'][[0,1,2,4]]
  nsubs=len(m)
  lnMstar=np.log(halo2star(mAcc))
  deltalnMstar=np.random.normal(0, sigmaMstar, nsubs)
  Mstar=np.exp(lnMstar+deltalnMstar)
  f=Mstar>MstarMin
  xbin=np.logspace(-2, np.log10(2), 20)
  x=xbin[:-1]*xbin[1]/xbin[0]
  n,_=myhist(Rp[f], xbin, weights=count[f])
  med=myhist(Rp[f], xbin, weights=(np.log(m*1e10)*count)[f])[0]/n
  med=np.exp(med)
  mean=myhist(Rp[f], xbin, weights=(m*1e10*count)[f])[0]/n
  p=weighted_crosssection(Rp[f], xbin, m[f]*1e10, [15, 50, 85], count[f])
  print np.log10(Mstar[f]).mean()
  #plt.fill_between(x, med, mean, color=color, alpha=alpha, label=label, **kwarg)
  return x*NFWHalo(H['MHost'][...]).Rv/1e3, med, mean,p.T

def lens_MstarRat(H,MstarMin=1e10):
  count,m,mAcc,Rp=H['Common'][[0,1,2,4]]
  nsubs=len(m)
  lnMstar=np.log(halo2star(mAcc))
  deltalnMstar=np.random.normal(0, sigmaMstar, nsubs)
  Mstar=np.exp(lnMstar+deltalnMstar)
  f=Mstar>MstarMin
  xbin=np.logspace(-2, np.log10(2), 20)
  x=xbin[:-1]*xbin[1]/xbin[0]
  n,_=myhist(Rp[f], xbin, weights=count[f])
  ntmp,_=myhist(Rp[f]*NFWHalo(H['MHost'][...]).Rv/1e3, [0.1,0.3,0.6,0.9], weights=count[f])
  print ntmp
  mtmp,_=myhist(Rp[f]*NFWHalo(H['MHost'][...]).Rv/1e3, [0.1,0.3,0.6,0.9], weights=(np.log10(Mstar)*count)[f])
  print mtmp/ntmp
  med=myhist(Rp[f], xbin, weights=(np.log(m*1e10/Mstar)*count)[f])[0]/n
  med=np.exp(med)
  mean=myhist(Rp[f], xbin, weights=(m*1e10/Mstar*count)[f])[0]/n
  p=weighted_crosssection(Rp[f], xbin,(m*1e10/Mstar)[f], [15, 50, 85], count[f])
  #plt.fill_between(x, med, mean, color=color, alpha=alpha, label=label, **kwarg)
  return x*NFWHalo(H['MHost'][...]).Rv/1e3, med, mean,p.T

def lens(H):
  count,m,mAcc,Rp=H['Common'][[0,1,2,4]]
  mu=m/mAcc
  xbin=np.logspace(-2.3, np.log10(2), 20)
  x=xbin[:-1]*xbin[1]/xbin[0]
  n,_=myhist(Rp, xbin, weights=count)
  med=myhist(Rp, xbin, weights=np.log(mu)*count)[0]/n
  med=np.exp(med)
  mean=myhist(Rp, xbin, weights=mu*count)[0]/n
  #plt.fill_between(x, med, mean, color=color, alpha=alpha, label=label, **kwarg)
  return x, med, mean

plt.figure()
x,y1,y2=lens(Galaxy)
plt.fill_between(x,y1,y2,color='r', alpha=0.3)
h1=Ellipse((0,0),0,0,fill=True, color='r', label='galaxy')
x,y1,y2=lens(Cluster)
plt.fill_between(x,y1,y2, color='none', edgecolor='g', hatch='\\')
h2=Ellipse((0,0),0,0,fill=False, color='g', hatch='\\', label='cluster')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$R_{\rm p}/R_{200}$')
plt.ylabel(r'$m/m_{\rm acc}$')
plt.axis([1e-2,2,1e-3,1])
plt.legend((h1,h2), ('galaxy','cluster'),loc=2)
#plt.savefig(outdir+'/SubLens.pdf')

Galaxy1=h5py.File(datadir+'MockGalaxy-1.hdf5','r') #1e12
Cluster1=h5py.File(datadir+'MockCluster-1.hdf5','r') #1e14
Cluster2=h5py.File(datadir+'MockCluster-2.hdf5','r') #1e15
Cluster6=h5py.File(datadir+'MockCluster-6.hdf5','r') #1e16, Mmax=Mhalo/2
Cluster4=h5py.File(datadir+'MockCluster-4.hdf5','r') #1e13, Mmax=Mhalo/2
plt.figure()
x,y1,y2,p=lens_Mstar(Cluster6)
plt.plot(x,y1, 'r-')
plt.plot(x, y2, 'r--')
plt.fill_between(x,y1,y2, color='none', edgecolor='r', lw=0.0, hatch='\\')
x,y1,y2,p=lens_Mstar(Cluster1)
plt.plot(x,y1, 'g-', label='Median')
plt.plot(x, y2, 'g--', label='Mean')
#plt.fill_between(x,y1,y2, color='none', edgecolor='g', hatch='\\')
plt.fill_between(x, p[0], p[2], color='g', alpha=0.3)
y=np.array([11.28, 11.98, 12.51])
yerr=np.array([[0.44, 0.21, 0.14],[-0.46,-0.19,-0.14]])
ylim=np.array([y+yerr[0], y+yerr[1]])
yerr=np.array([10**y-10**ylim[1], 10**ylim[0]-10**y])
xerr=[0.1,0.15,0.15]
plt.errorbar([0.2, 0.45, 0.75], 10**y/fs,ms=12, fmt='kv', alpha=0.3, label='Li15 shifted')
plt.errorbar([0.2, 0.45, 0.75], 10**y, xerr=xerr, yerr=yerr, fmt='ko', ms=8, label='Li15 original')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$R_{\rm p}[Mpc/h]$')
plt.ylabel(r'$m[M_\odot/h]$')
plt.axis([0.08,1,1e10,1e13])
plt.legend(loc=4)
#plt.savefig(outdir+'/SubLensMstar.pdf')

plt.figure()
x,y1,y2,p=lens_MstarRat(Cluster6)
plt.plot(x,y1, 'r-',alpha=0.3)
plt.plot(x, y2, 'r--',alpha=0.3)
plt.fill_between(x,y1,y2, color='none', edgecolor='r', linewidth=0.0, hatch='\\',alpha=0.3)
x,y1,y2,p=lens_MstarRat(Cluster1)
plt.plot(x,y1, 'g-', label='Median')
plt.plot(x, y2, 'g--', label='Mean')
#plt.fill_between(x,y1,y2, color='none', edgecolor='g', hatch='\\')
plt.fill_between(x, p[0], p[2], color='g', alpha=0.3)
y=np.array([3.48, 14.97, 41.15])
yerr=np.array([[2.48,6.37,12.51], [4.48, 6.45, 12.55]])
xerr=[0.1,0.15,0.15]
fs=0.55
plt.errorbar([0.2, 0.45, 0.75], y/fs,ms=12, fmt='kv', alpha=0.3, label='Li15 shifted')
plt.errorbar([0.2, 0.45, 0.75], y, xerr=xerr, yerr=yerr, fmt='ko', ms=8, label='Li15 original')
#plt.plot([0.2,0.45,0.75], 10**np.array([11.28-10.82,11.98-10.86,12.51-10.92]), 'ks-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$R_{\rm p}[{\rm Mpc}/h]$')
plt.ylabel(r'$m/m_{\star}$')
plt.axis([0.08,1,1e-1,1e2])
plt.legend(loc=4)
#plt.savefig(outdir+'/SubLensMstarRat.pdf')