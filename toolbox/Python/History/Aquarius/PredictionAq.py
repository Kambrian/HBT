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

AqRv={'A':179.49, 'B':137.02, 'C':177.26, 'D':177.28, 'E':154.96, 'F':152.72}
#AqMv={'A':}
#AqmP={'A':}

PhRv={'A':1414.,'B':1526., 'C':1332., 'D':1386., 'E':1369., 'F':1509., 'G':1704.}
PhMv={'A':6.570e4,'B':8.255e4,'C':5.495e4, 'D':6.191e4, 'E':5.969e4, 'F':7.997e4, 'G':1.150e5}
PhmP={'A':5.084e-4, 'B':6.127e-4, 'C':4.605e-4, 'D':4.721e-4, 'E':4.425e-4, 'F':6.219e-4, 'G':8.599e-4}

PhRat={}  
PhHalo={}
for H in 'ABCDEFG':
  PhHalo[H]=np.loadtxt(datadir+'Phoenix/Ph'+H+'2Halo.dat').T
  PhRat[H]=np.loadtxt(datadir+'Phoenix/Ph'+H+'2rat1000.dat').T
  PhRat[H][0]/=PhRv[H]

AqRat={}  
AqHalo={}
for H in 'ABCDE':
  AqHalo[H]=np.loadtxt(datadir+'Aquarius/Aq'+H+'2Halo.dat').T
  AqRat[H]=np.loadtxt(datadir+'Aquarius/Aq'+H+'2rat1000.dat').T
  AqRat[H][0]/=AqRv[H]
  
Rref=1
nbin=30
Rconv=0.3
nMinInfall=0
iInfall=0

A=0.085
alpha=0.95
mustar=0.38
beta=1.2
sigma=1.1
fs=0.55 #1????

colors='rgbcmyk'

for i,H in enumerate('ABCDEFG'):
  plt.loglog(PhRat[H][0], PhRat[H][1]*4**i, 'o-', color=colors[i])
  f=(PhRat[H][1]>0)&(PhRat[H][0]<1)&(PhRat[H][0]>0.1)
  p=np.polyfit(np.log(PhRat[H][0][f]), np.log(PhRat[H][1][f]), 1) 
  plt.plot(PhRat[H][0], np.exp(np.polyval(p, np.log(PhRat[H][0])))*4**i, 'k--')
  

plt.figure()
for i,H in enumerate('ABCDEFG'):
  plt.loglog(PhRat[H][0], PhRat[H][1]*PhHalo[H][1]*PhRv[H]**3*4**i, 'o-', color=colors[i])
  f=(PhRat[H][1]>0)&(PhRat[H][0]<1)&(PhRat[H][0]>0.1)
  p=np.polyfit(np.log(PhRat[H][0][f]), np.log(PhRat[H][1][f]), 1) 
  #beta=p[0]
  print beta
  fMin=1000*PhmP[H]/PhMv[H]
  rho0=fs*np.exp(sigma**2*alpha**2/2)*((fMin/Mstar)**(-alpha)-(0.1/Mstar)**(-alpha))/alpha*PhHalo[H][1]/PhMv[H]*PhRv[H]**3 # dN/d^3(R/Rv)
  rho=rho0*(PhRat[H][0]/Rstar)**(alpha*beta)
  plt.plot(PhRat[H][0], rho*4**i,'--', color=colors[i])

plt.ylim([1,1e7])  

#plt.figure()
#for i,H in enumerate('ABCDE'):
  #plt.loglog(AqRat[H][0], AqRat[H][1]*AqHalo[H][1]*AqRv[H]**3*4**i, 'o-', color=colors[i])
  #f=(AqRat[H][1]>0)&(AqRat[H][0]<1)&(AqRat[H][0]>0.1)
  #p=np.polyfit(np.log(AqRat[H][0][f]), np.log(AqRat[H][1][f]), 1) 
  ##beta=p[0]
  #print beta
  #fMin=1000*AqmP[H]/AqMv[H]
  #rho0=fs*np.exp(sigma**2*alpha**2/2)*((fMin/Mstar)**(-alpha)-(0.1/Mstar)**(-alpha))/alpha*AqHalo[H][1]/AqMv[H]*AqRv[H]**3 # dN/d^3(R/Rv)
  #rho=rho0*(AqRat[H][0]/Rstar)**(alpha*beta)
  #plt.plot(AqRat[H][0], rho*4**i,'--', color=colors[i])
