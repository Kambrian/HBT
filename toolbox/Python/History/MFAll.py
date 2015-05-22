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
  

##plot profiles
Rref=1
nbin=20
iInfall=0

def getMF(H, xbin, iInfall, nMin=100):
  '''H: halodata'''
  x=xbin[:-1]*np.sqrt(xbin[1]/xbin[0])
  
  fclip=xbin[:-1]<nMin*H.mP/H.Mvir #clip below nMin
  
  m=H.massTVV[:,iInfall]
  data=m[(H.r<Rref*H.Rvir)]*H.mP/H.Mvir
  n,_=np.histogram(data, xbin)
  yInfall=n/np.log(xbin[1]/xbin[0])
  yInfall[fclip]=np.nan
  
  f=(n>4)&(~fclip)
  parInfall=powerlaw_fit(x[f], yInfall[f], w=None)#1./np.sqrt(n)[f])
  #print parInfall[0]
				   
  m=H.m
  data=m[(H.r<Rref*H.Rvir)&(m>0)]*H.mP/H.Mvir
  n,_=np.histogram(data, xbin)
  y=n/np.log(xbin[1]/xbin[0])
  y[fclip]=np.nan
  
  f=(n>4)&(~fclip)
  try:
	par=powerlaw_fit(x[f], y[f], w=None)#1./np.sqrt(n)[f])
  except:
	par=[None, None]
  #print par[0]
  
  return y,yInfall,x,par[0],parInfall[0]

def plotBand(x,data, label=None, color='r', ls='-', alpha=0.3, filt=None):
  y=np.nanmean(data,0)
  ystd=np.nanstd(data,0)
  f=(y>0)&(ystd>0)
  if filt!=None:
	f=f&filt
  l,=plt.loglog(x[f], y[f], linestyle=ls, color=color, label=label)
  plt.fill_between(x[f], (y+ystd)[f], (y-ystd)[f], color=color, alpha=0.3)
  return l

#convergence of Aq-A
yall=[]
yallInfall=[]
xall=[]
xbin=np.logspace(-7,-3, 20)
#xbin=np.logspace(-7,-2, 50)
parInfall=[]
for h in '12345':
  H=HaloData('AqA'+h)
  y,y0,x,p,pi=getMF(H,xbin,iInfall, nMin=100)
  yall.append(y*x)
  yallInfall.append(y0*x)
  xall.append(x)
  parInfall.append(pi)
yall=np.array(yall)
yallInfall=np.array(yallInfall)
yrat=yall/yallInfall
#print parInfall

fmt=['m','r','g','b','c']
plt.figure()
for i in range(5):
  y=yallInfall[i]
  plt.loglog(xall[i][y>0], y[y>0], fmt[i]+'-', lw=i+1, label='AqA%d'%(i+1), alpha=1-i*0.2)
plt.plot(x, x*(x/parInfall[0][0])**parInfall[0][1], 'k', label=r'$\alpha=%.2f$'%(-parInfall[0][1]))  
plt.ylim([1e-2,0.3])
plt.legend(loc=3, fontsize=15)
plt.xlabel(r'$m_{\rm acc}/M_{200}$')
plt.ylabel(r'$(m_{\rm acc}/M_{200})\mathrm{d}N/\mathrm{d}\ln m_{\rm acc}$')
#plt.savefig(outdir+'MFInfall.AqA.pdf')

plt.figure()
for i in range(5):
  y=yrat[i]
  plt.loglog(xall[i][y>0], y[y>0], fmt[i]+'-', lw=i+1, label='AqA%d'%(i+1))
plt.legend(loc=3)
plt.xlabel(r'$m/M_{200}$')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}N_{\rm acc}$')
plt.ylim([1e-2,0.3])
#plt.savefig(outdir+'MFRat.AqA.eps')

## Mass Functions
plt.figure()
xbin=np.logspace(-7,-2,nbin)
##Aq
yall=[]
yallInfall=[]
par=[]
parInfall=[]
for h in 'ABCDEF':
  H=HaloData('Aq'+h+'2')
  y,y0,x,p,pi=getMF(H,xbin,iInfall)
  yall.append(y*x)
  yallInfall.append(y0*x)
  par.append(p)
  parInfall.append(pi)
yall=np.array(yall)
yallInfall=np.array(yallInfall)
yratAq=yall/yallInfall
par=np.array(par)
parInfall=np.array(parInfall)
print 'AqPars:', par.mean(0), par.std(0)
print 'AqPars Infall:', parInfall.mean(0), parInfall.std(0)

plotBand(x, yall, 'Aquarius', 'r', '-', filt=(x>1e-6))
plotBand(x, yallInfall, None, 'r', '--', filt=(x>1e-6)) 
plt.plot(x[x>1e-6], (x*(x/parInfall.mean(0)[0])**parInfall.mean(0)[1])[x>1e-6], 'k-')
plt.yscale('log', nonposy='clip')

##Ph
yall=[]
yallInfall=[]
par=[]
parInfall=[]
for h in 'ABCDEFGI':
  H=HaloData('Ph'+h+'2')
  y,y0,x,p,pi=getMF(H,xbin,iInfall)
  yall.append(y*x)
  yallInfall.append(y0*x)
  par.append(p)
  parInfall.append(pi)
yall=np.array(yall)
yallInfall=np.array(yallInfall)
yratPh=yall/yallInfall
par=np.array(par)
parInfall=np.array(parInfall)
print 'PhPars:', par.mean(0), par.std(0)
print 'PhPars Infall:', parInfall.mean(0), parInfall.std(0)

plotBand(x, yall, 'Phoenix', 'g', '-', 0.2, filt=(x>1e-6))  
plotBand(x, yallInfall, None, 'g', '--', 0.2, filt=(x>1e-6)) 
plt.plot(x[x>1e-6], (x*(x/parInfall.mean(0)[0])**parInfall.mean(0)[1])[x>1e-6], 'k-')
plt.yscale('log', nonposy='clip')

legend=plt.legend(loc=1,fontsize=12)

plt.xlabel(r'$m/M_{200}$')
plt.ylabel(r'$(m/M_{200})\mathrm{d}N/\mathrm{d}\ln m$')
#plt.savefig(outdir+'MF.All.pdf')

plt.figure()
plotBand(x, yratAq, 'Aquarius', 'r', '-', filt=(x>1e-6))
plotBand(x, yratPh, 'Phoenix', 'g', '-', 0.2, filt=(x>1e-6))
plt.yscale('log', nonposy='clip')
plt.legend(loc=3)
plt.xlabel(r'$m/M_{200}$')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}N_{\rm acc}$')
plt.ylim([1e-2,0.3])

#plt.savefig(outdir+'MFRat.All.pdf')
