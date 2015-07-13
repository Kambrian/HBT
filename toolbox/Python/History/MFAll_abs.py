#absolute mass on x axis
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
iInfall=0

def getMF(H, xbin, iInfall, nMin=100):
  '''H: halodata'''
  x=xbin[:-1]*np.sqrt(xbin[1]/xbin[0])
  
  fclip=xbin[:-1]<nMin*H.mP #clip below nMin
  
  m=H.massTVV[:,iInfall]
  data=m[(H.r<Rref*H.Rvir)&(m>0)]*H.mP
  n,_=np.histogram(data, xbin)
  yInfall=n/np.log(xbin[1]/xbin[0])/H.Mvir
  yInfall[fclip]=np.nan
  
  f=(n>1)&(~fclip)
  parInfall=powerlaw_fit(x[f], yInfall[f], w=None)#1./np.sqrt(n)[f])
  #print parInfall[0]
				   
  m=H.m
  data=m[(H.r<Rref*H.Rvir)&(m>0)]*H.mP
  n,_=np.histogram(data, xbin)
  y=n/np.log(xbin[1]/xbin[0])/H.Mvir
  y[fclip]=np.nan
  
  f=(n>1)&(~fclip)
  try:
	par=powerlaw_fit(x[f], y[f], w=None)#1./np.sqrt(n)[f])
  except:
	par=[[np.nan,np.nan], [np.nan, np.nan], np.nan]
  #print par[0]
  
  return y,yInfall,x,par[0], parInfall[0]

def batchMF(halolist,xbin):
  '''process a list of haloes'''
  yall=[]
  yallInfall=[]
  par=[]
  parInfall=[]
  for h in halolist:
	y,y0,x,p,pi=getMF(HaloData(h),xbin,iInfall)
	yall.append(y*x)
	yallInfall.append(y0*x)
	par.append(p)
	parInfall.append(pi)
  yall=np.array(yall)
  yallInfall=np.array(yallInfall)
  yrat=yall/yallInfall
  par=np.array(par)
  parInfall=np.array(parInfall)
  return {'x':x,'y':yall,'yi':yallInfall,'yr':yrat,'p':par,'pi':parInfall}

def plotBand(x,data, label=None, color='r', ls='-', alpha=0.3, filt=None):
  y=np.nanmean(data,0)
  ystd=np.nanstd(data,0)
  f=(y>0)&(ystd>0)
  if filt!=None:
	f=f&filt
  l,=plt.loglog(x[f], y[f], linestyle=ls, color=color, label=label)
  plt.fill_between(x[f], (y+ystd)[f], (y-ystd)[f], color=color, alpha=0.3)
  return l

xbin=np.logspace(-5,2,20)
#xbin=np.logspace(-7,-2, 50)
AqA=batchMF(['AqA'+h for h in '12345'], xbin)
xbin=np.logspace(-5,0,30)
Aq=batchMF(['Aq'+h+'2' for h in 'ABCDEF'], xbin)
xbin=np.logspace(-2,3,30)
Ph=batchMF(['Ph'+h+'2' for h in 'ABCDEF'], xbin)

print 'AqPars:', Aq['p'].mean(0), Aq['p'].std(0)
print 'AqPars Infall:', Aq['pi'].mean(0), Aq['pi'].std(0)
print 'PhPars:', Ph['p'].mean(0), Ph['p'].std(0)
print 'PhPars Infall:', Ph['pi'].mean(0), Ph['pi'].std(0)

#convergence of Aq-A
fmt=['m','r','g','b','c']
plt.figure()
x=AqA['x']
for i in range(5):
  y=AqA['y'][i]
  plt.loglog(x[y>0], y[y>0], fmt[i]+'-', lw=i+1, label='AqA%d'%(i+1), alpha=1-i*0.2)
plt.plot(x, x*AqA['p'][0][2]*x**AqA['p'][0][1], 'k-', label=r'$\alpha=%.2f$'%(-AqA['p'][0][1]))  
#plt.ylim([1e-2,0.3])
plt.legend(loc=3, fontsize=15)
plt.xlabel(r'$m_{\rm acc}[10^{10}M_{\odot}/h]$')
plt.ylabel(r'$(m_{\rm acc}/M_{200})\mathrm{d}N/\mathrm{d}\ln m_{\rm acc}$')
plt.savefig(outdir+'MF.AqA.pdf')

plt.figure()
for i in range(5):
  y=AqA['yr'][i]
  plt.loglog(x[y>0], y[y>0], fmt[i]+'-', lw=i+1, label='AqA%d'%(i+1))
plt.legend(loc=3)
plt.xlabel(r'$m_{\rm acc}[10^{10}M_{\odot}/h]$')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}N_{\rm acc}$')
#plt.ylim([1e-2,0.3])
#plt.savefig(outdir+'MFRat.AqA.eps')

## Mass Functions
plt.figure()
##Aq
x=Aq['x']
plotBand(x, Aq['y'], 'Aquarius', 'r', '-')#, filt=(x>1e-6))
plotBand(x, Aq['yi'], None, 'r', '--')#, filt=(x>1e-6)) 
plt.plot(x, (x*Aq['pi'].mean(0)[2]*(x)**Aq['pi'].mean(0)[1]), 'k-')
plt.yscale('log', nonposy='clip')
##Ph
x=Ph['x']
plotBand(x, Ph['y'], 'Phoenix', 'g', '-', 0.2)#, filt=(x>1e-6))  
plotBand(x, Ph['yi'], None, 'g', '--', 0.2)#, filt=(x>1e-6)) 
plt.plot(x, (x*Ph['pi'].mean(0)[2]*x**Ph['pi'].mean(0)[1]), 'k-')
plt.yscale('log', nonposy='clip')

legend=plt.legend(loc=1,fontsize=12)

plt.xlabel(r'$m[10^{10}M_{\odot}/h]$')
plt.ylabel(r'$(m/M_{200})\mathrm{d}N/\mathrm{d}\ln m$')
#plt.savefig(outdir+'MF.All.pdf')

plt.figure()
plotBand(Aq['x'], Aq['yr'], 'Aquarius', 'r', '-', filt=None)
plotBand(Ph['x'], Ph['yr'], 'Phoenix', 'g', '-', 0.2, filt=None)
plt.yscale('log', nonposy='clip')
plt.legend(loc=3)
plt.xlabel(r'$m[10^{10}M_{\odot}/h]$')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}N_{\rm acc}$')
#plt.ylim([1e-2,0.3])

#plt.savefig(outdir+'MFRat.All.pdf')
x=np.hstack([Aq['x'],Ph['x']])
y=np.hstack([Aq['yi'].mean(0), Ph['yi'].mean(0)])
parsInfall=powerlaw_fit(x[y>0],(y/x)[y>0],)
#parsInfall[0][0]-=0.002
plt.figure()
plotBand(Aq['x'], Aq['yi'], 'Aquarius', 'r', '--', 0.2)
plt.plot(Aq['x'], 0.089*Aq['x']**(1-0.95), 'k-')
plotBand(Ph['x'], Ph['yi'], 'Phoenix', 'g', '--', 0.2)
plt.plot(Ph['x'], 0.080*Ph['x']**(1-0.95), 'k-')
plt.plot(x, x*(x/(parsInfall[0][0]))**parsInfall[0][1], 'k:',label=r'$\alpha=%.2f$'%(-parsInfall[0][1]))

y=np.hstack([Aq['y'].mean(0), Ph['y'].mean(0)])
a=np.nanmean(y/x**(1+parsInfall[0][1]));pars=[[np.nan,parsInfall[0][1],a],];
pars=powerlaw_fit(x[y>0],(y/x)[y>0])
plotBand(Aq['x'], Aq['y'], None, 'r', '--', 0.2)
plotBand(Ph['x'], Ph['y'], None, 'g', '--', 0.2)
plt.plot(Aq['x'], 0.0077*Aq['x']**(1-0.95), 'k-')
plt.plot(Ph['x'], 0.0080*Ph['x']**(1-0.95), 'k-')
plt.plot(x, x*pars[0][2]*x**pars[0][1], 'k:')#, label=r'$\alpha=%.2f$'%(-pars[0][1]))
plt.yscale('log', nonposy='clip')

plt.legend(loc=3)
plt.xlabel(r'$m[10^{10}M_{\odot}/h]$')
plt.ylabel(r'$(m/M_{200})\mathrm{d}N/\mathrm{d}\ln m$')
#plt.savefig(outdir+'MFabs.All2.pdf')

#plt.plot(x, x*AqA['pi'][0][2]*x**AqA['pi'][0][1], 'ko', label=r'$\alpha=%.2f$'%(-AqA['pi'][0][1]))  
#plt.plot(x, 0.073*x**-0.96*x, 'k:',label=r'$\alpha=%.2f$'%(0.96))