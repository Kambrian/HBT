import numpy as np
import h5py
from scipy.interpolate import interp1d, UnivariateSpline

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'

def powerlaw_fit(x,y, w=None):
  '''y=(x/x0)^beta
  or y=yp*(x/xp)^beta, where xp is the minimum variance point of y,
  or y=A*x^beta.
  w is the weight on y, e.g., 1/sigma_y
  pivotfit: pivot or amp; if pivot, fit (x/x0)^beta; otherwise fit A*x^beta'''
  f=(~np.isnan(y))&(~np.isnan(x))
  weight=None
  if w is not None:
	weight=(w*y)[f]
  par,V=np.polyfit(np.log(x)[f], np.log(y)[f], 1, w=weight, cov=True)
  beta=par[0]
  x0=np.exp(-par[1]/par[0])
  ebeta=np.sqrt(V[0,0])
  ex0=x0*np.sqrt(V[1,1]/par[0]**2+V[0,0]*par[1]**2/par[0]**4-2*par[1]/par[0]**3*V[0,1])
  xp=np.exp(-V[0,1]/V[0,0])
  yp=(xp/x0)**beta
  A=np.exp(par[1])
  eA=A*np.sqrt(V[1,1])
  return [x0,beta,A],[ex0,ebeta,eA],[xp,yp]

def fmtexp10(x, fmt='%1.0f'):
  ''' print x as a*10^b
  precision of a specified in fmt 
  '''
  j=np.floor(np.log10(x))
  i=x/10.**j
  if float(fmt%i)==1:
	s=r'10^{%d}'%j
  else:
	s=fmt%i+r'\times 10^{%d}'%j
  return s

def einasto(x, alpha=0.16):
  '''Einasto Profile: rho=rhos*exp{-2/alpha*[(r/rs)^alpha-1]},
  % with typical alpha~0.16
  % input: x=r/rs
  % output: rho/rhos
  % rhos, rs usually quoted as rho_{-2}, r_{-2}.
  '''
  return np.exp(-2./alpha*(x**alpha-1.))

class HaloData(object):
  def __init__(self,halo='AqA2', variant='DTree'):
	if variant!=None:
	  if halo[:2]=='Ph':
		mbdfile=h5py.File(datadir+'/Phoenix/'+halo+'.'+variant+'.hdf5','r')
	  else:
		mbdfile=h5py.File(datadir+'/Aquarius/'+halo+'Mbd.'+variant+'.hdf5','r')
	else:
	  mbdfile=h5py.File(datadir+'/Aquarius/'+halo+'Mbd.hdf5','r')
	self.name=halo
	self.mP=mbdfile['/PartMass'][...]
	self.m=mbdfile['/mass'][...]
	self.x=mbdfile['/x'][...]
	self.flag=mbdfile['/DirectInfall'][...]
	self.snapTVV=mbdfile['/snapTVV'][...]
	self.massTVV=mbdfile['/massTVV'][...]
	mbdfile.close()
	if not ((halo=='AqA4')&(variant==None)):
	  self.x*=1e3 #convert from Mpc/h to kpc/h
	self.r=np.sqrt(np.sum(self.x**2,1))
	
	cenid=np.argmax(self.m)
	Mcen=self.m[cenid]*self.mP
	#self.r=np.sqrt(np.sum((self.x-self.x[cenid])**2,1))#shift center to mstbnd
	self.rmin=self.r[cenid] #distance between mostbound particle and core of main sub
	self.r=np.delete(self.r, cenid)
	self.m=np.delete(self.m, cenid)
	self.x=np.delete(self.x, cenid, axis=0)
	self.flag=np.delete(self.flag, cenid, axis=0)
	self.snapTVV=np.delete(self.snapTVV,cenid,axis=0)
	self.massTVV=np.delete(self.massTVV,cenid,axis=0)
	HaloVir={'AqA':[134.3,179.41, 16.2], 'AqB':[59.82, 137.02, 9.7],'AqC':[129.5,177.26, 15.2], 'AqD':[129.5,177.28, 9.4],'AqE':[86.52,154.96, 8.3],'AqF':[82.82,152.72,9.8],
		  'PhA':[6.570e4, 1414., 6.0],'PhB':[8.255e4,1526.,4.2],'PhC':[5.495e4,1332.,5.1],'PhD':[6.191e4,1386.,4.1],'PhE':[5.969e4,1369.,5.2],'PhF':[7.997e4,1509.,4.6],'PhG':[1.150e5,1704.,3.3],'PhH':[1.136e5, 1686.,4.7], 'PhI':[2.411e5, 2185.,4.9]} #C200 definition
	#HaloVir={'AqA':[134.3,179.41, 16.10], 'AqB':[59.82, 137.02, 8.16],'AqC':[129.5,177.26, 12.34], 'AqD':[129.5,177.28, 8.73],'AqE':[86.52,154.96, 8.67],'AqF':[82.82,152.72,9.8],
		  #'PhA':[6.570e4, 1414., 6.0],'PhB':[8.255e4,1526.,4.2],'PhC':[5.495e4,1332.,5.1],'PhD':[6.191e4,1386.,4.1],'PhE':[5.969e4,1369.,5.2],'PhF':[7.997e4,1509.,4.6],'PhG':[1.150e5,1704.,3.3],'PhH':[1.136e5, 1686.,4.7], 'PhI':[2.411e5, 2185.,4.9]} #C200 definition
	#HaloVir={'AqA':[134.3,179.41, 14.13], 'AqB':[59.82, 137.02, 8.04],'AqC':[129.5,177.26, 12.10], 'AqD':[129.5,177.28, 8.47],'AqE':[86.52,154.96, 8.44],'AqF':[82.82,152.72,9.8],
		  #'PhA':[6.570e4, 1414., 6.0],'PhB':[8.255e4,1526.,4.2],'PhC':[5.495e4,1332.,5.1],'PhD':[6.191e4,1386.,4.1],'PhE':[5.969e4,1369.,5.2],'PhF':[7.997e4,1509.,4.6],'PhG':[1.150e5,1704.,3.3],'PhH':[1.136e5, 1686.,4.7], 'PhI':[2.411e5, 2185.,4.9]} #C200 definition
	self.Mvir=HaloVir[halo[:3]][0]
	self.Rvir=HaloVir[halo[:3]][1]
	self.Cvir=HaloVir[halo[:3]][2]
	print Mcen/self.Mvir
	try:
	  if halo[:2]=='Ph':
		self.DM=np.loadtxt(datadir+'/Phoenix/'+halo+'Halo.dat')
	  else:
		self.DM=np.loadtxt(datadir+'/Aquarius/'+halo+'Halo.dat')
	except:
	  try:
		self.DM=np.loadtxt(datadir+'/Aquarius/'+halo+'density.cen_mstbnd.dat')
	  except:
		print "WARNING: No DM profile found for", halo
		pass
	
  def get_sub_density(self,flags, xbin, Rref):
	'''x and Rref in units of Rvir'''
	xcen=xbin[1:]/np.sqrt(xbin[1]/xbin[0])
	vol=np.diff(xbin**3)*4.*np.pi/3.
	count,tmp=np.histogram(self.r[flags]/self.Rvir, xbin)
	density=count/vol
	density_err=np.sqrt(count)/vol
	sp=UnivariateSpline(xcen[count>0], density[count>0], w=1/density_err[count>0], k=1)
	densityRef=sp(Rref)
	density/=densityRef
	density_err/=densityRef
	#density=sp(xcen)/densityRef
	return xcen,density,densityRef,density_err
  
  def get_host_density(self, xbin, Rref):
	'''x and Rref in units of Rvir'''
	xcen=xbin[1:]/np.sqrt(xbin[1]/xbin[0])
	x=self.DM[:,0]/self.Rvir
	y=self.DM[:,1]
	#sp=UnivariateSpline(np.log(self.DM[:,0]/self.Rvir), np.log(self.DM[:,1]), k=1, s=0.)
	sp=interp1d(np.log(self.DM[:,0]/self.Rvir), np.log(self.DM[:,1]), bounds_error=False)
	eval_sp=lambda x: np.exp(sp(np.log(x)))
	densityRef=eval_sp(Rref)
	density=eval_sp(xcen)/densityRef
	if self.DM.shape[1]==3:
	  sp_err=interp1d(np.log(self.DM[:,0]/self.Rvir), np.log(self.DM[:,2]), bounds_error=False)
	  eval_sperr=lambda x: np.exp(sp_err(np.log(x)))
	  density_err=eval_sperr(xcen)/densityRef
	else:
	  density_err=0.
	return xcen,density,densityRef,density_err

  def getMF(self, xbin=None, iInfall=0, nMin=100, Rref=1.):
	'''H: halodata'''
  
	if xbin is None:
	  xbin=np.logspace(np.log10(20*self.mP), np.log10(self.Mvir/100.), 30)
	  
	x=xbin[:-1]*np.sqrt(xbin[1]/xbin[0])
	
	fclip=xbin[:-1]<nMin*self.mP #clip below nMin
	
	m=self.massTVV[:,iInfall]
	data=m[(self.r<Rref*self.Rvir)]*self.mP
	n,_=np.histogram(data, xbin)
	yInfall=n/np.log(xbin[1]/xbin[0])/self.Mvir
	yInfall[fclip]=np.nan
	
	f=(n>1)&(~fclip)
	parInfall=powerlaw_fit(x[f], yInfall[f], w=None)#1./np.sqrt(n)[f])
	#print parInfall[0]
					
	m=self.m
	data=m[(self.r<Rref*self.Rvir)&(m>0)]*self.mP
	n,_=np.histogram(data, xbin)
	y=n/np.log(xbin[1]/xbin[0])/self.Mvir
	y[fclip]=np.nan
	
	f=(n>1)&(~fclip)
	try:
	  par=powerlaw_fit(x[f], y[f], w=None)#1./np.sqrt(n)[f])
	except:
	  par=[[np.nan,np.nan], [np.nan, np.nan], np.nan]
	#print par[0]
	
	return y,yInfall,x,par[0], parInfall[0]

  def get_surv_rate(self, iInfall=0, nMinInfall=10000, rmin=0.5, rmax=0.8):
	f=(self.massTVV[:,iInfall]>nMinInfall)&(self.r/self.Rvir>rmin)&(self.r/self.Rvir<rmax)
	SurvRate=(f&(self.m>0)).sum()/f.sum()
	return SurvRate

  def get_prof_rat(self,xbin=None,fMin=None,Rref=1.):
	if fMin is None:
	  fMin=1000.*self.mP/self.Mvir
	if xbin is None:
	  xbin=np.logspace(-3,np.log10(2), 30)
	  
	rHalo,denHalo,denHaloRef,denHaloErr=self.get_host_density(xbin, Rref)
	rSub,denSub,denSubRef,denSubErr=self.get_sub_density((self.m>fMin*self.Mvir/self.mP), xbin, Rref)
	plt.plot(rSub[denSub>0], (denSub/denHalo)[denSub>0], 'ro', label=self.name)
	f=(rSub<Rref)&(denSub>0)
	pars=powerlaw_fit(rSub[f], (denSub/denHalo)[f])
	return rSub, denSub/denHalo, [pars[0][1], pars[1][1]]
  
  def get_pars(self):
	ParMF,ParMFInfall=self.getMF()[3:]
	SurvFrac=self.get_surv_rate()
	ParProf=self.get_prof_rat()[2]
	
	
def cross_section(x,y,xbin,q=[10,30,50,70,90]):
  ''' return cross section (y percentiles specified by q) at xbin '''
  count,xbin=np.histogram(x,xbin)
  nbin=len(xbin)-1
  bin=np.digitize(x,xbin)-1
  p=np.zeros([nbin, len(q)])
  xmid=np.zeros(nbin)
  for i in xrange(nbin):
	data=y[bin==i]
	if len(data):
	  p[i]=np.percentile(data,q)
	  xmid[i]=np.median(x[bin==i])
	else:
	  p[i]=np.nan+np.zeros(len(q))
	  xmid[i]=(xbin[i]+xbin[i+1])/2.
  #p=np.array([np.percentile(y[bin==i],q) for i in xrange(nbin)])
  #xmid=np.array([np.median(x[bin==i]) for i in xrange(nbin)])
  #t=xbin[:-1]
  #h1=plt.plot(t, p[:,0], 'k-', label='Median')
  #h2=plt.plot(t, p[:,1:3], 'r:', label='25;75')
  #h3=plt.plot(t, p[:,3:], 'g--', label='10;90')
  return p.T,xmid #p.T of shape [npercentile, nbin]