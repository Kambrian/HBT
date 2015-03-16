import numpy as np
import h5py
from scipy.interpolate import interp1d, UnivariateSpline

rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'


def fmtexp10(x, fmt='%1.0f'):
  ''' print x as a*10^b
  precision of a specified in fmt 
  '''
  j=np.floor(np.log10(x))
  i=x/10.**j
  if i==1:
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
  def __init__(self,halo='A2', variant='DTree'):
	if variant!=None:
	  mbdfile=h5py.File(datadir+halo+'Mbd.'+variant+'.hdf5','r')
	else:
	  mbdfile=h5py.File(datadir+halo+'Mbd.hdf5','r')
	self.name=halo
	self.mP=mbdfile['/PartMass'][...]
	self.m=mbdfile['/mass'][...]
	self.x=mbdfile['/x'][...]
	self.flag=mbdfile['/DirectInfall'][...]
	self.snapTVV=mbdfile['/snapTVV'][...]
	self.massTVV=mbdfile['/massTVV'][...]
	mbdfile.close()
	if not ((halo=='A4')&(variant==None)):
	  self.x*=1e3 #convert from Mpc/h to kpc/h
	self.r=np.sqrt(np.sum(self.x**2,1))
	
	cenid=np.argmax(self.m)
	print self.m[cenid]*self.mP/134.22
	#self.r=np.sqrt(np.sum((self.x-self.x[cenid])**2,1))#shift center to mstbnd
	self.rmin=self.r[cenid] #distance between mostbound particle and core of main sub
	self.r=np.delete(self.r, cenid)
	self.m=np.delete(self.m, cenid)
	self.x=np.delete(self.x, cenid, axis=0)
	self.flag=np.delete(self.flag, cenid, axis=0)
	self.snapTVV=np.delete(self.snapTVV,cenid,axis=0)
	self.massTVV=np.delete(self.massTVV,cenid,axis=0)
	
	if halo[0]=='A':
	  #C200 definition
	  self.Mvir=134.22
	  self.Rvir=179.38
	elif halo[0]=='B':
	  self.Mvir=60.9
	  self.Rvir=137.86
	try:
	  self.DM=np.loadtxt(datadir+halo+'density.cen_mstbnd.dat')
	except:
	  pass
	
  def get_sub_density(self,flags, xbin, Rref):
	'''x and Rref in units of Rvir'''
	xcen=xbin[1:]/np.sqrt(xbin[1]/xbin[0])
	vol=np.diff(xbin**3)
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
	sp_err=interp1d(np.log(self.DM[:,0]/self.Rvir), np.log(self.DM[:,2]), bounds_error=False)
	eval_sperr=lambda x: np.exp(sp_err(np.log(x)))
	density_err=eval_sperr(xcen)/densityRef
	return xcen,density,densityRef,density_err

def cross_section(x,y,xbin,q=[10,30,50,70,90]):
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
  return p.T,xmid