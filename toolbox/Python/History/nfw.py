import numpy as np
from scipy.optimize import fsolve

OmegaM=0.3
OmegaL=0.7
G=43.0071
H0=100. #km/s/(Mpc/h)
RhoCrit=3*H0**2/8/np.pi/G
#z=0.
#scaleF=1./(1+z);
#Hz=H0 * sqrt(OmegaM/scaleF**3+ (1 -OmegaM -OmegaL) / scaleF**2 +OmegaL)
#Hratio=Hz/HUBBLE0
#OmegaZ=OmegaM/scaleF**3./Hratio**2

def NFWFunc(x):
  return np.log(1+x)-x/(1+x)

class NFWHalo:
  '''NFW halo'''
  def __init__(self,m=100.,c=None,rhos=None,rs=None, DeltaC=200.):
	'''initialize '''
	self.M=m
	if c is None:
	  self.C=5.74*(self.M/2e2)**-0.097 #Duffy08, C200, at z=0, all
	  #self.C=6.67*(self.M/2e2)**-0.092 #relaxed
	  #self.C=5.71*(self.M/2e2)**-0.084/2**0.47 #z=1, all
	else:
	  self.C=c
	self.Rhos=DeltaC/3.*self.C**3/NFWFunc(self.C)*RhoCrit
	self.Rv=(self.M/DeltaC/RhoCrit/(4*np.pi/3))**(1./3)
	self.Rs=self.Rv/self.C
	if rhos is not None:
	  self.Rhos=rhos
	  self.Rs=rs
	  #TODO: complete this
	self.Ms=4*np.pi*self.Rhos*self.Rs**3
	self.Pots=4*np.pi*self.Rhos*self.Rs**2
	
  def mass(self,r):
	'''cumulative mass profile'''
	return self.Ms*NFWFunc(r/self.Rs)
  
  def density(self,r):
	'''density'''
	x=r/self.Rs
	return self.Rhos/x/(1+x)**2
  
  def density_cum(self,r):
	'''cumulative density, inside r'''
	x=r/self.Rs
	return 3*self.Rhos*NFWFunc(x)/x**3
  
  def strip_func(self, sat, r, k=1):
	''' m/m_0 for a subhalo (sat) inside the current halo'''
	x=np.array(r,ndmin=1)
	y=np.zeros_like(x)
	for i,R in enumerate(x):
	  rhs=(2+k)*self.density_cum(R)-3*self.density(R)
	  func=lambda a:np.log(sat.density_cum(np.exp(a)*sat.Rs))-np.log(rhs)
	  result=fsolve(func, 1.)
	  y[i]=sat.mass(np.exp(result[0])*sat.Rs)/sat.M
	if np.isscalar(r):
	  return y[0]
	return y
