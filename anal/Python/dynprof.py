from numpy import *

def spherical_basis(r):
  """ get the spherical basis vector at radial vector r[n,3]
  return: er: radial;
	  et: azimuthal,(0~2*pi)
	  ef: elevation (0~pi)
  """ 
  dxy2=sum(r[:2]**2,1)
  dr=sqrt(sum(r**2,1))
  er=r/dr
  er[dr==0]=[1.,0.,0.] #arbitrary
  modet=sqrt(dxy2)
  et=zeros_like(r)
  et[:,0]=-r[:,1]/modet
  et[:,1]=r[:,0]/modet
  et[modet==0]=[1.,0.,0.]
  ef=zeros_like(r)
  modef=sqrt(dxy2+(dxy2/r[:,2])**2)
  ef[:,0]=-r[:,0]/modef
  ef[:,1]=-r[:,1]/modef
  ef[:,2]=dxy2/r[:,2]/modef
  ef[modet==0]=[0.,1.,0.]
  ef[r[:,2]==0]=[0.,0.,1.]
  return er,et,ef

import h5py
f=h5py.File('catalogue_bhb.hdf5','r')
x=f['/Stars/Coordinates'][:]
v=f['/Stars/Velocity'][:]
p=f['/Stars/Potential'][:]
idmin=argmin(p)
er,et,ef=spherical_basis(x-x[idmin])
vr=sum(v*er,axis=1)
vt=sum(v*et,1)
vf=sum(v*ef,1)
r=sqrt(sum((x-x[idmin])**2,1))
rbins=logspace(-2,1,20)
n,tmp=histogram(r,rbins)
vrsum,tmp=histogram(r,rbins,weights=vr)
vtsum,tmp=histogram(r,rbins,weights=vt)
vfsum,tmp=histogram(r,rbins,weights=vf)
vr2sum,tmp=histogram(r,rbins,weights=vr*vr)
vt2sum,tmp=histogram(r,rbins,weights=vt*vt)
vf2sum,tmp=histogram(r,rbins,weights=vf*vf)
sig2r=vr2sum/n-(vrsum/n)**2
sig2t=vt2sum/n-(vtsum/n)**2
sig2f=vf2sum/n-(vfsum/n)**2
beta=1-(sig2t+sig2f)/2/sig2r
semilogy(rbins[1:],beta,'.')