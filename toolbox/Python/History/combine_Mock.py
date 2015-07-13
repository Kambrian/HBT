''' generate Monte-Carlo samples of subhaloes '''
import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import h5py

rootdir='/gpfs/data/jvbq85/SubProf/'
#rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'

halo=sys.argv[1] #Cluster or Galaxy

d1=[]
d2=[]
d3=[]
d4=[]
for ithread in xrange(1,6):
  infile=h5py.File(datadir+'Mock'+halo+'-%d.hdf5'%ithread,'r')
  m=infile['MHost'][...]
  sigmaC=infile['Maccio'].attrs['sigma_lnC']
  rowsCommon=infile['Common'].attrs['rows']
  rowsMaccio=infile['Maccio'].attrs['rows']
  d1.append(infile['Common'][...])
  d2.append(infile['Maccio'][...])
  d3.append(infile['Ludlow'][...])
  d4.append(infile['Ludlow1'][...])
  infile.close()

outfile=h5py.File(datadir+'Mock'+halo+'New.hdf5','w')
outfile.create_dataset('MHost',data=m)
comm=np.hstack(d1)
comm[0]=comm[0]/5. #stack haloes, renormalize
dset=outfile.create_dataset('Common', data=comm)
dset.attrs['rows']=rowsCommon

dset=outfile.create_dataset('Maccio',data=np.hstack(d2))
dset.attrs['rows']=rowsMaccio
dset.attrs['sigma_lnC']=sigmaC
dset=outfile.create_dataset('Ludlow',data=np.hstack(d3))
dset=outfile.create_dataset('Ludlow1',data=np.hstack(d4))

outfile.close()
