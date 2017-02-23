''' append mass-concentration using only the bound mass'''
import sys,os
sys.path.append(os.path.abspath('..'))
import numpy as np
import matplotlib.pyplot as plt
import h5py
from nfw import NFWHalo,mean_concentration
from MbdIO import cross_section
import emcee
from myutils import myhist

#plt.ion()

#rootdir='/gpfs/data/jvbq85/SubProf/'
rootdir='/work/Projects/SubProf/'
datadir=rootdir+'data/'
outdir=rootdir+'plots/'
outfile=h5py.File(datadir+'MockCluster.hdf5','r+')
common=outfile['/Common']
m,mAcc=common[1],common[2]
nsub=len(m)
offset=0.66 #cluster
#offset=1.29 #galaxy
#============================gen concentration==========================
sigmaC=0.3#log10-scatter=0.13, ln-scatter=0.3 
deltalnC=np.random.normal(0, sigmaC, nsub)
#=========================combine data===========

lncMaccio=np.log(mean_concentration(m,'MaccioW1'))
lncLudlow=np.log(mean_concentration(m,'Ludlow'))
def eval_luminosityBnd(lncbar):
  cBnd=np.exp(lncbar+deltalnC-offset)
  #============================evaluate truncation and luminosity=========
  halos=[NFWHalo(mAcc[i], cBnd[i]) for i in range(nsub)]
  rt=np.array([halos[i].radius(m[i]) for i in range(nsub)])
  Lt=np.array([halos[i].luminosity(rt[i]) for i in range(nsub)])
  Lv=np.array([halos[i].Ls for i in range(nsub)])
  return np.array([cBnd,rt,Lt,Lv])
Maccio=eval_luminosityBnd(lncMaccio)
Ludlow=eval_luminosityBnd(lncLudlow)

#============================save========================================
dset=outfile.create_dataset('MaccioShift',data=Maccio)
dset.attrs['rows']='cAcc,rt,Lt,Lv'
dset.attrs['sigma_lnC']=sigmaC

dset=outfile.create_dataset('LudlowShift',data=Ludlow)
dset.attrs['rows']='cAcc,rt,Lt,Lv'
dset.attrs['sigma_lnC']=sigmaC

outfile.close()