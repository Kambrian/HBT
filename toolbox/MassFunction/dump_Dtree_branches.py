import numpy, h5py
import pp

def GetBranch(treefile, branchID):
  MaxSnap=67
  snap0=treefile['/treeIndex/finalSnapshot'][branchID]
  offset=treefile['/treeIndex/firstNode'][branchID]
  len=treefile['/treeIndex/numberOfNodes'][branchID]
  nbound=treefile['/haloTrees/particleNumber'][offset:offset+len]
  vmaxPhysical=treefile['/haloTrees/maximumCircularVelocity'][offset:offset+len]
  Nbound=0
  Vmax=0.
  IsCentral=0
  HostId=-1;
  if snap0==MaxSnap:
    Nbound=nbound[0]
    Vmax=vmaxPhysical[0]
    IsCentral=treefile['/haloTrees/isFoFCentre'][offset]
    HostId=treefile['/haloTrees/fofIndex'][offset]-1
  NboundPeak=nbound.max()
  VmaxPeak=vmaxPhysical.max()
  PID=treefile['/haloTrees/mostBoundID'][offset]
  return PID, IsCentral, Nbound, NboundPeak, Vmax, VmaxPeak, branchID, HostId
    
def DumpParticles(ifile):
  print ifile
  rootdir='/gpfs/data/jch/Millennium2/Trees/trees/treedir_067/'
  #if not os.path.exists(outdir):
      #os.mkdir(outdir)
  treefile=h5py.File(rootdir+'tree_067.%d.hdf5'%ifile,'r')
  nbranch=treefile['/treeIndex/firstNode'].shape[0]
  branches=numpy.array([GetBranch(treefile, i) for i in xrange(nbranch)], dtype=[('MostBoundID', 'i8'), ('IsCentral', 'i4'), ('Nbound', 'i4'), ('NboundPeak', 'i4'), ('Vmax', 'f4'), ('VmaxPeak', 'f4'), ('BranchID', 'i4'), ('HostID', 'i4')])
  treefile.close()

  outdir='/gpfs/data/jvbq85/HBT/data/Millennium2/subcat/analysis/DTrees/'
  outfile=h5py.File(outdir+'DTreeMostBoundParticles.%d.hdf5'%ifile, 'w')
  outfile.create_dataset('Branches', data=branches)
  outfile.close()
  return ifile

job_server = pp.Server(ncpus=12) 
jobs=[job_server.submit(DumpParticles, (i,), (GetBranch,), ('h5py', 'numpy')) for i in xrange(128)]
for job in jobs:
  job()

  
  