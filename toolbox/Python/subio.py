from math import *
import numpy as np
import ctypes,os,ConfigParser,h5py

HBTInt=ctypes.c_longlong
HBTReal=ctypes.c_double
#=============C complex datatypes=====================
HBTxyz=HBTReal*3
HBTInt_p=ctypes.POINTER(HBTInt)
HBTReal_p=ctypes.POINTER(HBTReal)
HBTxyz_p=ctypes.POINTER(HBTxyz)

class ParticleData(ctypes.Structure):
  _fields_=[('PID', HBTInt_p),
	    ('Pos', HBTxyz_p),
	    ('Vel', HBTxyz_p),
	    ('Nsnap', HBTInt)
	    ]  


Particle_p=ctypes.POINTER(Particle_t)  

class Tracer_t(ctypes.Structure):
  pass
Tracer_p=ctypes.POINTER(Tracer_t)
Tracer_t._fields_=[('lnL', ctypes.c_double),
				   ('nP', ctypes.c_int),
				   ('mP', ctypes.c_double),
				   ('P', Particle_p),
				   ('nbin_r', ctypes.c_int),
				   ('FlagRLogBin', ctypes.c_int),
				   ('RadialCount', ctypes.POINTER(ctypes.c_int)),
				   ('rmin', ctypes.c_double),
				   ('rmax', ctypes.c_double),
				   ('nView', ctypes.c_int),
				   ('ViewType', ctypes.c_char),
				   ('Views', Tracer_p)
				  ]
	
class NFWHalo_t(ctypes.Structure):
  """all properties are physical"""
  _fields_=[('z', ctypes.c_double),
	    ('M', ctypes.c_double),
	    ('c', ctypes.c_double),
	    ('Rv', ctypes.c_double), 
	    ('Rs', ctypes.c_double),
	    ('Rhos', ctypes.c_double),
	    ('Pots', ctypes.c_double),#-4*pi*G*rhos*rs^2, the potential at r=0
	    ('Ms', ctypes.c_double), #4*pi*rhos*rs^3
	    ('virtype', ctypes.c_int)
	    ]  

#=======================load the library==========================
lib=ctypes.CDLL("../libdyn.so")
#general
lib.MaxNPar=10
lib.ParType=ctypes.c_double*lib.MaxNPar
lib.MODEL_TOL_BIN=ctypes.c_double.in_dll(lib,'MODEL_TOL_BIN')
lib.MODEL_TOL_BIN_ABS=ctypes.c_double.in_dll(lib,'MODEL_TOL_BIN_ABS')
lib.MODEL_TOL_REL=ctypes.c_double.in_dll(lib,'MODEL_TOL_REL')
lib.SubSampleSize=ctypes.c_int.in_dll(lib,'SubSampleSize')
lib.NumRadialCountBin=ctypes.c_int.in_dll(lib,'NumRadialCountBin')
lib.alloc_integration_space.restype=None
lib.alloc_integration_space.argtypes=[]
lib.free_integration_space.restype=None
lib.free_integration_space.argtypes=[]
lib.like_to_chi2.restype=ctypes.c_double
lib.like_to_chi2.argtypes=[ctypes.c_double, ctypes.c_int]

class Catalog:
	def __init__(self):
		pass

class GroupCatalog(Catalog):
	def __init__(self):
		Catalog.__init__(self)
		GroupCatalog.Ngroups=0
		
class SubCatalog(Catalog):
