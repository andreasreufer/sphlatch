#!/usr/bin/env python

from h5part import H5PartDump
import numpy as np
np.seterr(all='ignore')

import matplotlib as mp
import matplotlib.pyplot as plt

mp.rc('text.latex', preamble = '\usepackage{amssymb}, \usepackage{wasysym}')

class Disk(object):
  def __getData(self, step):
    if step._v_pathname == '/Step#0':
      return

    i = self.snames.index(step._v_pathname)
    #self.Ltot_mat00[i,:,:] = self.disk001.Ltot_mat00[:,:]
    #self.Ltot_mat01[i,:,:] = self.disk001.Ltot_mat01[:,:]
    #self.Ltot_mat02[i,:,:] = self.disk001.Ltot_mat02[:,:]
    #self.Ltot_mat04[i,:,:] = self.disk001.Ltot_mat04[:,:]
    #self.Ltot_mat05[i,:,:] = self.disk001.Ltot_mat05[:,:]

    self.mimp_mat01[i] = sum( step.disk001.mimp_mat01[:] )
    self.mtar_mat01[i] = sum( step.disk001.mtar_mat01[:] )
    
    self.mimp_mat02[i] = sum( step.disk001.mimp_mat02[:] )
    self.mtar_mat02[i] = sum( step.disk001.mtar_mat02[:] )
    
    self.mimp_mat04[i] = sum( step.disk001.mimp_mat04[:] )
    self.mtar_mat04[i] = sum( step.disk001.mtar_mat04[:] )

    self.t[i] =  float( step._v_attrs["time"] )

  def __getStepNum(self,sname):
    return int( sname.replace("/Step#","").replace("/disk001",""))


  def __init__(self,fname,loadsat=False):
    dump = H5PartDump(fname)
    self.loadsat = loadsat
    
    snames = dump.getStepNames()
    snames.sort(key=self.__getStepNum)
    nos    = len(snames)
    self.snames = snames
    self.nos    = nos
 
    self.mimp_mat01   = np.zeros([nos])
    self.mtar_mat01   = np.zeros([nos])
    
    self.mimp_mat02   = np.zeros([nos])
    self.mtar_mat02   = np.zeros([nos])
    
    self.mimp_mat04   = np.zeros([nos])
    self.mtar_mat04   = np.zeros([nos])

    self.t            = np.zeros([nos])
    
    dump.forEachStep(self.__getData)
    return

    self.m   = np.zeros([nos,maxnoc])
    self.t   = np.zeros([nos,1])
    self.noc = np.zeros([nos])
    
    self.pos = np.zeros([nos,maxnoc,3])
    self.vel = np.zeros([nos,maxnoc,3])
    self.L   = np.zeros([nos,maxnoc,3])
    
    self.P   = np.zeros([nos,maxnoc,3])
    self.rho = np.zeros([nos,maxnoc])
    self.rc  = np.zeros([nos,maxnoc])

    dump.forEachStep(self.__getMasses)

    m = self.m
    t = self.t
    noc = self.noc
    
    self.dt = dt

    # get the mass derivative at half steps
    dmdth = ( m[1:,:] - m[0:-1,:] ) / dt

    # interpolate at full steps
    dmdt = np.zeros([nos,maxnoc])
    dmdt[1:-1,:] = 0.5*(dmdth[1:,:] + dmdth[0:-1,:])
    dmdt[0,:] = dmdt[1,:]
    dmdt[-1,:] = dmdt[-2,:]

    dnoc = noc[1:] - noc[:-1]

    self.dmdt = dmdt
    self.dnoc = dnoc

  def __del__(self):
    pass


maxmat  = 33

def getDisk(pdumpfn, sim):
  pdumpf = H5PartDump(pdumpfn)
  sname  = pdumpf.getStepNames()[0]
  pdump = pdumpf.getStep(sname)

  m       = pdump.m[:,0]
  mat     = pdump.mat[:,0]
  ecc     = pdump.ecc[:,0]
  orbit   = pdump.orbit[:,0]
  clumpid = pdump.clumpid[:,0]

  #maxnoc = clumpid.max()
  maxnoc = int(pdumpf.getAttr(sname,"noclumps"))

  isdisk = ( orbit == 2 )

  res = sim.results
  for cid in range(1,maxnoc+1):
    for mati in range(0,maxmat):
      mcur = m.compress( isdisk & (clumpid == cid) & (mat == mati) ).sum()
      res.mdiskmatm[cid,mati] = mcur
      res.mdiskmatv[cid,mati] = 0.
      
      res.ediskmatm[cid,mati] = ( m * ecc ).compress( isdisk & (clumpid == cid) & (mat == mati) ).sum() / mcur
      res.ediskmatv[cid,mati] = 0.
      
  
  pdumpf.close()
  return

