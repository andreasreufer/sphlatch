#!/usr/bin/env python

from h5part import H5PartDump
import numpy as np

class SimClumps(object):
  def __init__(self,fname):
    dump = H5PartDump(fname)
    
    snames = dump.getStepNames()
    snames.sort(key=self.__getStepNum)
    
    nos = len(snames)
    maxnoc = 0
    
    for sname in snames:
      step = dump.getStep(sname)
      maxnoc = max(maxnoc,step.m.shape[0])

    m = np.zeros([nos,maxnoc])
    t = np.zeros([nos,1])
    noc = np.zeros([nos])

    cs = 0
    for sname in snames:
      step = dump.getStep(sname)
      # number of clumps + escaping matter
      cnoc = step.m.shape[0]
      m[cs,0:cnoc] = step.m[:,0]
      t[cs] = float( step._v_attrs["time"] )
      noc[cs] = cnoc - 1
      
      cs += 1

    self.t      = t
    self.m      = m
    self.noc    = noc
    self.maxnoc = maxnoc

    # get the mass derivative at half steps
    dmdth = ( m[1:,:] - m[0:-1,:] ) / ( t[1:] - t[0:-1] )

    # interpolate at full steps
    dmdt = np.zeros([nos,maxnoc])
    dmdt[1:-1,:] = 0.5*(dmdth[1:,:] + dmdth[0:-1,:])

    self.dmdt = dmdt
    self.dnoc = noc[1:] - noc[:-1]

  def __del__(self):
    pass

  

  def bodiesBound(self,nob):
    pass

  def systemBound(self):
    return True
  
  def __getStepNum(self,sname):
    return int(sname.replace("/Step#",""))



