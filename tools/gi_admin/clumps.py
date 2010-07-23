#!/usr/bin/env python

from h5part import H5PartDump
import numpy as np

import matplotlib as mp
#mp.use('Agg')
import matplotlib.pyplot as plt

class ClumpsPlotConfig(object):
  def __init__(self):
    self.mmin = 0.0008
    self.mmax = 2.0

    self.ME = 5.97360e27

    self.noplotclumps = 2
    self.avgpts = 10
    self.plotesc = True


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
    self.nos    = nos
    self.maxnoc = maxnoc

    # get the mass derivative at half steps
    dmdth = ( m[1:,:] - m[0:-1,:] ) / ( t[1:] - t[0:-1] )

    # interpolate at full steps
    dmdt = np.zeros([nos,maxnoc])
    dmdt[1:-1,:] = 0.5*(dmdth[1:,:] + dmdth[0:-1,:])
    dmdt[0,:] = dmdt[1,:]
    dmdt[-1,:] = dmdt[-2,:]

    dnoc = noc[1:] - noc[:-1]

    self.dmdt = dmdt
    
  def __del__(self):
    pass

  def bodiesBound(self,nob):
    pass

  def systemBound(self):
    return True
  
  def __getStepNum(self,sname):
    return int(sname.replace("/Step#",""))

  def dmdtAboveTau(self,tau, mmin):
    taum = abs( (self.m / self.dmdt) )
    return ( np.isnan(taum) | (taum > tau) | (self.m < mmin) )

  def _getMeanMasses(self,nop):
    if nop > self.nos:
       nop = nos
    mmean = np.mean( self.m[-nop-1:-1,:] , axis=0)
    mvar  = np.sqrt( np.var( self.m[-nop-1:-1,:] , axis=0) )
    return (mmean, mvar)

  def plot(self,pname,cfg,tit):
    fig = plt.figure()
    fig.clear()
    ax = plt.axes()
    
    nopc = cfg.noplotclumps
    ME   = cfg.ME
    
    (mmean, mvar) = self._getMeanMasses(cfg.avgpts)

    annstr = ""
    for cuc in range(1,nopc+1):
      ax.semilogy(self.t, self.m[:,cuc] / ME, '-',ms=1.0)
      annstr += str(cuc) + ".   " + ( "%.3e" % (mmean[cuc] / ME) ) + " +/- " \
          +  ( "%.3e" % (mvar[cuc] / ME) ) + " Me\n"

    if cfg.plotesc:
      ax.semilogy(self.t, self.m[:,0] / ME, 'r-')
      annstr += "esc. " + ( "%.3e" % (mmean[0] / ME) ) + " +/- " \
          +  ( "%.3e" % (mvar[0] / ME) ) + " Me\n"

    ax.text(0.5, 0.4, annstr, transform = ax.transAxes)
    ax.set_ylim(cfg.mmin, cfg.mmax)
    ax.set_title(tit)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("clump mass [Me]")
    fig.savefig(pname)

