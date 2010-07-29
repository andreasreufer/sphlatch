#!/usr/bin/env python

from h5part import H5PartDump
import numpy as np

import matplotlib as mp
#mp.use('Agg')
import matplotlib.pyplot as plt

class ClumpsPlotConfig(object):
  def __init__(self):
    self.mmin = 0.9e-5
    self.mmax = 2.5
    self.tmin =  -5.
    self.tmax =  55.

    self.ME = 5.97360e27

    self.noplotclumps =  8
    self.avgpts       = 10
    self.plotesc      = True


class SimClumps(object):
  def __countStepsClumps(self, step):
    self.nos += 1
    self.maxnoc = max(self.maxnoc, step.m.shape[0])
  
  def __getMasses(self, step):
    i = self.snames.index(step._v_pathname)
    cnoc = step.m.shape[0]
    self.m[i,0:cnoc] = step.m[:,0]
    self.t[i] = float( step._v_attrs["time"] )
    self.noc[i] = cnoc - 1
    self.pos[i,0:cnoc,:] = step.pos[:,:]
  
  def __getStepNum(self,sname):
    return int(sname.replace("/Step#",""))

  def __init__(self,fname):
    dump = H5PartDump(fname)
    
    self.nos = 0
    self.maxnoc = 0
    dump.forEachStep(self.__countStepsClumps)
    
    snames = dump.getStepNames()
    snames.sort(key=self.__getStepNum)
    self.snames = snames
    
    nos = self.nos
    maxnoc = self.maxnoc

    self.m   = np.zeros([nos,maxnoc])
    self.t   = np.zeros([nos,1])
    self.noc = np.zeros([nos])
    
    self.pos = np.zeros([nos,maxnoc,3])

    dump.forEachStep(self.__getMasses)

    m = self.m
    t = self.t
    noc = self.noc

    # get the mass derivative at half steps
    dmdth = ( m[1:,:] - m[0:-1,:] ) / ( t[1:] - t[0:-1] )

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

  def bodiesBound(self,nob):
    pass

  def systemBound(self,noclmps):
    return True
  
  def _getMeanMasses(self,nop):
    if nop > self.nos:
       nop = self.nos
    mmean = np.mean( self.m[-nop-1:-1,:] , axis=0)
    mvar  = np.sqrt( np.var( self.m[-nop-1:-1,:] , axis=0) )
    return (mmean, mvar)

  def _getExtents(self,noc):
    pass

  def _getdnopdtcol(self,sim):
    tcol   = sim.tcol
    mpp    = sim.mtot / sim.nop
    self.dnopdtcol = self.dmdt * tcol / mpp

  def plotTrajectories(self,pname,cfg,sim):
    fig = plt.figure()
    fig.clear()
    ax = plt.axes()

    simkey = sim.params.key
    
    nopc = min( cfg.noplotclumps, self.maxnoc - 1) 
    for cuc in range(1,nopc+1):
      ax.plot(self.pos[:,cuc,0], self.pos[:,cuc,1], '-')
    
    ax.set_title(simkey)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.savefig(pname)

  def plotMass(self,pname,cfg,sim):
    fig = plt.figure()
    fig.clear()

    ax1 = plt.subplot(211)

    simkey = sim.params.key
    tcol   = sim.tcol
    
    self._getdnopdtcol(sim)

    nopc = min( cfg.noplotclumps, self.maxnoc - 1) 
    ME   = cfg.ME
    
    (mmean, mvar) = self._getMeanMasses(cfg.avgpts)

    annstr = ""
    for cuc in range(1,nopc+1):
      ax1.semilogy(self.t / tcol, self.m[:,cuc] / ME, '-',ms=1.0)
      annstr += str(cuc) + ".   " + ( "%.3e" % (mmean[cuc] / ME) ) + " +/- " \
          +  ( "%.3e" % (mvar[cuc] / ME) ) + " Me\n"

    if cfg.plotesc:
      ax1.semilogy(self.t / tcol, self.m[:,0] / ME, 'r-')
      annstr += "esc. " + ( "%.3e" % (mmean[0] / ME) ) + " +/- " \
          +  ( "%.3e" % (mvar[0] / ME) ) + " Me\n"

    ax1.text(0.5, 0.1, annstr, transform = ax1.transAxes)

    ax1.set_ylim(cfg.mmin, cfg.mmax)
    ax1.set_xlim(cfg.tmin, cfg.tmax)
    ax1.set_title(simkey)
    ax1.set_xlabel("time [tcol]     tcol = " + str(tcol) + "s")
    ax1.set_ylabel("clump mass [Me]")
    
    ax2 = plt.subplot(212)
    for cuc in range(1,nopc+1):
      ax2.semilogy(self.t / tcol, self.dnopdtcol[:,cuc], '-',ms=1.0)
    
    if cfg.plotesc:
      ax2.semilogy(self.t / tcol, self.dnopdtcol[:,0], 'r-')
    
    ax2.set_xlim(cfg.tmin, cfg.tmax)
    ax2.set_ylim(0.9, 1.1*sim.nop)
    ax2.set_xlabel("time [tcol]     tcol = " + str(tcol) + "s")
    ax2.set_ylabel("no parts change per tcol [nop / tcol]")
    fig.savefig(pname)


