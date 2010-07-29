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
    self.tmin = -1.e4
    self.tmax =  1.e5

    self.ME = 5.97360e27

    self.noplotclumps = 8
    self.avgpts = 10
    self.plotesc = True


class SimClumps(object):
  def _countStepsClumps(self, step):
    self.nos += 1
    self.maxnoc = max(self.maxnoc, step.m.shape[0])
  
  def _getMasses(self, step):
    i = self.snames.index(step._v_pathname)
    cnoc = step.m.shape[0]
    self.m[i,0:cnoc] = step.m[:,0]
    self.t[i] = float( step._v_attrs["time"] )
    self.noc[i] = cnoc - 1

  def __init__(self,fname):
    dump = H5PartDump(fname)
    
    self.nos = 0
    self.maxnoc = 0
    dump.forEachStep(self._countStepsClumps)
    
    snames = dump.getStepNames()
    snames.sort(key=self.__getStepNum)
    self.snames = snames
    
    nos = self.nos
    maxnoc = self.maxnoc

    self.m   = np.zeros([nos,maxnoc])
    self.t   = np.zeros([nos,1])
    self.noc = np.zeros([nos])

    dump.forEachStep(self._getMasses)

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

  def systemBound(self):
    return True
  
  def __getStepNum(self,sname):
    return int(sname.replace("/Step#",""))

  def dmdtTau(self):
    return abs(self.m / self.dmdt) 

  #def dmdtAboveTau(self,tau, mmin):
  #  taum = abs( (self.m / self.dmdt) )
  #  return ( np.isnan(taum) | (taum > tau) | (self.m < mmin) )

  def _getMeanMasses(self,nop):
    if nop > self.nos:
       nop = self.nos
    mmean = np.mean( self.m[-nop-1:-1,:] , axis=0)
    mvar  = np.sqrt( np.var( self.m[-nop-1:-1,:] , axis=0) )
    return (mmean, mvar)

  def plot(self,pname,cfg,tit):
    fig = plt.figure()
    fig.clear()
    ax = plt.axes()
    
    nopc = min( cfg.noplotclumps, self.maxnoc - 1) 
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
    ax.set_xlim(cfg.tmin, cfg.tmax)
    ax.set_title(tit)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("clump mass [Me]")
    fig.savefig(pname)

  def plotdnop(self,pname,cfg,tit):
    fig = plt.figure()
    fig.clear()
    ax = plt.axes()
   
    #nopc = min( cfg.noplotclumps, self.maxnoc - 1) 
    nopc = min( 4, self.maxnoc - 1) 
    for cuc in range(1,nopc+1):
      ax.semilogy(self.t, self.dnopdtscl[:,cuc], '-',ms=1.0)
    
    #if cfg.plotesc:
    #  ax.semilogy(self.t, self.dnopdtscl[:,0], 'r-')

    ax.set_xlim(cfg.tmin, cfg.tmax)
    ax.set_ylim(0.9, 4.e5)
    ax.set_title(tit)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("no parts change per tscl [nop / tscl]")
    fig.savefig(pname)

