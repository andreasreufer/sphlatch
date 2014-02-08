#!/usr/bin/env python

from h5part import H5PartDump
import numpy as np
np.seterr(all='ignore')

from numpy import cross, dot, transpose, sqrt, arctan, arccos, arcsin, sin, cos, arctan2, power, pi, isnan, nan_to_num, mean, std

import matplotlib as mp
import matplotlib.pyplot as plt

from matplotlib.patches import Wedge

from plot_helpers import *
from const_cgs import *

mp.rc('text.latex', preamble = '\usepackage{amssymb}, \usepackage{wasysym}')
#mp.rcParams['font.serif'] = 'Minion Pro'
#mp.rcParams['font.family'] = 'serif'
mp.rcParams['font.size'] = 10.

rad2deg = 360/(2.*pi)

class nrange(object):
  def __init__(self, lower = None, upper = None):
    self.lower = lower
    self.upper = upper
  
  def intersection(self, aRange):
    if self.upper < aRange.lower or aRange.upper < self.lower:
      return None
    else:
      return nrange(max(self.lower,aRange.lower), \
                min(self.upper,aRange.upper))
  def getrange(self):
    return range(self.lower, self.upper)

def getDeriv1D(x, t):
  dt = t[1:] - t[0:-1]
  dxdth = ( x[1:] - x[0:-1] ) / dt
  dxdt  = np.zeros([x.shape[0]])
  dxdt[1:-1] = 0.5*(dxdth[1:] + dxdth[0:-1] )
  dxdt[0]  = dxdt[1]
  dxdt[-1] = dxdt[-2]
  return dxdt

def getDeriv2D(x, t):
  dt = t[1:] - t[0:-1]
  dxdth        = ( x[1:,:] - x[0:-1,:] ) / dt
  dxdt         = np.zeros([x.shape[0], x.shape[1]])
  dxdt[1:-1,:] = 0.5*(dxdth[1:,:] + dxdth[0:-1,:] )
  dxdt[0,:]    = dxdt[1,:]
  dxdt[-1,:]   = dxdt[-2,:]
  return dxdt

def getStableTimespan(t, cond, tmin):
  if len(cond) < 2 or not cond.any():
    return []
    
  # cond has at least 2 elements and at least 1 of them true
  ei = len(cond) - 1
    
  while ( ei > 1 ):
    while not cond[ei] and ei > 1:
      ei -= 1

    if not cond[ei]:
      return []

    si = ei
    while si > 0 and cond[si-1]:
      si -= 1

    if ( t[ei] - t[si] ) > tmin:
      return range(si, ei+1)
    else:
      ei = si - 1
  return []

def getStableMean(y, t, tau, tmin):
  dt = t[1:] - t[0:-1]
  dydth = ( y[1:] - y[0:-1] ) / dt
  dydt  = np.zeros( [y[:].shape[0]] )

  dydt[1:-1] = 0.5*( dydth[1:] - dydth[0:-1] )
  dydt[0] = 0.
  dydt[-1] = 0.

  cond = ( abs(y / dydt) > tau ) & ( t > 0 )

  xrng  = getStableTimespan( t, cond, tmin)
  ymean = mean(y[xrng])
  ydev  = std(y[xrng])

  return ( ymean, ydev, xrng )

def plotStableMean(y, t, tau, tmin, x, ax, str, col, logged=False):
  (ym, yv, yx) = getStableMean(y, t, tau, tmin)
  if len(yx) > 0:
    ystr = ( str % ym ) + r" +/- " + ( str % yv )
    if logged:
      ax.semilogy(x[yx], ym*np.ones(len(yx)), 'k--', lw=0.2)
    else:
      ax.plot(x[yx], ym*np.ones(len(yx)), 'k--', lw=0.2)
    ax.annotate(ystr, ( x[yx[0]], ym ), color=col )
  return (ym, yv, yx)
    


def clumpColor(i):
  if i == 0:
    return 'grey'
  if i == 1:
    return 'k'
  else:
    carr = ['r', 'g', 'b', 'c', 'm', 'y']
    return carr[(i-2) % 6]



class ClumpsPlotConfig(object):
    def __init__(self):
      self.mmin = 0.9e-5
      self.mmax = 2.5
      self.tmin =  -5.
      self.tmax =  55.
    
      self.ME = 5.97360e27
      self.ML = 7.34770e25

      self.noplotclumps =  8
      self.avgpts       = 10
      self.plotesc      = True


class MoonsPlotConfig(object):
    def __init__(self):
      self.mmin = 0.9e-5
      self.mmax = 2.5
      self.tmin =  -5.
      self.tmax = 105.
                          
      self.ME = 5.97360e27
      self.ML = 7.34770e25

      self.noplotclumps =  8
      self.avgpts       = 10
      self.plotesc      = True


class SimClumps(object):
  def __countStepsClumps(self, step):
    self.nos += 1
    self.maxnoc = max(self.maxnoc, step.m.shape[0])

    if self.loadsat:
      self.maxmat = step.mdiskmat.shape[1]
      self.satloaded = True
      
    if self.loadeng:
      self.maxmat = step.mmat.shape[1]
      self.engloaded = True
    
    if self.loadcrot:
      self.maxmat = step.mparentmat.shape[1]

  def __getBasicData(self, step):
    if step._v_pathname == '/Step#0' and not self.forceload:
      return

    i = self.snames.index(step._v_pathname)
    cnoc = step.m.shape[0]
    
    self.m[i,0:cnoc] = step.m[:,0]
    self.t[i] = float( step._v_attrs["time"] )
    self.noc[i] = cnoc - 1
    self.pos[i,0:cnoc,:] = step.pos[:,:]
    self.vel[i,0:cnoc,:] = step.vel[:,:]
    self.L[i,0:cnoc,:]   = step.L[:,:]
    self.P[i,0:cnoc,:]   = step.P[:,:]
    
    self.rho[i,0:cnoc] = step.rho[:,0]
    self.rc[i,0:cnoc]  = step.rc[:,0]

  def __getSatData(self, step):
    i = self.snames.index(step._v_pathname)
    cnoc = step.m.shape[0]
    cnom = step.mclmpmat[0,]
    
    self.rclmp[i,0:cnoc] = step.rclmp[:,0]
    self.posclmp[i,0:cnoc,:] = step.posclmp[:,:]
    
    self.mclmp[i,0:cnoc] = step.mclmp[:,0]
    self.mdisk[i,0:cnoc] = step.mdisk[:,0]
    self.mimp[i,0:cnoc]  = step.mimp[:,0]
    self.mesc[i,0:cnoc]  = step.mesc[:,0]
    
    self.Lclmp[i,0:cnoc,:] = step.Lclmp[:,:]
    self.Ldisk[i,0:cnoc,:] = step.Ldisk[:,:]
    self.Limp[i,0:cnoc,:]  = step.Limp[:,:]
    self.Lesc[i,0:cnoc,:]  = step.Lesc[:,:]

    self.mclmpmat[i,0:cnoc,:] = step.mclmpmat[:,:]
    self.mdiskmat[i,0:cnoc,:] = step.mdiskmat[:,:]
    self.mescmat[i,0:cnoc,:]  = step.mescmat[:,:]
    self.mimpmat[i,0:cnoc,:]  = step.mimpmat[:,:]

  def __getEnergies(self, step):
    i = self.snames.index(step._v_pathname)
    cnoc = step.m.shape[0]
    self.ETkin[i] = step._v_attrs["ekin"][0]
    self.ETpot[i] = step._v_attrs["epot"][0]
    self.ETthm[i] = step._v_attrs["ethm"][0]
    self.ETthmAdd[i] = step._v_attrs["ethermadded"][0]
      
    self.Ekin[i,0:cnoc]  = step.Ekin[:,0]
    self.Epot[i,0:cnoc]  = step.Epot[:,0]
    self.Erot[i,0:cnoc]  = step.Erot[:,0]
    self.Ethm[i,0:cnoc]  = step.Ethm[:,0]
      
    self.mmat[i,0:cnoc,:]  = step.mmat[:,:]

    self.I[i,0:cnoc]     = step.I[:,0]
    self.nop[i,0:cnoc]  = step.nop[:,0]
    if step._v_attrs._v_attrnames.count("noclumps") > 0:
      self.nofc[i] = step._v_attrs["noclumps"][0]
      
  def __getClmpeng(self, step):
    i = self.snames.index(step._v_pathname)
    cnoc = step.m.shape[0]
    
    self.Erotclmp[i,0:cnoc]  = step.Erotclmp[:,0]
    self.Epotclmp[i,0:cnoc]  = step.Epotclmp[:,0]
    self.Iclmp[i,0:cnoc]     = step.Iclmp[:,0]
  
  def __getClmprot(self, step):
    i = self.snames.index(step._v_pathname)
    cnoc = step.m.shape[0]
    
    self.Lparent[i,0:cnoc,:] = step.Lparent[:,:]
    self.Iparent[i,0:cnoc]  = step.Iparent[:,0]
    self.mparentmat[i,0:cnoc,:]  = step.mparentmat[:,:]

  def __getStepNum(self,sname):
    return int(sname.replace("/Step#",""))
  
  def __getdiskData(self, step):
    i = self.snames.index(step._v_pathname)

  def __init__(self,sim,fname="clumps.h5part",loadsat=False,loadeng=False,loadclmp=False,loadcrot=False,forceload=False):
    dump = H5PartDump(sim.dir + fname)
    self.loadsat = loadsat
    self.loadeng = loadeng
    self.loadclmp = loadclmp
    self.loadcrot = loadcrot
    self.forceload = forceload
    self.sim = sim
      
    self.satloaded = False
    self.engloaded = False
    
    self.nos = 0
    self.maxnoc = 0
    dump.forEachStep(self.__countStepsClumps)
    
    snames = dump.getStepNames()
    snames.sort(key=self.__getStepNum)
    self.snames = snames
    
    nos = self.nos
    maxnoc = self.maxnoc

    self.m   = np.zeros([nos,maxnoc])
    self.t   = np.zeros([nos])
    self.noc = np.zeros([nos])
    
    self.pos = np.zeros([nos,maxnoc,3])
    self.vel = np.zeros([nos,maxnoc,3])
    self.L   = np.zeros([nos,maxnoc,3])
    
    self.P   = np.zeros([nos,maxnoc,3])
    self.rho = np.zeros([nos,maxnoc])
    self.rc  = np.zeros([nos,maxnoc])

    dump.forEachStep(self.__getBasicData)

    self.mLR = np.zeros([nos])
    self.mSR = np.zeros([nos])
    self.mXR = np.zeros([nos])
      
    if maxnoc > 1:
      self.mLR = self.m[:,1]

    if maxnoc > 2:
      self.mSR = self.m[:,2]
      self.mXR = self.m[:,2:].sum(axis=1)
    
    if loadsat:
      maxmat = self.maxmat
      self.posclmp = np.zeros([nos,maxnoc,3])
      self.rclmp = np.zeros([nos,maxnoc])
      self.mclmp = np.zeros([nos,maxnoc])
      self.mdisk = np.zeros([nos,maxnoc])
      self.mimp  = np.zeros([nos,maxnoc])
      self.mesc  = np.zeros([nos,maxnoc])
      
      self.Lclmp = np.zeros([nos,maxnoc,3])
      self.Ldisk = np.zeros([nos,maxnoc,3])
      self.Limp  = np.zeros([nos,maxnoc,3])
      self.Lesc  = np.zeros([nos,maxnoc,3])
      
      self.mclmpmat = np.zeros([nos,maxnoc,maxmat])
      self.mdiskmat = np.zeros([nos,maxnoc,maxmat])
      self.mescmat =  np.zeros([nos,maxnoc,maxmat])
      self.mimpmat =  np.zeros([nos,maxnoc,maxmat])
      
      dump.forEachStep(self.__getSatData)

      self.Lbound  = self.Lclmp + self.Ldisk + self.Limp

    self.rmaxbound   = dump.getAttr( dump.getStepNames()[0], "rmaxbound")
    self.rmaxunbound = dump.getAttr( dump.getStepNames()[0], "rmaxunbound")
    
    if loadeng:
      maxmat = self.maxmat
      self.ETkin = np.zeros([nos])
      self.ETpot = np.zeros([nos])
      self.ETthm = np.zeros([nos])
    
      self.ETthmAdd = np.zeros([nos])
      
      self.Ekin  = np.zeros([nos,maxnoc])
      self.Epot  = np.zeros([nos,maxnoc])
      self.Erot  = np.zeros([nos,maxnoc])
      self.Ethm  = np.zeros([nos,maxnoc])

      self.I     = np.zeros([nos,maxnoc])
      self.nop   = np.zeros([nos,maxnoc])
    
      self.mmat =  np.zeros([nos,maxnoc,maxmat])
      self.mmat[:,0,:] = self.sim.mtot - self.mmat[:,1:,:].sum(axis=1)
      
      self.nofc = np.zeros([nos])

      dump.forEachStep(self.__getEnergies)
      
      self.ETtot = self.ETkin + self.ETpot + self.ETthm 
      #self.ETtot = self.ETkin + self.ETpot + self.ETthm + self.Erot
      self.Tot  = self.Erot / (-self.Epot)
    
    if loadclmp:
      self.Epotclmp  = np.zeros([nos,maxnoc])
      self.Erotclmp  = np.zeros([nos,maxnoc])
      self.Iclmp     = np.zeros([nos,maxnoc])

      dump.forEachStep(self.__getClmpeng)
    
    if loadcrot:
      maxmat = self.maxmat
      self.Lparent = np.zeros([nos,maxnoc,3])
      self.Iparent = np.zeros([nos,maxnoc])
      self.mparentmat = np.zeros([nos,maxnoc,maxmat])
      
      dump.forEachStep(self.__getClmprot)

    self.m[:,0] = self.sim.mtot - self.m[:,1:].sum(axis=1)
    m = self.m
    
    t = self.t
    noc = self.noc
    
    self.xi = ( self.mLR - self.sim.tarb.m ) / self.sim.impb.m

    dt = t[1:] - t[0:-1]
    self.dt = dt

  def __del__(self):
    pass

  def _getNoClumps(self,mmin=0.):
    nocraw = ( self.m[:,1:] > mmin ).sum(axis=1)
    tiok    = np.logical_not( ( ( self.m[:,1:] > mmin ) & \
        ( self.r[:,1:] > self.rmaxbound ) ).any(axis=1) )
    return ( nocraw.compress(tiok), tiok )

  def storeResults(self):
    res = self.sim.results
    (noc, tiok) = self._getNoClumps()
    res.noclumps = noc[-1]

  def getStableMassTimespan(self, tau, tmin, cids=(1,2)):
    m   = self.m
    t   = self.t
    nos = self.nos
    maxnoc = self.maxnoc

    mnrng = nrange(0, nos)
    
    # try to find a common range
    for cid in cids:
      # more clumps requested to analyse, than exist
      if cid >= maxnoc:
        return range(nos-1, nos)

      (mim, miv, rngi) = getStableMean( m[:,cid], t, tau, tmin)
      if len(rngi) > 0:
        mnrng = mnrng.intersection(nrange( rngi[0], rngi[-1]))
      else:
        mnrng = nrange(0,0)

    mrng = mnrng.getrange()

    # if common range is empty, fall back to last data point
    if len(mrng) == 0:
      mrng = range(nos-1, nos)

    return mrng

  def getTimespan(self, tmin, tmax):
    t   = self.t
    nos = self.nos

    tok = ( ( ( t >= tmin ) & ( t <= tmax ) ).nonzero() )[0]
    if tok.shape[0] > 0:
      rng = range( np.min( tok ), np.max( tok ) )
    else:
      rng = range(nos-1, nos)
    return rng



  def plotLRmoons(self,pname,cfg):
    sim = self.sim
    fig = plt.figure()
    fig.clear()

    simkey = sim.params.key
    tcol   = sim.tcol
    
    t   = self.t
    m   = self.m
    
    #self.dmdt = getDeriv2D(m, t)

    nopc = min( cfg.noplotclumps, self.maxnoc - 1) 
    ME   = cfg.ME
    ML   = cfg.ML
    
    annstr = ""
    
    nos = self.nos
    dt  = self.dt
    x   = ( t / 3600 )

    
    # plot composition of the disk
    msi02 = self.mdiskmat[:,1,1]
    mwatr = self.mdiskmat[:,1,2]
    miron = self.mdiskmat[:,1,5]

    ax1 = plt.subplot(221)
    ax1.semilogy( x, (msi02+1.) / ML, '-', label='SiO2', c="red", nonposy='clip') 
    ax1.semilogy( x, (mwatr+1.) / ML, '-', label='H2O', c="darkcyan") 
    ax1.semilogy( x, (miron+1.) / ML, '-', label='iron', c="mediumblue" ) 
    
    tau  = 10.*tcol
    tmin =  3.*tcol

    self.dmsiO2 = plotStableMean(msi02 / ML, t, ax1, x, "$%6.4f \mathrm{M_L}$", "red", tau, tmin, logged=True)
    
    self.dmiron = plotStableMean(miron / ML, t, ax1, x, "$%6.4f \mathrm{M_L}$", "mediumblue", tau, tmin, logged=True)

    self.dmwatr = plotStableMean(mwatr / ML, t, ax1, x, "$%6.4f \mathrm{M_L}$", "darkcyan", tau, tmin, logged=True)


    paramstxt = '$m_T =' + ("%3.2f" % sim.params.mtar) + r' M_E, ' +\
        'm_I = ' + ("%3.2f" % sim.params.mimp) + 'M_E, ' +\
        r'v_{imp} = ' + ("%3.2f" % sim.params.vimprel) + r'v_{esc}, ' +\
        ("%2.0f" % sim.params.impa) + r'^\circ$'

    ax1.set_ylim( 1.e-2, 10.0 )
    #ax1.set_xlim( cfg.tmin, cfg.tmax)
    ax1.set_title(paramstxt)

    # plot 
    ax2 = plt.subplot(222)
    ax2.plot(x, self.rclmp[:,1], 'r-', ms=1.0)
    ax2.plot(x, self.rc[:,1],    'b-', ms=1.0)
    #ax2.set_xlim( cfg.tmin, cfg.tmax)
    #ax2.set_ylim( 5.5e8, 8.0e8 )


    # plot mass categories
    mclmp = self.mclmp[:,1]
    mimp  = self.mimp[:,1]
    mdisk = self.mdisk[:,1]
    mesc  = self.mesc[:,1]

    ax3 = plt.subplot(223)
    ax3.semilogy( x, mclmp / ME, 'r-' ) 
    ax3.semilogy( x, mimp  / ME, 'g-' ) 
    ax3.semilogy( x, mdisk / ME, 'b-' ) 
    ax3.semilogy( x, mesc  / ME, 'k:' ) 
    #ax3.set_xlim( cfg.tmin, cfg.tmax)
    ax3.set_ylim( 1.e-5, 1.0 )
    
    Lclmp = self.Lclmp[:,1,2]
    Limp  = self.Limp[:,1,2]
    Ldisk = self.Ldisk[:,1,2]
    Lesc  = self.Lesc[:,1,2]

    LboundLRz = self.Lbound[:,1,2]
    self.LboundLRz = LboundLRz

    # plot distribution of angular momentum
    ax4y1 = ( Lclmp                       ) / lem
    ax4y2 = ( Lclmp + Limp                ) / lem
    ax4y3 = ( Lclmp + Limp + Ldisk        ) / lem
    ax4y4 = ( Lclmp + Limp + Ldisk + Lesc ) / lem

    ax4 = plt.subplot(224)
    ax4.plot( x, ax4y4, 'k:' )
    ax4.fill_between(x, 0, ax4y1, facecolor='r', lw=0.)
    ax4.fill_between(x, ax4y1, ax4y2, facecolor='g', lw=0.)
    ax4.fill_between(x, ax4y2, ax4y3, facecolor='b', lw=0.)
    
    #self.dLboundLRzdt = plotStableMean(LboundLRz / lem, t, ax4, x, "$%6.3f \mathrm{lem}$", "black", tau, tmin)
    
    #ax4.set_xlim( cfg.tmin, cfg.tmax)
    ax4.set_ylim( 0., 3. )
    
    fig.savefig(pname)

  def getResults(self,tmin=0.90,tmax=0.99):
    sim = self.sim
    res = sim.results
    maxnoc = self.maxnoc
    t = self.t

    nomat  = self.mmat.shape[2]

    if maxnoc <= 1:
      res.valnop  = -1
      return
    
    #tau = 10.*sim.tcol
    #tmin= 5.*sim.tcol
    #rng = self.getStableMassTimespan(tau, tmin)
    
    rng = self.getTimespan(tmin*sim.tstop, tmax*sim.tstop)
    print "results range: ", np.min(t[rng])/sim.tcol, "-", np.max(t[rng])/sim.tcol, " tcol"
    
    self.rng = rng
    if ( self.nofc[rng[:]] < 0 ).any():
      print "missing clump search in results range, invalidating results ..."
      res.problem = True
      return

    self.getDynamics2D()

    res.resize(maxnoc,nomat)
    res.mm = np.mean(self.m[rng,:], axis=0)
    res.mv = np.std(self.m[rng,:], axis=0)
    
    res.mmatm = np.mean(self.mmat[rng,:,:], axis=0)
    res.mmatv = np.std(self.mmat[rng,:,:], axis=0)

    res.U0 = self.Ethm[0,:]
    res.Um = np.mean(self.Ethm[rng,:], axis=0)
    res.Uv = np.std(self.Ethm[rng,:], axis=0)
    res.Ua = np.mean(self.ETthmAdd[rng])
    
    res.Epot0 = self.Epot[0,:]
    res.Epotm = np.mean(self.Epot[rng,:], axis=0)
    res.Epotv = np.std(self.Epot[rng,:], axis=0)
    
    res.Ekin0 = self.Ekin[0,:]
    res.Ekinm = np.mean(self.Ekin[rng,:], axis=0)
    res.Ekinv = np.std(self.Ekin[rng,:], axis=0)
    
    res.Erot0 = self.Erot[0,:]
    res.Erotm = np.mean(self.Erot[rng,:], axis=0)
    res.Erotv = np.std(self.Erot[rng,:], axis=0)

    res.Etot0 = self.ETtot[0]
    res.dEtotm = np.mean( ( self.ETtot[:] - self.ETtot[0] ) / self.ETtot[0] )
    res.dEtotv = np.var( ( self.ETtot[:] - self.ETtot[0] ) / self.ETtot[0] )

    res.vinfm = np.mean(self.vinf[rng,:,:], axis=0)
    res.vinfv = np.std(self.vinf[rng,:,:], axis=0)

    res.dblvarthetam = np.mean(self.dblvartheta[rng,:], axis=0)
    res.dblvarthetav =  np.std(self.dblvartheta[rng,:], axis=0)
    
    res.rhom = np.mean(self.rho[rng,:], axis=0)
    res.rhov = np.std(self.rho[rng,:], axis=0)

    for i in range(0,nomat):
      res.compm[:,i] = np.mean(self.mmat[rng,:,i] / self.m[rng,:], axis=0)
      res.compv[:,i] = np.std(self.mmat[rng,:,i] / self.m[rng,:], axis=0)
    
    res.Lm = np.mean(self.L[rng,:,:], axis=0)
    res.Lv = np.std(self.L[rng,:,:], axis=0)
    
    res.Im = np.mean(self.I[rng,:], axis=0)
    res.Iv = np.std(self.I[rng,:], axis=0)
    
    res.nopm = np.mean(self.nop[rng,:], axis=0)
    res.nopv = np.std(self.nop[rng,:], axis=0)

    res.valid = True
    res.valtmin = self.t[ rng[0] ]
    res.valtmax = self.t[ rng[-1] ]
    res.valnop  = len(rng)

  # works only in the XY-plane
  def getDynamics2D(self):
    sim = self.sim
    nos = self.nos
    t   = self.t
    dt  = self.dt
    
    tau = 10.*sim.tcol
    tmin= 5.*sim.tcol

    mLR = self.mLR
    
    th = self.t / 3600.
    
    nos  = self.m.shape[0]
    norc = self.m.shape[1]
    
    self.r = np.zeros([nos, norc])
    
    # introduce new variables
    self.dblvartheta = np.zeros([nos, norc])
    self.vrelinf     = np.zeros([nos, norc, 3])
    self.Pinf        = np.zeros([nos, norc, 3])
    self.vinf        = np.zeros([nos, norc, 3])
   

    for i in range(2,norc):
      mi = self.m[:,i]
      mu = sim.G*(mi + mLR)
      rpos  = self.pos[:,1,:] - self.pos[:,i,:]
      rvel  = self.vel[:,1,:] - self.vel[:,i,:]

      r     = np.sqrt( (rpos*rpos).sum(axis=1) )
      v     = np.sqrt( (rvel*rvel).sum(axis=1) )

      # v relative to LR
      vrelinf  = np.sqrt( v*v - 2*sim.G*mLR / r )
      self.r[:,i] = r

      rv    = ( rpos * rvel ).sum(axis=1)
      beta  = arcsin( rv / (r*v) )

      alpha = arctan2( rvel[:,1], rvel[:,0] )

      k1 = r*v*v / mu

      theta = arctan2( (k1*sin(beta)*cos(beta)),(k1*cos(beta)*cos(beta) - 1.) )
      e = sqrt( (k1-1.)*(k1-1)*cos(beta)*cos(beta)  + sin(beta)*sin(beta) )
      
      vartheta = arctan( 1./sqrt( e*e - 1.) )
      
      alphaneginf = sim.gi.alphaneginf
      alphaposinf    = alpha - theta + beta + vartheta

      dblvartheta = alphaposinf - alphaneginf
      
      self.vrelinf[:,i,0] = nan_to_num( -vrelinf*cos( alphaposinf ) )
      self.vrelinf[:,i,1] = nan_to_num( -vrelinf*sin( alphaposinf ) )

      self.dblvartheta[:,i] = dblvartheta

    # assuming a center of mass system, we get vinf for the largest remnant
    # from momentum conservation, get vinf for
    Prelinf = np.zeros([nos, norc, 3])
    Prelinf[:,:,0] = self.vrelinf[:,:,0]*self.m
    Prelinf[:,:,1] = self.vrelinf[:,:,1]*self.m

    mclmps = self.m[:,1:].sum(axis=1)
    self.vinf[:,1,0] =  - Prelinf.sum(axis=1)[:,0] / mclmps
    self.vinf[:,1,1] =  - Prelinf.sum(axis=1)[:,1] / mclmps
    self.Pinf[:,1,0] = self.vinf[:,1,0] * self.m[:,1]
    self.Pinf[:,1,1] = self.vinf[:,1,1] * self.m[:,1]
    
    for i in range(2,norc):
      self.vinf[:,i,0] = self.vrelinf[:,i,0] + self.vinf[:,1,0]
      self.vinf[:,i,1] = self.vrelinf[:,i,1] + self.vinf[:,1,1]
      self.Pinf[:,i,0] = self.vinf[:,i,0] * self.m[:,i]
      self.Pinf[:,i,1] = self.vinf[:,i,1] * self.m[:,i]



  def plotSummary(self,pname):
    sim = self.sim
    res = self.sim.results
    
    self.getDynamics2D()
    
    
    fig = plt.figure( figsize=(11.69, 8.27), dpi=50 )
    fig.clear()

    t  = self.t
    th = self.t / 3600
    tc = self.t / sim.tcol
    
    trhmin = res.valtmin / 3600.
    trhmax = res.valtmax / 3600.

    tcol = sim.tcol

    res  = sim.results
    rval = res.valid
    rrng = (res.valtmin / 3600., res.valtmax / 3600.)
    
    leng = self.loadeng

    mtar    = self.sim.tarb.m
    mtarmat = self.sim.tarb.mmat
    mimp    = self.sim.impb.m
    mimpmat = self.sim.impb.mmat

    vesc = sim.gi.vesc
    
    nos   = self.m.shape[0]
    norc  = self.m.shape[1]
    nomat = self.mmat.shape[1]
    
    
    # composition color map 
    cmap = np.zeros( [32,3] )
    cmap[0] = cnameToRGB("plum")
    cmap[1] = cnameToRGB("red") 
    cmap[2] = cnameToRGB("darkcyan") 
    cmap[4] = cnameToRGB("red") 
    cmap[5] = cnameToRGB("mediumblue") 

    # FIXME: get individual energie & composition

    # accretion efficiency
    ax1 = plt.subplot(3,4,1)
    
    ax1.set_title("accretion efficiency")
    ax1.plot(th, self.xi, 'g', lw=1.0)
    
    if rval:
      xires = (res.mm[1] - mtar)/mimp
      ax1.hlines(xires, rrng[0], rrng[1], linestyles='-', lw=0.25)

    if leng:
      for i in (1,2,5):
        xii = ( self.mmat[:,1,i] - mtarmat[i] ) / mimpmat[i]
        ax1.plot(th, xii, '-', color=cmap[i], lw=0.5)
    

    # accretion efficiency
    ax2 = plt.subplot(3,4,2)
    
    ax2.set_title("striping efficiency")
    ax2.plot(th, (mimp - self.m[:,2])/mimp, 'g', lw=1.0)
    
    if leng:
      for i in (1,2,5):
        chii = ( mimpmat[i] - self.mmat[:,2,i]) / mimpmat[i]
        ax2.plot(th, chii, '-', color=cmap[i] , lw=0.5)
      
    
    
    Ered = sim.gi.Ered
    ax3 = plt.subplot(3,4,3)
    ax3.set_title(r" $E_{pot}, \Delta U / \frac{1}{2} \mu v_{imp}^2$")

    dUT = self.ETthm -self.ETthm[0] 
    dETpot = self.ETpot-self.ETpot[0]

    ax3.plot( th,    dUT / Ered, 'r', lw=0.5)
    ax3.plot( th, dETpot / Ered, 'g', lw=0.5)
    ax3.plot( th, self.ETthmAdd / Ered, 'r--', lw=0.5)

    ax4 = plt.subplot(3,4,4)
    ax4.set_title(r"rot. stability $- E_{rot} / E_{pot}$")
    ax5 = plt.subplot(3,4,5)
    ax5.set_title(r"$M / M_{imp}$")
    ax6 = plt.subplot(3,4,6)
    ax6.set_title(r"$v_{rel, \inf} / v_{esc}$")
    ax7 = plt.subplot(3,4,7)
    ax7.set_title(r"deflection angle $\vartheta$")
    
    vrelinfsc = np.sqrt( (self.vrelinf * self.vrelinf).sum(axis=2) )
    vrelminf = np.sqrt( sim.params.vimprel*sim.params.vimprel - 1. )
    
    ax6.hlines(vrelminf, th[0], th[-1], linestyles='--', lw=0.25)
    ax7.hlines(2.*sim.gi.vartheta*rad2deg, th[0], th[-1], linestyles='--', lw=0.25)
    # LR fraction
    ax8 = plt.subplot(3,4,8, frameon=False)
    ax8.set_title(r"resulting fragments")
    
    # relative energy per mass
    ax9 = plt.subplot(3,4,9)
    ax9.set_title(r"rel. internal energy $U_i / M_i$")

    ax10 = plt.subplot(3,4,10)
    ax10.set_title(r"ang. momentum $L_z$")
    
    ax11 = plt.subplot(3,4,11)
    ax11.set_title(r"cumulative mass $M_i$")
    
    ax12 = plt.subplot(3,4,12)
    ax12.set_title(r"$\Delta E_{tot} / E_{tot 0}$")
    ax12.plot( th, ( self.ETtot[:] - self.ETtot[0] ) / self.ETtot[0], '-k', lw=0.5 )
    ax12.hlines(0., th[0], th[-1], linestyles='--', lw=0.25)
    
    refp = []
    cmpp = []
    
    tht1 = 0.; tht2 = 0.
    thi1 = 0.; thi2 = 0.

    # the target body composition
    if rval:
      for j in range(0,32):
        tht2 = tht1 + 360.*(mtarmat[j] / mtar)
        thi2 = thi1 + 360.*(mimpmat[j] / mimp)
      
        rt = 1.0
        ri = sqrt(mimp / mtar)
        rm = sqrt( (mtar+mimp) / mtar )
        rd = 0.05

        refp.append( Wedge( (-rm-rd, rm+rd), rt, tht1, tht2, fill=False, color='k', lw=0.1) )
      
        refp.append( Wedge( (rm+rd,rm+rd), ri, thi1, thi2, fill=False, color='k', lw=0.1) )

        r0 = sqrt( res.mmatm[0,j] / mtarmat[j] )
        r1 = sqrt( res.mmatm[1,j] / mtarmat[j] )

        cmpp.append( Wedge( (-rm-rd, rm+rd), r1, tht1, tht2, fill=True, color=cmap[j], lw=0.) )
        ax8.text( -rm-rd, rm+rd-0.5, r"$ M_{LR} $")

        cmpp.append( Wedge( (rm+rd,rm+rd), r0, thi1, thi2, fill=True, color=cmap[j], lw=0.) )
        ax8.text( rm+rd, rm+rd-0.6*ri, r"$ M_{esc} $")
        
        # secondary remnants
        for k in range(2, min(7,norc) ):
          pltpos = ( (k-2.)*(ri+rd)*1.5 - 2*rm + ri, -(rm+rd) )
          if mtarmat[j] > 0.:
            rk = sqrt( res.mmatm[k,j] / mtarmat[j] )
            cmpp.append( Wedge( pltpos, rk, thi1, thi2, fill=True, color=cmap[j], lw=0.) )
            refp.append( Wedge( pltpos, ri, thi1, thi2, fill=False, color='k', lw=0.1) )
            ax8.text(pltpos[0], pltpos[1]-0.6*ri, 
                r"$ M_{" + str(k) + r"} $")


        thi1 = thi2
        tht1 = tht2

      for patch in cmpp:
        ax8.add_patch(patch)
      for patch in refp:
        ax8.add_patch(patch)

      ax8.set_xlim( -2*(rm+rd), 2*(rm+rd) )
      ax8.set_ylim( -2*(rm+rd), 2*(rm+rd) )

    for i in range(0,norc):
      cc = clumpColor(i)
      
      ax5.semilogy(th, (self.m[:,i]+1.) / mimp, lw=0.5, color=cc)

      if (rval and res.mm[i] > 1.e-2*sim.mtot) or not rval:
        T = -self.Erot[:,i] / self.Epot[:,i]
        Tval = (T > 0.)
        ax4.semilogy(th.compress(Tval), T.compress(Tval), lw=0.5,color=cc)
        
        ax6.plot(th, vrelinfsc[:,i] / vesc, lw=0.5,color=cc)

        ax7.plot(th, self.dblvartheta[:,i] * rad2deg, lw=0.5, color=cc)

    
        # stuff which is not plotted for escaping mass
        if not i == 0:
          Urel= self.Ethm[:,i] / self.m[:,i]
          Urelval = (self.m[:,i] > 0)
          ax9.plot(th.compress(Urelval), Urel.compress(Urelval), lw=0.5,color=cc) 

          Lzi    = self.L[:,i,2]
          Lzipos = Lzi > 0.
          Lzineg = Lzi < 0.
          ax10.semilogy(th.compress(Lzipos), Lzi.compress(Lzipos), lw=0.5,color=cc)
          ax10.semilogy(th.compress(Lzineg), -Lzi.compress(Lzineg), lw=0.5,color=cc, ls='--')



    mprev  = np.zeros([nos])
    mnext  = np.zeros([nos])
    
    for i in range(norc-1,-1,-1):
      cc = clumpColor(i)
      mnext = mprev + self.m[:,i]
      ax11.fill_between(th, mprev, mnext, facecolor=cc, lw=0.)
      mprev = mnext

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7):
      ymin = ax.axis()[2]
      ymax = ax.axis()[3]
      ax.vlines(trhmin, ymin, ymax, linestyles='-', lw=0.25)
      ax.vlines(trhmax, ymin, ymax, linestyles='-', lw=0.25)
    
    ax1.set_ylim(-1.0, 1.1)
    ax2.set_ylim(-0.1, 1.1)
    ax4.set_ylim( 1.e-3, 5.)
    ax5.set_ylim(1.e-3, 1.1*(sim.mtot / sim.impb.m) )
    ax6.set_ylim(-0.1*sim.params.vimprel, 1.2*sim.params.vimprel)
    
    paramstxt = '$m_T =' + ("%3.2f" % sim.params.mtar) + r' M_E, ' +\
        'm_I = ' + ("%4.3f" % sim.params.mimp) + 'M_E, ' +\
        r'v_{imp} = ' + ("%3.2f" % sim.params.vimprel) + r'v_{esc}, ' +\
        ("%2.0f" % sim.params.impa) + r'^\circ$'
    
    norma = plt.axes( [0.0, 0.0, 1.0, 1.0], frameon=False)
    norma.text( 0.03, 0.95, paramstxt)

    fig.savefig(pname)
        
  def plotDynamics(self,pname,cfg):
    sim = self.sim
    fig = plt.figure()
    fig.clear()

    nopc = 2
    ME   = cfg.ME
    ML   = cfg.ML
    
    annstr = ""
    
    nos = self.nos
    t   = self.t
    dt  = self.dt
    x   = ( t / 3600 )[:,0]
    
    # plot composition of the disk
    ax1 = plt.subplot(221)
    
    rpos12 = self.pos[:,1,:] - self.pos[:,2,:]
    r12    = np.sqrt( ( rpos12 * rpos12 ) .sum(axis=1) )
    
    mUbd   = self.m[:,0]
    mLR    = self.m[:,1]
    mSR    = self.m[:,2]
    
    mtar   = sim.tarb.m
    mimp   = sim.impb.m

    accreff = (mLR - mtar)/mimp

    v2SR = np.sum( self.vel[:,2,:] * self.vel[:,2,:],axis=1)
    vesc   = sim.gi.vesc

    G = sim.G
    k1 = 1. + (mSR / mLR)

    Etot = (-G*mLR*mSR / r12) + 0.5*mSR*k1*v2SR
    vSRinf = np.sqrt( 2.*Etot / (mSR * k1 ) )
    vinf   = vSRinf*k1
    
    tau = 10.*sim.tcol
    tmin= 5.*sim.tcol
    
    ( accreffm, accreffv, accreffx) = getStableMean(accreff, t[:,0], tau, tmin)
    sim.results.accreffm = accreffm
    sim.results.accreffv = accreffv
    print "accreff: ",accreffm, accreffv, len(accreffx)

    self.accreffm = accreffm
    self.accreffv = accreffv
    self.accreffx = accreffx

    ( vinfm, vinfv, vinfx ) = getStableMean(vinf, t[:,0], tau, tmin)
    sim.results.vinfclmpm.append( vinfm )
    sim.results.vinfclmpv.append( vinfv )
    print "vinf: ",vinfm, vinfv, len(vinfx)

    self.vinfm = vinfm
    self.vinfv = vinfv
    self.vinfx = vinfx


    for cuc in range(1,nopc+1):
      P = np.sqrt( ( self.P[:,cuc,:] * self.P[:,cuc,:] ).sum(axis=1) )
      P0 = P[0]
      ax1.semilogx( r12, P / P[0], ms=1.0)
    ax1.set_title("P [P0] vs. r [cm]")
    ax1.semilogx(r12[0], 1, 'kx', ms=10.)
    
    # plot mass categories
    ax2 = plt.subplot(223)
    for cuc in range(1,nopc+1):
      Ekin = 0.5*( self.vel[:,cuc,:] * self.vel[:,cuc,:] ).sum(axis=1)
      ax2.semilogx(r12, Ekin / Ekin[0], ms=1.0)
    ax2.semilogx(r12[0], 1, 'kx', ms=10.)
    ax2.set_title("Ekin/M [E0/M] vs. r [cm]")
    

    # plot 
    ax3 = plt.subplot(222)
    ax3.set_title("v_inf / vesc vs. t [h]")
    ax3.plot(x[0], sim.gi.vinf / vesc, 'kx',ms=10.)
    ax3.plot(x, vinf / vesc, '-', ms=1.0)
    ax3.plot(x[vinfx], vinfm/vesc*np.ones([len(vinfx)]), 'k--',ms=1.0,lw=0.5)
    ax3.set_xlim( x[0], x[-1])
    ax3.set_ylim( 0., 2.0 )


    ax4 = plt.subplot(224)
    for cuc in range(1,nopc+1):
      M0 = self.m[0,cuc]
      ax4.plot(x, self.m[:,cuc] / M0, '-',ms=1.0)

    ax4.plot(x, accreff, 'k--',ms=1.0)
    ax4.plot(x[accreffx], accreffm*np.ones([len(accreffx)]), 'k-',ms=1.0,lw=0.5)
    ax4.set_title("M [M0] vs. t [h], accr. eff.")
    ax4.set_xlim( x[0], x[-1])
    ax4.set_ylim( -0.5, 2.0)

    paramstxt = '$m_T =' + ("%3.2f" % sim.params.mtar) + r' M_E, ' +\
        'm_I = ' + ("%3.2f" % sim.params.mimp) + 'M_E, ' +\
        r'v_{imp} = ' + ("%3.2f" % sim.params.vimprel) + r'v_{esc}, ' +\
        ("%2.0f" % sim.params.impa) + r'^\circ$'
    
    norma = plt.axes( [0.0, 0.0, 1.0, 1.0], frameon=False)
    norma.text( 0.03, 0.95, paramstxt)

    self.vinf = vinf
    
    fig.savefig(pname)

