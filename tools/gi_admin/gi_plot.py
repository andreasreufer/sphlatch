#!/usr/bin/env python
import sys
import numpy as np

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.mpl    as mpl

from fewbody import FewBodies
from h5part import H5PartDump

#mp.rc('text', usetex=True)
mp.rc('text.latex', preamble = '\usepackage{amssymb}, \usepackage{wasysym}')

eVinK = 11604.505

def plotDummy(norma, physa, plt, pdump, cdump, cfg):
  pass

def plotGrid(norma, physa, plt, pdump, cdump, cfg):
  physa.grid(True, lw=0.1, color='darkgrey')
    
def plotSimParams(norma, physa, plt, pdump, cdump, cfg):
  paramstxt = '$m_T =' + ("%3.2f" % cfg.mtar) + r' M_E, ' +\
      'm_I = ' + ("%3.2f" % cfg.mimp) + 'M_E, ' +\
      r'v_{imp} = ' + ("%3.2f" % cfg.vimp) + r'v_{esc}, ' +\
      ("%2.0f" % cfg.impa) + r'^\circ, '+\
      'T = ' + ("%4.0f" % cfg.T ) + 'K$'

  norma.text( cfg.parm_vc[0], cfg.parm_vc[1], paramstxt, \
      color=cfg.parm_fc, size=cfg.parm_txts )

def cnameToRGB(str):
  return np.array( col.hex2color( col.cnames[str] ) )

def cnameToRGBA(str, alpha=1.0):
  rgba = np.zeros(4)
  rgba[0:3] = np.array( col.hex2color( col.cnames[str] ) )
  rgba[3] = alpha
  return rgba

def colorBodyAndMat(pdump, cdump, filt, cfg):
  cmap = np.zeros((2,6,3))
  cmap[0,2] = cnameToRGB("blue")
  cmap[0,4] = cnameToRGB("red")
  cmap[0,5] = cnameToRGB("darkgrey")
  cmap[1,2] = cnameToRGB("aqua")
  cmap[1,4] = cnameToRGB("orange")
  cmap[1,5] = cnameToRGB("grey")

  mat = pdump.mat[:,0].compress(filt)[:]
  bod = np.int32( (pdump.id[:].compress(filt)[:] > cfg.idfilt) )

  color = cmap[bod, mat]
  pt2size = cfg.pt2size

  return (color, pt2size)

def colorTemp(pdump, cdump, filt, cfg):
  T = (pdump.T[:,0] * eVinK)
  Tmin = cfg.Tmin
  Tmax = cfg.Tmax

  color = (T - Tmin) / (Tmax - Tmin)
  pt2size = cfg.pt2size
  return (color, pt2size)

def colorDens(pdump, cdump, filt, cfg):
  rho = (pdump.rho[:,0])
  rhomin = cfg.rhomin
  rhomax = cfg.rhomax

  color = (rho - rhomin) / (rhomax - rhomin)
  pt2size = cfg.pt2size
  return (color, pt2size)

def colorClumps(pdump, cdump, filt, cfg):
  cid = pdump.clumpid[:,0].compress(filt)[:]

  maxnoc = 10
  maxcid = maxnoc - 1
  cid[cid > maxcid] = maxcid
  
  cmap = np.zeros((maxnoc,4))

  cmap[0] = cnameToRGBA("lime") # escaping stuff
  cmap[1] = cnameToRGBA("red") #
  cmap[2] = cnameToRGBA("blue") # 
  cmap[3] = cnameToRGBA("yellow") # 
  cmap[4] = cnameToRGBA("orchid") # 
  cmap[5] = cnameToRGBA("aqua") # 
  cmap[6] = cnameToRGBA("orange") # 
  cmap[7] = cnameToRGBA("royalblue") # 
  cmap[8] = cnameToRGBA("fuchsia") # 
  cmap[9] = cnameToRGBA("darkgrey") # beyond maximum clump number

  color = cmap[cid]
  pt2size = cfg.pt2size
  return (color, pt2size)


class GIplotConfig(object):
  def __init__(self):
    # the plot dimensions
    self.xinch = 3.2
    self.yinch = 1.8
    self.dpi = 400

    self.ax = [-1.0e10, 1.0e10, -1.5e10, 1.5e10]

    # define the X/Y axes of the plot as dimension numbers
    self.X = 1
    self.Y = 0
    self.Z = 2

    self.filt = "clumpnegz"

    self.ds = 5.e7
    self.pt2size = 4.00

    self.axisbg = "black"

    self.plotclmp = True
    self.plottraj = False
    self.clmptrajdt = 3600
    
    self.time = 0.
    self.G    = 6.67428e-8

    self.Mearth = 5.9736e27
    self.Mmin    = -float('inf')
    self.MminPlt = 1.e-4*self.Mearth
    self.MminLbl = 5.e-3*self.Mearth

    self.colorFunc = colorBodyAndMat

    self.clpc_ec = 'white'
    self.clpc_fc = 'none'
    self.clpc_lw = 0.2
    self.clpc_al = 0.3

    self.clpp_ec = 'white'
    self.clpp_fc = 'none'
    self.clpp_lw = 0.2
    self.clpp_al = 0.6
    
    self.clpc_txtc = 'white'
    self.clpc_txts = 4

    self.scal_fc = 'white'
    self.scal_vc = [0.1, 0.1, 0.3, 0.1]
    self.scal_ut = 'm'
    self.scal_sc = 0.01
    self.scal_txts = 4
    
    self.time_fc = 'white'
    self.time_vc = [0.1, 0.8]
    self.time_ut = 'h'
    self.time_sc = 0.000277777777
    self.time_txts = 4

    self.parm_fc = 'white'
    self.parm_vc = [0.5, 0.8]
    self.parm_txts = 4
    self.parm_txt = ""
    self.copy_vc = [0.6, 0.1]
    self.copy_txt = "Andreas Reufer, University of Bern"

    self.idfilt = 2.e6

    self.imgext = ".png"

    self.prePlot = plotGrid
    self.postPlot = plotDummy


class GIplot(object):
  def __init__(self, cfg=GIplotConfig()):
    self.fig = plt.figure( figsize=(cfg.xinch, cfg.yinch) )
    
    self.physa = plt.axes( [0.0, 0.0, 1.0, 1.0], axisbg=cfg.axisbg)
    self.norma = plt.axes( [0.0, 0.0, 1.0, 1.0], frameon=False)
    self.norma.set_xticks([])
    self.norma.set_yticks([])

    self.cfg = cfg
    
    ax = cfg.ax
    xcen = ( ax[0] + ax[1] ) / 2.
    xscl = ( ax[1] - ax[0] ) / cfg.xinch
    ycen = ( ax[2] + ax[3] ) / 2.
    yscl = ( ax[3] - ax[2] ) / cfg.yinch

    scl = max( xscl, yscl )
    corrax = [ xcen - 0.5*scl*cfg.xinch, xcen + 0.5*scl*cfg.xinch, \
        ycen - 0.5*scl*cfg.yinch, ycen + 0.5*scl*cfg.yinch ]

    self.scl = scl
    self.corrax = corrax
  
  def plotParticles(self,pfile,cfile, ifile):
    fig = self.fig

    cfg = self.cfg
    pax = self.physa
    nax = self.norma
    
    #print "loading particles ..."
    pdumpf = H5PartDump(pfile)
    sname = (pdumpf.getStepNames())[0]

    pdump = pdumpf.getStep(sname)
    nop = (pdump.m.shape)[0]

    #print "load time and G ..."
    G     = pdumpf.getAttr(sname,"gravconst")
    time  = pdumpf.getAttr(sname,"time")
    
    #print "loading clumps ..."
    cdumpf = H5PartDump(cfile)
    cdump = cdumpf.getStep(sname)
    noc = (cdump.m.shape)[0]
    
    clumps = FewBodies( cdump, G, cfg.Mmin)
    cposz = np.zeros( max( cdump.id[:,0] ) + 1)
    for i in range(noc):
      cposz[ int(cdump.id[i]) ] = cdump.pos[i,cfg.Z]

    #print "selecting particles ..."
    h = pdump.h
    cid = pdump.clumpid
    
    hmed = np.median(h)
    self.hmed = hmed

    # auto determine the scatter ptsize^2 (no idea why we need a fudge factor)
    cfg.pt2size = 0.008*np.power( hmed*cfg.dpi / self.scl , 2.)

    filt = np.ones(nop, dtype=np.bool)
    if cfg.filt == "negzonly":
      filt = pdump.pos[:,cfg.Z] < 0.
    
    if cfg.filt == "clumpcut":
      for i in range(nop):
        cidi = int( pdump.clumpid[i,:] )
        if (cidi == 0 or cidi >= noc):
          filt[i] = True
        else:
          cidi = int( pdump.clumpid[i,:] )
          zrel = pdump.pos[i,cfg.Z] - cposz[cidi]
          filt[i] = ( ( zrel > -hmed ) & ( zrel < hmed ) )
    
    if cfg.filt == "clumpnegz":
      for i in range(nop):
        cidi = int( pdump.clumpid[i,:] )
        if (cidi == 0 or cidi >= noc):
          filt[i] = True
        else:
          cidi = int( pdump.clumpid[i,:] )
          filt[i] = ( pdump.pos[i,cfg.Z] < cposz[cidi] )
    
    if cfg.filt == "slice":
      filt = ( pdump.pos[:,cfg.Z] > -hmed ) & ( pdump.pos[:,cfg.Z] <  hmed )

        
    #print "filtering particles ..."
    pos = pdump.pos[:,:].compress(filt, axis=0)[:,:]
    mat = pdump.mat[:,0].compress(filt)[:]
    bod = np.int32( (pdump.id[:].compress(filt)[:] > cfg.idfilt) )

    #print "z-sorting points ...   "
    sidx = pos[:,cfg.Z].argsort()
    
    if cfg.plottraj:
      #print "integrate clumps ...   "
      clumps.integrate(dt, cfg.ds)

    #print "pre-plot ... "
    cfg.prePlot(self.norma, self.physa, plt, pdump, cdump, cfg)
    
    #print "plotting clumps with trajectories ... "
    for i in range(1,noc):
      curtraj = clumps.traj[i]
      if ( clumps.m[i] > cfg.MminPlt ):
        pax.add_patch( mp.patches.Circle((curtraj[0,cfg.X], curtraj[0,cfg.Y]),\
          radius=clumps.rc[i], ec=cfg.clpc_ec, fc=cfg.clpc_fc, \
          lw=cfg.clpc_lw, alpha=cfg.clpc_al) )
    
        if ( clumps.m[i] > cfg.MminLbl ):
          pax.text( curtraj[0,cfg.X], curtraj[0,cfg.Y], \
              '$\mathrm{'+ '%1.4f' % ( clumps.m[i] / cfg.Mearth ) +' M_{E}}$', \
              size=cfg.clpc_txts, color=cfg.clpc_txtc)
    
        if cfg.plottraj:
          pax.add_patch( \
              mp.patches.Circle((curtraj[-1,cfg.X], curtraj[-1,cfg.Y]), \
              radius=clumps.rc[i], ec=cfg.clpp_ec, fc=cfg.clpp_fc, \
              lw=cfg.clpp_lw, alpha=cfg.clpp_al) )
          pax.plot( curtraj[:,cfg.X], curtraj[:,cfg.Y], color=cfg.clpt_fc, \
              lw=cfg.clpt_lw, alpha=cfg.clpt_al)

    #print "get point colors and size"
    (pcol, pt2size) = cfg.colorFunc(pdump, cdump, filt, cfg)
    
    #print "plotting points ...   "
    sct = pax.scatter( pos[sidx,cfg.X], pos[sidx,cfg.Y], pt2size, \
        pcol[sidx,:], lw=0, vmin=0., vmax=1.)

    #self.sidx = sidx
    #self.pos  = pos

    pax.axis("scaled")
    pax.axis(self.corrax)
   
    # plot a scale
    #print "plot a scale    ...   "
    nax.plot( [cfg.scal_vc[0], cfg.scal_vc[2]], \
        [cfg.scal_vc[1], cfg.scal_vc[3]], lw=.3, color=cfg.scal_fc )
    
    nax.axis([0., 1., 0., 1.])
    
    #print "plot time ..."
    timetxt = '$t = %6.2f' % (time*cfg.time_sc) + ' ' + cfg.time_ut + '$'
    nax.text( cfg.time_vc[0], cfg.time_vc[1], timetxt, \
        color=cfg.time_fc, size=cfg.time_txts)
    
    scldst = cfg.scal_vc[2] - cfg.scal_vc[0]
    scltxt = '$%6.1e' % (scldst*self.scl*cfg.yinch*cfg.scal_sc) \
        + "} "+ cfg.scal_ut + '$'
    scltxt = scltxt.replace('e+0', 'e+')
    scltxt = scltxt.replace('e+', '*10^{')
    nax.text( cfg.scal_vc[0], cfg.scal_vc[1] + 0.01, scltxt, \
        color=cfg.scal_fc, size=cfg.scal_txts)

    #print "parameters and copyright ..."
    self.norma.text( cfg.parm_vc[0], cfg.parm_vc[1], cfg.parm_txt, \
        color=cfg.parm_fc, size=cfg.parm_txts )
    
    self.norma.text( cfg.copy_vc[0], cfg.copy_vc[1], cfg.copy_txt, \
        color=cfg.parm_fc, size=cfg.parm_txts )
    
    #print "post plot stuff ..."
    cfg.postPlot(self.norma, self.physa, plt, pdump, cdump, cfg)

    # save the figure
    #print "save figure     ...   "
    plt.rc('savefig', dpi=self.cfg.dpi)
    plt.savefig(ifile)

  




