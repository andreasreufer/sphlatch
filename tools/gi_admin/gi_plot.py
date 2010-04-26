#!/usr/bin/env ipython

import tables as pt
import numpy as np
import matplotlib as mp
mp.use('Agg')
import pylab as pl
import sys
import fewbody


#pl.rc('text', usetex=True)
pl.rc('text.latex', preamble = '\usepackage{amssymb}, \usepackage{wasysym}')



class GIplotConfig(object):
  def __init__(self):
    # the plot dimensions
    self.xinch = 3
    self.yinch = 2

    # define the X/Y axes of the plot as dimension numbers
    self.X = 1
    self.Y = 0
    self.Z = 2

    # FIXME: find this value automatically
    self.zmin = -7.e7
    self.zmax =  7.e7

    self.ds = 5.e7
    self.ptsize = 0.10

    cmap = np.zeros((2,6), float)
    cmap[0,2] = 0.30 # ice
    cmap[0,4] = 0.95 # dun
    cmap[0,5] = 0.05 # iron
    cmap[1,2] = 0.35
    cmap[1,4] = 0.90
    cmap[1,5] = 0.10
    self.cmap = cmap

    self.axisbg = "k"

    self.plottraj = False
    self.dpi = 400

    self.time = 0.
    self.G    = 6.67428e-8

    self.Mearth = 5.9736e27
    self.Mmin = 1.e-4*self.Mearth
      
    self.clpc_ec = 'yellow'
    self.clpc_fc = 'none'
    self.clpc_lw = 0.3
    self.clpc_al = 0.3

    self.clp_txtc = 'yellow'
    self.clp_txts = 4

    self.clpp_ec = 'green'
    self.clpp_fc = 'none'
    self.clpp_lw = 0.3
    self.clpp_al = 0.6

    self.idfilt = 2.e6


class GIplot(object):
  def __init__(self, pfile, cfile, log=sys.stdout, cfg=GIplotConfig):
    self.pfile = pfile
    self.cfile = cfile
    self.cfg = cfg

    pl.figure( figsize=(cfg.xinch, cfg.yinch) )
    self.a = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg=cfg.axisbg)

    #dumph = pt.openFile(sys.argv[1], "r")
    dumph = pt.openFile(pfile, "r")
    dump  = dumph.root.current

    attrlist = dump._v_attrs
    attrnams = attrlist._v_attrnames

    self.time = cfg.time
    if attrnams.__contains__("time"):
      self.time = attrlist["time"]

    self.G = cfg.G
    if attrnams.__contains__("gravconst"):
      self.G    = attrlist["gravconst"]
    dumph.close()
    
    #print "load clumps ...        "
    self.cpos = {}
    clmph = pt.openFile(cfile, "r")
    self.clumps = fewbody.FewBodies( clmph.root.current, self.G, cfg.Mmin)
    cdump = clmph.root.current
    for i in range(cdump.id.shape[0]):
      self.cpos[int(cdump.id[i])] = cdump.pos[i,:]
    clmph.close()

  def plotClumps(self, dt):
    cfg    = self.cfg
    clumps = self.clumps

    #print "integrate clumps ...   "
    if cfg.plottraj:
      clumps.integrate(dt, cfg.ds)
    
    #print "plotting clumps with trajectories ... "
    for i in range(clumps.noc):
      curtraj = clumps.traj[i]
      self.a.add_patch( mp.patches.Circle((curtraj[0,X], curtraj[0,Y]), \
          radius=clumps.rc[i], ec=cfg.clpc_ec, fc=cfg.clpc_fc, \
          lw=cfg.clpc_lw, alpha=cfg.clpc_al) )
      
      self.a.text( curtraj[0,X], curtraj[0,Y], \
          '$\mathrm{'+ '%1.4f' % ( clumps.m[i] / Mearth ) +' M_{E}}$', \
          size=cfg.clpc_txts, color=cfg.clpc_txtc)
    
      if cfg.plottraj:
        self.a.add_patch( mp.patches.Circle((curtraj[-1,X], curtraj[-1,Y]), \
            radius=clumps.rc[i], ec=cfg.clpp_ec, fc=cfg.clpp_fc, \
            lw=cfg.clpp_lw, alpha=cfg.clpp_al) )
        self.a.plot( curtraj[:,X], curtraj[:,Y], color=cfg.clpt_fc, \
            lw=cfg.clpt_lw, alpha=cfg.clpt_al)


  def plotParticles(self):
    cfg = self.cfg
    cpos = self.cpos
    
    dumph = pt.openFile(self.pfile, "r")
    dump  = dumph.root.current
    
    #print "selecting particles ..."
    h = dump.h
    cid = dump.clumpid
    nop = dump.pos.shape[0]
    filt = np.zeros(nop, dtype=np.bool)
    for i in range(nop):
      if (dump.clumpid[i,:] == 0):
        filt[i] = True
      else:
        cidi = int( dump.clumpid[i,:] )
        zrel = dump.pos[i,Z] - cpos[cidi][Z]
        filt[i] = ( ( zrel > zmin ) & ( zrel < zmax ) )
        #filt[i] = ( ( dump.pos[i,Z] > ymin ) & ( dump.pos[i,Z] < ymax ) )

    #print "filtering particles ..."
    pos = dump.pos[:,:].compress(filt, axis=0)[:,:]
    mat = dump.mat[:,0].compress(filt)[:]
    bod = np.int32( (dump.id[:].compress(filt)[:] > cfg.idfilt) )
    dumph.close()

    #print "z-sorting points ...   "
    color = cmap[bod, mat]
    sidx = pos[:,Z].argsort()
    
    #print "plotting points ...   "
    self.a.scatter( pos[sidx,X], pos[sidx,Y], cfg.ptsize, color[sidx], lw=0)
    pl.clim(0.0, 1.0)

  
  def storePlot(self, file, ax):
    xcen = ( ax[0] + ax[1] ) / 2.
    xscl = ( ax[1] - ax[0] ) / xinch
    ycen = ( ax[2] + ax[3] ) / 2.
    yscl = ( ax[3] - ax[2] ) / yinch

    scl = max( xscl, yscl )
    corrax = [ xcen - 0.5*scl*xinch, xcen + 0.5*scl*xinch, \
        ycen - 0.5*scl*yinch, ycen + 0.5*scl*yinch ]

    self.a.axis("scaled")
    self.a.axis(corrax)

    pl.rc('savefig', dpi=self.cfg.dpi)
    pl.savefig(file)



