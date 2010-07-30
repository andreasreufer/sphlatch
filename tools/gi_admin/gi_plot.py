#!/usr/bin/env python
import sys
import numpy as np

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

from fewbody import FewBodies
from h5part import H5PartDump

#mp.rc('text', usetex=True)
mp.rc('text.latex', preamble = '\usepackage{amssymb}, \usepackage{wasysym}')

class GIplotConfig(object):
  def __init__(self):
    # the plot dimensions
    self.xinch = 3
    self.yinch = 2
    self.dpi = 400

    # 
    self.ax = [-1.0e10, 1.0e10, -1.5e10, 1.5e10]

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

    self.plotclmp = True
    self.plottraj = False
    self.clmptrajdt = 3600
    
    self.time = 0.
    self.G    = 6.67428e-8

    self.Mearth = 5.9736e27
    self.Mmin    = 0.e-0*self.Mearth
    self.MminLbl = 1.e-3*self.Mearth

    self.clpc_ec = 'yellow'
    self.clpc_fc = 'none'
    self.clpc_lw = 0.3
    self.clpc_al = 0.3

    self.clpc_txtc = 'yellow'
    self.clpc_txts = 4

    self.clpp_ec = 'green'
    self.clpp_fc = 'none'
    self.clpp_lw = 0.3
    self.clpp_al = 0.6

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
    self.parm_vc = [0.8, 0.8]
    self.parm_txts = 4
    self.copy_vc = [0.6, 0.1]
    self.copy_txt = '$\mathrm{Andreas Reufer, University of Bern}$'

    self.idfilt = 2.e6

    self.imgext = ".png"


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

  def plotTime(self, time):
    cfg    = self.cfg
    timetxt = '$t = %6.2f' % (time*cfg.time_sc) + ' ' + cfg.time_ut + '$'
    self.norma.text( cfg.time_vc[0], cfg.time_vc[1], timetxt, \
        color=cfg.time_fc, size=cfg.time_txts)

  def plotSimParams(self, params):
    cfg    = self.cfg
    pos = cfg.parm_vc
    col = cfg.parm_fc
    siz = cfg.parm_txts

    parmtxt = '$ haba haba2 haba3 $'
    self.norma.text( pos[0], pos[1], parmtxt, \
        color=cfg.parm_fc, size=cfg.parm_txts )
    
    poscopy = cfg.copy_vc
    self.norma.text( poscopy[0], poscopy[1], cfg.copy_txt, \
        color=cfg.parm_fc, size=cfg.parm_txts )

  def plotClumps(self,cfile):
    cfg    = self.cfg
    clumps = self.clumps
    pax = self.physa

    #print "integrate clumps ...   "
    if cfg.plottraj:
      clumps.integrate(dt, cfg.ds)
    
    #print "plotting clumps with trajectories ... "
    for i in range(clumps.noc):
      curtraj = clumps.traj[i]
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

  def plotParticles(self,pfile,cfile, ifile):
    fig = self.fig

    cfg = self.cfg

    print "loading particles ..."
    pdumpf = H5PartDump(pfile)
    sname = (pdumpf.getStepNames())[0]

    pdump = pdumpf.getStep(sname)
    nop = (pdump.m.shape)[0]

    print "load time and G ..."
    G     = pdumpf.getAttr(sname,"gravconst")
    time  = pdumpf.getAttr(sname,"time")

    print "loading clumps ..."
    cdumpf = H5PartDump(cfile)
    cdump = cdumpf.getStep(sname)
    noc = (cdump.m.shape)[0]
    
    clumps = FewBodies( cdump, G, cfg.Mmin)
    cposz = np.zeros( max( cdump.id[:,0] ) + 1)
    for i in range(noc):
      cposz[ int(cdump.id[i]) ] = cdump.pos[i,cfg.Z]
    
    print "selecting particles ..."
    h = pdump.h
    cid = pdump.clumpid

    filt = np.zeros(nop, dtype=np.bool)
    for i in range(nop):
      cidi = int( pdump.clumpid[i,:] )
      if (cidi == 0 or cidi >= noc):
        filt[i] = True
      else:
        cidi = int( pdump.clumpid[i,:] )
        zrel = pdump.pos[i,cfg.Z] - cposz[cidi]
        filt[i] = ( ( zrel > cfg.zmin ) & ( zrel < cfg.zmax ) )
        
    print "filtering particles ..."
    pos = pdump.pos[:,:].compress(filt, axis=0)[:,:]
    mat = pdump.mat[:,0].compress(filt)[:]
    bod = np.int32( (pdump.id[:].compress(filt)[:] > cfg.idfilt) )

    print "z-sorting points ...   "
    color = cfg.cmap[bod, mat]
    sidx = pos[:,cfg.Z].argsort()
    
    print "plotting points ...   "
    self.physa.scatter( pos[sidx,cfg.X], pos[sidx,cfg.Y], cfg.ptsize, \
        color[sidx], lw=0, vmin=0., vmax=1.)

    self.physa.axis("scaled")
    self.physa.axis(self.corrax)
   
    # plot a scale
    print "plot a scale    ...   "
    self.norma.plot( [cfg.scal_vc[0], cfg.scal_vc[2]], \
        [cfg.scal_vc[1], cfg.scal_vc[3]], lw=.3, color=cfg.scal_fc )
    
    self.norma.axis([0., 1., 0., 1.])
    
    scldst = cfg.scal_vc[2] - cfg.scal_vc[0]
    scltxt = '$%6.1e' % (scldst*self.scl*cfg.yinch*cfg.scal_sc) \
        + "} "+ cfg.scal_ut + '$'
    scltxt = scltxt.replace('e+0', 'e+')
    scltxt = scltxt.replace('e+', '*10^{')
    self.norma.text( cfg.scal_vc[0], cfg.scal_vc[1] + 0.01, scltxt, \
        color=cfg.scal_fc, size=cfg.scal_txts)

    # save the figure
    print "save figure     ...   "
    plt.rc('savefig', dpi=self.cfg.dpi)
    plt.savefig(ifile)
  




