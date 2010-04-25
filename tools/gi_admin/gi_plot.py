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

# the plot dimensions
xinch = 3
yinch = 2

# define the X/Y axes of the plot as dimension numbers
X = 1
Y = 0
Z = 2

Mearth = 5.9736e27
Mmin = 1.e-4*Mearth
ds = 5.e7
ptsize = 0.10

# FIXME: find this value automatically
zmin = -7.e7
zmax =  7.e7

cmap = np.zeros((2,6), float)
cmap[0,2] = 0.30 # ice
cmap[0,4] = 0.95 # dun
cmap[0,5] = 0.05 # iron
cmap[1,2] = 0.35
cmap[1,4] = 0.90
cmap[1,5] = 0.10


class GIplot(object):
  def __init__(self, pfile, cfile, log=sys.stdout):
    self.pfile = pfile
    self.cfile = cfile

    pl.figure( figsize=(xinch, yinch) )
    self.a = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg='k')

    #dumph = pt.openFile(sys.argv[1], "r")
    dumph = pt.openFile(pfile, "r")
    dump  = dumph.root.current

    attrlist = dump._v_attrs
    attrnams = attrlist._v_attrnames

    self.time = 0.
    if attrnams.__contains__("time"):
      self.time = attrlist["time"]

    self.G = 6.67428e-8
    if attrnams.__contains__("gravconst"):
      self.G    = attrlist["gravconst"]
    dumph.close()
    
    #print "load clumps ...        "
    self.cpos = {}
    clmph = pt.openFile(cfile, "r")
    self.clumps = fewbody.FewBodies( clmph.root.current, self.G, Mmin)
    #print "keep ",self.clumps.rc.shape[0]," clumps"
    cdump = clmph.root.current
    for i in range(cdump.id.shape[0]):
      self.cpos[int(cdump.id[i])] = cdump.pos[i,:]
    clmph.close()

  def plotClumps(self, dt):
    clumps = self.clumps

    #print "integrate clumps ...   "
    clumps.integrate(dt, ds)
    
    #print "plotting clumps with trajectories ... "
    for i in range(clumps.noc):
      curtraj = clumps.traj[i]
      self.a.add_patch( mp.patches.Circle((curtraj[0,X], curtraj[0,Y]), radius=clumps.rc[i], ec='yellow', fc='none', lw=0.3, alpha=0.3) )
      self.a.add_patch( mp.patches.Circle((curtraj[-1,X], curtraj[-1,Y]), radius=clumps.rc[i], ec='green', fc='none', lw=0.3, alpha=0.6) )
      self.a.plot( curtraj[:,X], curtraj[:,Y], color='white', lw=0.3, alpha=0.3)
      self.a.text( curtraj[0,X], curtraj[0,Y], '$\mathrm{'+ '%1.4f' % ( clumps.m[i] / Mearth ) +' M_{E}}$', size=4, color='yellow')


  def plotParticles(self):
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
    bod = np.int32( (dump.id[:].compress(filt)[:] > 2.e6) )
    dumph.close()

    #print "z-sorting points ...   "
    color = cmap[bod, mat]
    sidx = pos[:,Z].argsort()
    
    #print "plotting points ...   "
    self.a.scatter( pos[sidx,X], pos[sidx,Y], ptsize, color[sidx], lw=0)

  
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

    pl.rc('savefig', dpi=400)
    pl.savefig(file)



