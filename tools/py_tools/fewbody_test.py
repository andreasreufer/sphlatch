import tables
import fewbody

G =  6.67e-08
mmin = 5.e23

fileh = tables.openFile("dump0010494_T03.0000e+04_clumps.h5part", "r")
clumps = fewbody.FewBodies( fileh.root.current, G, mmin)
fileh.close()

import matplotlib as mp
import pylab as pl
import numpy as np

pl.figure( figsize=(4, 2) )
a = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg='k')

dt = 3600
ds = 1.e7

print 'plotting clumps ......'
for i in range(clumps.noc):
  pl.gca().add_patch( mp.patches.Circle((clumps.pos[i,2], clumps.pos[i,0]), radius=clumps.rc[i], ec='yellow', fc='none', lw=0.3, alpha=0.3) )
  pl.gca().add_patch( mp.patches.Arrow(clumps.pos[i,2], clumps.pos[i,0], dt*clumps.vel[i,2], dt*clumps.vel[i,0], ec='white', lw=0.3, alpha=0.3, hatch='/+*') )

#clumps.advance(100)

for i in range(clumps.noc):
  pl.gca().add_patch( mp.patches.Circle((clumps.pos[i,2], clumps.pos[i,0]), radius=clumps.rc[i], ec='green', fc='none', lw=0.3, alpha=0.6) )

