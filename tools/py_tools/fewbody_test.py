import tables
import fewbody

G =  6.67e-08
mmin = 5.e23

fileh = tables.openFile("dump0010494_T03.0000e+04_clumps.h5part", "r")
clumps = fewbody.FewBodies( fileh.root.current, G, mmin)
fileh.close()

import matplotlib as mp
mp.use('agg')
import pylab as pl
import numpy as np

#from mp.patches import Circle

#pl.scatter(clumps.pos[:,0], clumps.pos[:,2], 1.)

pl.figure( figsize=(4, 2) )
a = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg='k')

dt = 3600

print 'plotting clumps ......'
for i in range(clumps.noc):
  pl.gca().add_patch( mp.patches.Circle((clumps.pos[i,2], clumps.pos[i,0]), radius=clumps.rc[i], ec='yellow', fc='none', lw=0.3, alpha=0.3) )
  pl.gca().add_patch( mp.patches.Arrow(clumps.pos[i,2], clumps.pos[i,0], dt*clumps.vel[i,2], dt*clumps.vel[i,0], ec='white', lw=0.3, alpha=0.3, hatch='/+*') )
  #pl.gca()

for i in range(0,30):
  for j in range(0,72):
    clumps.advance(10)
  for i in range(1,clumps.noc):
    pl.gca().add_patch( mp.patches.Circle((clumps.pos[i,2], clumps.pos[i,0]), radius=clumps.rc[i], ec='white', fc='none', lw=0.1, alpha=0.3) )

ymin = -5.e9
ymax =  5.e9

cmap = np.zeros((2,6), float)
cmap[0,2] = 0.30 # ice
cmap[0,4] = 0.95 # dun
cmap[0,5] = 0.05 # iron
cmap[1,2] = 0.35
cmap[1,4] = 0.90
cmap[1,5] = 0.10

print 'plotting particles ...'

fileh = tables.openFile("dump0010494_T03.0000e+04.h5part", "r")
dump = fileh.root.current;
filt = ( ( dump.pos[:,1] > ymin ) & ( dump.pos[:,1] < ymax ) )
pos = dump.pos[:,:].compress(filt, axis=0)[:,:]
mat = dump.mat[:,0].compress(filt)[:]
bod = np.int32( (dump.id[:].compress(filt)[:] > 2.e6) )
color = cmap[bod, mat]
pl.scatter( pos[:,2], pos[:,0], 0.15, color, linewidth=0)
fileh.close()

pl.axis("scaled")
pl.axis([-8.e9, 4.e9, -4.e9, 2.e9])

pl.rc('savefig', dpi=600)
pl.savefig("out.png")


