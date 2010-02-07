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

pl.figure( figsize=(4, 2) )
a = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg='k')

dt = 3600*5
ds = 5.e7

clumps.integrate(dt, ds)

print 'plotting clumps ......'
for i in range(clumps.noc):
  
  curtraj = clumps.traj[i]

  # plot a circle at the first position
  a.add_patch( mp.patches.Circle((curtraj[0,2], curtraj[0,0]), radius=clumps.rc[i], ec='yellow', fc='none', lw=0.3, alpha=0.3) )
  
  # plot a circle at the last position
  a.add_patch( mp.patches.Circle((curtraj[-1,2], curtraj[-1,0]), radius=clumps.rc[i], ec='green', fc='none', lw=0.3, alpha=0.6) )
 
  # plot trajectory
  a.plot( curtraj[:,2], curtraj[:,0], color='white', lw=0.3, alpha=0.3)
  a.text( curtraj[0,2], curtraj[0,0], '$\mathrm{'+str( clumps.m[i] )+' g}$', size=3, color='white')

ymin = -5.e9
ymax =  0.e9

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
a.scatter( pos[:,2], pos[:,0], 0.15, color, linewidth=0)
fileh.close()

a.axis("scaled")
a.axis([-8.e9, 4.e9, -4.e9, 2.e9])


b = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg='none', frameon=False)
b.text(0.1, 0.05, 'HABA', color='b')
#b.text(0.1, 0.10, '$\mathrm{Giant~impact}$', color='b', size=8)
b.text(0.1, 0.10, "$\\mathrm{x [cm]}$", color='w', size=8)
b.axis([0.,1.,0.,1.])


pl.rc('savefig', dpi=600)
pl.savefig("integrate.png")



#clumps.advance(100)

#for i in range(clumps.noc):
  #pl.gca().add_patch( mp.patches.Circle((clumps.pos[i,2], clumps.pos[i,0]), radius=clumps.rc[i], ec='green', fc='none', lw=0.3, alpha=0.6) )

