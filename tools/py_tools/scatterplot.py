#!/usr/bin/env ipython

import tables as pt
import numpy as np
import matplotlib as mp
mp.use('Agg')
import pylab as pl
import sys
import fewbody

<<<<<<< HEAD
=======
#pl.rc('text', usetex=True)
pl.rc('text.latex', preamble = '\usepackage{amssymb}, \usepackage{wasysym}')

>>>>>>> sacerdos/develop
execfile("config.sh")

Mearth = 5.9736e27
Mmin = 1.e-4*Mearth
dt = 0.5*3600
ds = 5.e7
ptsize = 0.10

# FIXME: find this value automatically
ymin = -7.e7
ymax =  7.e7

cmap = np.zeros((2,6), float)
cmap[0,2] = 0.30 # ice
cmap[0,4] = 0.95 # dun
cmap[0,5] = 0.05 # iron
cmap[1,2] = 0.35
cmap[1,4] = 0.90
cmap[1,5] = 0.10


#def plottext():

#def plotdump():

dumph = pt.openFile(sys.argv[1], "r")
dump  = dumph.root.current;

attrlist = dump._v_attrs
attrnams = attrlist._v_attrnames

time = 0.
if attrnams.__contains__("time"):
  time = attrlist["time"]

G = 6.67428e-8
if attrnams.__contains__("gravconst"):
  G    = attrlist["gravconst"]

print "load clumps ...        "
cpos = {}
clmph = pt.openFile(sys.argv[2], "r")
clumps = fewbody.FewBodies( clmph.root.current, G, Mmin)
print "keep ",clumps.rc.shape[0]," clumps"
cdump = clmph.root.current
for i in range(cdump.id.shape[0]):
  cpos[int(cdump.id[i])] = cdump.pos[i,:]
clmph.close()

print "selecting particles ..."
h = dump.h
cid = dump.clumpid

nop = dump.pos.shape[0]
filt = np.zeros(nop, dtype=np.bool)
for i in range(nop):
  if (dump.clumpid[i,:] == 0):
    filt[i] = True
  else:
    cidi = int( dump.clumpid[i,:] )
    yrel = dump.pos[i,1] - cpos[cidi][1]
    filt[i] = ( ( yrel > ymin ) & ( yrel < ymax ) )
    #filt[i] = ( ( dump.pos[i,1] > ymin ) & ( dump.pos[i,1] < ymax ) )

print "filtering particles ..."
pos = dump.pos[:,:].compress(filt, axis=0)[:,:]
mat = dump.mat[:,0].compress(filt)[:]
bod = np.int32( (dump.id[:].compress(filt)[:] > 2.e6) )
dumph.close()

print "z-sorting points ...   "
color = cmap[bod, mat]
sidx = pos[:,1].argsort()

pl.figure( figsize=(3, 2) )
a = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg='k')
print "plotting points ...   "
a.scatter( pos[sidx,2], pos[sidx,0], ptsize, color[sidx], lw=0)

print "integrate clumps ...   "
clumps.integrate(dt, ds)


print "plotting clumps with trajectories ... "
for i in range(clumps.noc):
  curtraj = clumps.traj[i]
  a.add_patch( mp.patches.Circle((curtraj[0,2], curtraj[0,0]), radius=clumps.rc[i], ec='yellow', fc='none', lw=0.3, alpha=0.3) )
  a.add_patch( mp.patches.Circle((curtraj[-1,2], curtraj[-1,0]), radius=clumps.rc[i], ec='green', fc='none', lw=0.3, alpha=0.6) )
  a.plot( curtraj[:,2], curtraj[:,0], color='white', lw=0.3, alpha=0.3)
<<<<<<< HEAD
  a.text( curtraj[0,2], curtraj[0,0], '$\mathrm{'+ '%1.4f' % ( clumps.m[i] / Mearth ) +' M_{E}}$', size=4, color='yellow')

=======
  a.text( curtraj[0,2], curtraj[0,0], '$\mathrm{'+ '%1.4f' % ( clumps.m[i] / Mearth ) +' M_{earth}}$', size=3, color='yellow')
  
>>>>>>> sacerdos/develop


a.axis("scaled")
a.axis([-8.e9, 4.e9, -5.e9, 3.e9])

pl.rc('savefig', dpi=400)
pl.savefig(sys.argv[3])


