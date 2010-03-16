#!/usr/bin/env ipython

import matplotlib as mp
mp.use('Agg')
import tables as pt
import numpy as np
import pylab as pl
import sys
import orbit


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

dumph.close()

Mearth = 5.9736e27
Mmin = 1.e-4*Mearth
Mmin = 1.e23

#print "load clumps ...        "
#cpos = {}
clmph = pt.openFile(sys.argv[2], "r")

rotYZtoZY = np.array([[1.,0.,0.],[0.,0.,1.],[0.,1.,0.]])
clmph.root.current.pos = np.dot( clmph.root.current.pos[:,:], rotYZtoZY )
clmph.root.current.vel = np.dot( clmph.root.current.vel[:,:], rotYZtoZY )

#satSystem = orbit.SatSystem(clmph.root.current, G, Mmin)
clmph.close()

#print "keep ",clumps.rc.shape[0]," clumps"
#cdump = clmph.root.current
#for i in range(cdump.id.shape[0]):
#  cpos[int(cdump.id[i])] = cdump.pos[i,:]
#clmph.close()

#for sat in satSystem.sats:
#  a.add_patch( mp.patches.Circle(

#print "plotting clumps with trajectories ... "
#for i in range(clumps.noc):
#  curtraj = clumps.traj[i]
#  a.add_patch( mp.patches.Circle((curtraj[0,2], curtraj[0,0]), radius=clumps.rc[i], ec='yellow', fc='none', lw=0.3, alpha=0.3) )
#  a.add_patch( mp.patches.Circle((curtraj[-1,2], curtraj[-1,0]), radius=clumps.rc[i], ec='green', fc='none', lw=0.3, alpha=0.6) )
#  a.plot( curtraj[:,2], curtraj[:,0], color='white', lw=0.3, alpha=0.3)
#  a.text( curtraj[0,2], curtraj[0,0], '$\mathrm{'+ '%1.4f' % ( clumps.m[i] / Mearth ) +' M_{\earth}}$', size=3, color='yellow')
#  
  #a.text( curtraj[0,2], curtraj[0,0], ' '+ '%1.4f' % ( clumps.m[i] / Mearth ) +' Me', size=3, color='yellow')


pl.figure( figsize=(3, 2) )
#a = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg='k')
a = pl.axes( [0.1, 0.1, 0.9, 0.9])
pl.grid(True)

a.axis("scaled")
a.axis([-8.e9, 4.e9, -5.e9, 3.e9])

#for cursat in satSystem.sats:
#  cursat.plotorbit(a)
#  a.add_patch( mp.patches.Circle((cursat.r0m[0], cursat.r0m[1]), radius=cursat.rc, ec='green', fc='none', lw=1.0, alpha=0.6) )


M = 1.
m = 0.
rM = np.array( [0.0, 0., 0.] )
#rm = np.array( [1.5, 0., 0.0] )
#rm = np.array( [0.7071, 0.0000, -0.7071] )
#rm = np.array( [0.7071, 0.7071, -0.0000] )
#rm = np.array( [-0.2071, 0.0000, 0.7071] )
rm = np.array( [ 1.0000, 0.0000, 0.0000] )

vM = np.array( [0., 0., 0.] )
#vm = np.array( [0., 1.00, 0.] )
#vm = np.array( [0., -1.00, 0.] )
#vm = np.array( [-0.7071,  0.7071, 0.] )
vm = np.array( [0., 1.224745, 0.] )

G = 1.
rc = 0.1

a.axis([-3.0, 3.0, -2.0, 2.0])
mysat = orbit.Satellite(rm, vm, m, rM, vM, M, rc, G)
(rcart, r, R) = mysat.plotorbit(a)
a.add_patch( mp.patches.Circle((mysat.r0m[0], mysat.r0m[1]), radius=mysat.rc, ec='green', fc='none', lw=0.3, alpha=0.6) )
a.plot(rM[0], rM[1], 'x')
print "mysat i = ",(mysat.i / np.pi)," e = ",(mysat.e), " argp = ",(mysat.argp/np.pi), " longa = ",(mysat.longa / np.pi )

pl.rc('savefig', dpi=400)
pl.savefig(sys.argv[3])


