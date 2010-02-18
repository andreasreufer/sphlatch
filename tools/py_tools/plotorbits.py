#!/usr/bin/env python

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

Mmin = 1.e27

print "load clumps ...        "
cpos = {}
clmph = pt.openFile(sys.argv[2], "r")
sats = orbit.SatSystem(clmph.root.current, G, Mmin)

#print "keep ",clumps.rc.shape[0]," clumps"
#cdump = clmph.root.current
#for i in range(cdump.id.shape[0]):
#  cpos[int(cdump.id[i])] = cdump.pos[i,:]
#clmph.close()


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
a = pl.axes( [0.0, 0.0, 1.0, 1.0], axisbg='k')

a.axis("scaled")
a.axis([-8.e9, 4.e9, -5.e9, 3.e9])

pl.rc('savefig', dpi=400)
pl.savefig(sys.argv[3])


