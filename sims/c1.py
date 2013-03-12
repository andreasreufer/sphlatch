#!/usr/bin/env ipython

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

from sim_admin import SimAdmin, SimParam, SimSetConfig
from simulation import resolvePath
from gi_plot import *
from gi_viz  import GIviz, GIvizConfig
from clumps import ClumpsPlotConfig
import os
import stat
import time
import numpy as np
import commands
import socket
import datetime

nan = float('nan')
lem = 3.40000e41

sscfg = SimSetConfig()
sscfg.name = "c1"

simadm = SimAdmin(sscfg)
sims = simadm._sims

eVinK = 11604.505

nofill = [0.5, 0.5, 0.5, 0.0]
fig = plt.figure()
na = plt.axes( [0.0, 0.0, 1.0, 1.0] )
ax = fig.add_subplot(111)


#for sim in sims.values():
#for sim in simadm.getSimsByFilter( SimParam(nan, 1.000, nan, nan) ):
for sim in simadm.getSimsByFilter( SimParam(nan, 0.100, nan, nan) ):
  res = sim.results

  if res.noc > 1:
    mdisk = sim.results.mdiskmatm[1,1] / ML
    Ltot  = sim.results.Lm[1,2] / lem

    if mdisk > 1.0:
      print (sim.params.key, mdisk, Ltot)

    mrk = "x"
    if sim.params.mimp == 0.02:
      mrk = "d"
    elif sim.params.mimp == 0.01:
      mrk = "o"
    elif sim.params.mimp == 0.035:
      mrk = "p"


    ax.scatter( Ltot, mdisk, marker=mrk)

#ax.scatter( 1.0, 2.0, marker='x' )
ax.set_xlabel(r"$L / L_{E-M}$",color="black")
ax.set_ylabel(r"$M_{disk} / M_L$",color="black")

#ax.axis( [0., 4.1, 0., 5.5 ] )
ax.axis( [0., 0.41, 0., 0.5 ] )

#ax.xaxis.set_major_formatter(math_formatter)
ax.yaxis.set_major_formatter(math_formatter)

ax.grid(True)

fig.savefig("out.pdf")

