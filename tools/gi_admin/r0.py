#!/usr/bin/env ipython

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

from sim_admin import SimAdmin, SimParam, SimSetConfig
from gi_plot import GIplot, GIplotConfig, colorTemp, colorDens, colorClumps, plotGrid, plotSimParams
from gi_viz  import GIviz, GIvizConfig
from clumps import ClumpsPlotConfig
import os
import numpy as np

nan = float('nan')

ssconf = SimSetConfig()
ssconf.name = "r0"
ssconf.dir = "./" + ssconf.name
#ssconf.bodiesdb = "~/pmc_colls/bodies"
ssconf.bodiesdb = "./bodies"
ssconf.T = 1500.

ssconf.mimp = np.array( [0.1, 1.0] )
ssconf.mtar = np.array( [0.1, 1.0] )
ssconf.vimp = np.array( [1.0, 1.1, 1.15, 1.2, 1.25, 1.3, 1.5] )
ssconf.impa = np.array( [0, 15, 30, 45, 60, 75, 80, 85] )

ssconf.attr.append( ("UMIN",1.e9) )
ssconf.auxfiles = "~/pmc_colls/aux/aneos.input"
ssconf.srcdir   = "~/repos/sphlatch/apps/simple_sph"
ssconf.maketarg = "simple_sph_GHUaC_"
ssconf.binary   = "simple_sph_GHUaC_"

ssconf.subcmd  = "qsub -M andreas.reufer@space.unibe.ch -b y -cwd -l h_cpu=240:00:00 -l h_vmem=512M -pe smp $NOCPUS -R y -N $SIMNAME ./$BINARY initial.h5part $SAVETIME $STOPTIME $RUNARGS"
ssconf.delcmd  = "qdel $JOBID"
ssconf.nocpus  = 4
ssconf.runargs = "dump 8"

vizcfg = GIvizConfig()
vizcfg.subcmd = "qsub -M andreas.reufer@space.unibe.ch -b n -cwd -l h_cpu=4:00:00 -l h_vmem=512M -N $JOBNAME $JOBCMD"

axnorm = [-1.6e10, 1.6e10, -0.9e10, 0.9e10]
axzoom = [-3.2e9 , 3.2e9 , -1.8e9 , 1.8e9 ]

cfgmatXY = GIplotConfig()
cfgmatXY.ax = axnorm
cfgmatXY.postPlot  = plotSimParams
cfgmatXY.imgdir = ssconf.dir + "/viz_matXY"

cfgmatYZ = GIplotConfig()
cfgmatYZ.ax = axnorm
cfgmatYZ.X = 1
cfgmatYZ.Y = 2
cfgmatYZ.Z = 0
cfgmatYZ.iname = "matYZ_"
cfgmatYZ.postPlot  = plotSimParams
cfgmatYZ.imgdir = ssconf.dir + "/viz_matYZ"

cfgcidXY = GIplotConfig()
cfgcidXY.ax = axnorm
cfgcidXY.iname = "cidXY_"
cfgcidXY.colorFunc = colorClumps
cfgcidXY.postPlot  = plotSimParams
cfgcidXY.imgdir = ssconf.dir + "/viz_cidXY"

cfgcidYZ = GIplotConfig()
cfgcidYZ.ax = axnorm
cfgcidYZ.X = 1
cfgcidYZ.Y = 2
cfgcidYZ.Z = 0
cfgcidYZ.iname = "cidYZ_"
cfgcidYZ.colorFunc = colorClumps
cfgcidYZ.postPlot  = plotSimParams
cfgcidYZ.imgdir = ssconf.dir + "/viz_cidYZ"

cfgTnorm = GIplotConfig()
cfgTnorm.ax = axnorm
cfgTnorm.Tmin = 1000
cfgTnorm.Tmax = 4000
cfgTnorm.colorFunc = colorTemp
cfgTnorm.iname = "Tnorm_"
cfgTnorm.postPlot  = plotSimParams
cfgTnorm.imgdir = ssconf.dir + "/viz_Tnorm"

cfgTzoom = GIplotConfig()
cfgTzoom.ax = axzoom
cfgTzoom.Tmin = 1000
cfgTzoom.Tmax = 4000
cfgTzoom.colorFunc = colorTemp
cfgTzoom.postPlot  = plotSimParams
cfgTzoom.iname = "Tzoom_"
cfgTzoom.postPlot  = plotSimParams
cfgTzoom.imgdir = ssconf.dir + "/viz_Tzoom"

cfgDensNorm = GIplotConfig()
cfgDensNorm.ax = axnorm
cfgDensNorm.rhomin = 0.
cfgDensNorm.rhomax = 6.
cfgDensNorm.colorFunc = colorDens
cfgDensNorm.iname = "densNorm_"
cfgDensNorm.postPlot  = plotSimParams
cfgDensNorm.imgdir = ssconf.dir + "/viz_DensNorm"


GIvizs = [GIviz(vizcfg, cfgmatXY), \
   GIviz(vizcfg, cfgmatYZ), \
   GIviz(vizcfg, cfgcidXY), \
   GIviz(vizcfg, cfgcidYZ), \
   GIviz(vizcfg, cfgTnorm), \
   GIviz(vizcfg, cfgTzoom), \
   GIviz(vizcfg, cfgDensNorm) ]

simadm = SimAdmin(ssconf)
sims = simadm._sims

simA = sims["mtar001.000_mimp001.000_impa85.0_vimp1.15"]
simB = sims["mtar000.100_mimp000.100_impa00.0_vimp1.25"]

simsA = simadm.getSims( SimParam(1.0, 1.0, nan, 1.15) )

def vizSim(sim):
  for viz in GIvizs:
    viz.vizSim(sim)

for sim in sims.values():
  if sim.nodumps > 0:
    print sim.params.key, sim.nodumps
    vizSim(sim)


