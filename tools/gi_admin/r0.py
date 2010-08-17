#!/usr/bin/env ipython

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

from sim_admin import SimAdmin, SimParam, SimSetConfig
from gi_plot import GIplot, GIplotConfig, colorScalarLinear, colorScalarLog10, colorClumps, colorBodyAndPhaseANEOS, plotGrid, plotSimParams
from gi_viz  import GIviz, GIvizConfig
from clumps import ClumpsPlotConfig
import os
import numpy as np

nan = float('nan')

ssconf = SimSetConfig()
ssconf.name = "r0"
ssconf.dir = "./" + ssconf.name
ssconf.bodiesdb = "./bodies"
ssconf.T = 1500.

ssconf.mimp = np.array( [0.1, 1.0] )
ssconf.mtar = np.array( [0.1, 1.0] )
ssconf.vimp = np.array( [1.0, 1.1, 1.15, 1.2, 1.25, 1.3, 1.5] )
ssconf.impa = np.array( [0, 15, 30, 45, 60, 75, 80, 85] )

ssconf.attr.append( ("UMIN",1.e9) )
ssconf.auxfiles = "./aux/aneos.input"
ssconf.srcdir   = "~/repos/sphlatch/apps/simple_sph"
ssconf.maketarg = "simple_sph_GHUaC_"
ssconf.binary   = "simple_sph_GHUaC_"

ssconf.subcmd  = "qsub -M andreas.reufer@space.unibe.ch -b y -cwd -l h_cpu=240:00:00 -l h_vmem=512M -pe smp $NOCPUS -R y -N $SIMNAME ./$BINARY initial.h5part $SAVETIME $STOPTIME $RUNARGS"
ssconf.delcmd  = "qdel $JOBID"
ssconf.nocpus  = 4
ssconf.runargs = "dump 8"

axnorm = [-1.6e10, 1.6e10, -0.9e10, 0.9e10]
axzoom = [-3.2e9 , 3.2e9 , -1.8e9 , 1.8e9 ]

eVinK = 11604.505

cfg_mat_XY_n = GIplotConfig()
cfg_mat_XY_n.ax = axnorm
cfg_mat_XY_n.postPlot  = plotSimParams
cfg_mat_XY_n.imgdir = ssconf.dir + "/viz_mat_XY_n"

cfg_mat_YZ_n = GIplotConfig()
cfg_mat_YZ_n.ax = axnorm
cfg_mat_YZ_n.X = 1
cfg_mat_YZ_n.Y = 2
cfg_mat_YZ_n.Z = 0
cfg_mat_YZ_n.postPlot  = plotSimParams
cfg_mat_YZ_n.imgdir = ssconf.dir + "/viz_mat_YZ_n"


cfg_phs_XY_n = GIplotConfig()
cfg_phs_XY_n.ax = axnorm
cfg_phs_XY_n.colorFunc = colorBodyAndPhaseANEOS
cfg_phs_XY_n.postPlot  = plotSimParams
cfg_phs_XY_n.imgdir = ssconf.dir + "/viz_phs_XY_n"

cfg_phs_XY_z = GIplotConfig()
cfg_phs_XY_z.ax = axzoom
cfg_phs_XY_z.colorFunc = colorBodyAndPhaseANEOS
cfg_phs_XY_z.postPlot  = plotSimParams
cfg_phs_XY_z.imgdir = ssconf.dir + "/viz_phs_XY_z"


cfg_cid_XY_n = GIplotConfig()
cfg_cid_XY_n.ax = axnorm
cfg_cid_XY_n.colorFunc = colorClumps
cfg_cid_XY_n.postPlot  = plotSimParams
cfg_cid_XY_n.imgdir = ssconf.dir + "/viz_cid_XY_n"

cfg_cid_YZ_n = GIplotConfig()
cfg_cid_YZ_n.ax = axnorm
cfg_cid_YZ_n.X = 1
cfg_cid_YZ_n.Y = 2
cfg_cid_YZ_n.Z = 0
cfg_cid_YZ_n.colorFunc = colorClumps
cfg_cid_YZ_n.postPlot  = plotSimParams
cfg_cid_YZ_n.imgdir = ssconf.dir + "/viz_cid_YZ_n"


cfg_T___XY_n = GIplotConfig()
cfg_T___XY_n.ax        = axnorm
cfg_T___XY_n.scal_min  = 1000. / eVinK
cfg_T___XY_n.scal_max  = 9000. / eVinK
cfg_T___XY_n.scal_name = "T"
cfg_T___XY_n.cbar_plot = True
cfg_T___XY_n.cbar_sc   = eVinK
cfg_T___XY_n.cbar_ut   = "K"
cfg_T___XY_n.cbar_fmt  = '%5.0f'
cfg_T___XY_n.colorFunc = colorScalarLog10
cfg_T___XY_n.postPlot  = plotSimParams
cfg_T___XY_n.imgdir    = ssconf.dir + "/viz_T___XY_n"
cfg_T___XY_n.verbose   = True

cfg_T___XY_z = GIplotConfig()
cfg_T___XY_z.ax        = axzoom
cfg_T___XY_z.scal_min  = 1000. / eVinK
cfg_T___XY_z.scal_max  = 9000. / eVinK
cfg_T___XY_z.scal_name = "T"
cfg_T___XY_z.cbar_plot = True
cfg_T___XY_z.cbar_sc   = eVinK
cfg_T___XY_z.cbar_ut   = "K"
cfg_T___XY_z.cbar_fmt  = '%5.0f'
cfg_T___XY_z.colorFunc = colorScalarLog10
cfg_T___XY_z.postPlot  = plotSimParams
cfg_T___XY_z.imgdir    = ssconf.dir + "/viz_T___XY_z"


cfg_rho_XY_n = GIplotConfig()
cfg_rho_XY_n.ax        = axnorm
cfg_rho_XY_n.cmap      = "RdYlGn"
cfg_rho_XY_n.scal_min  = 0.
cfg_rho_XY_n.scal_max  = 6.
cfg_rho_XY_n.scal_name = "rho"
cfg_rho_XY_n.cbar_plot = True
cfg_rho_XY_n.cbar_sc   = 1.
cfg_rho_XY_n.cbar_ut   = "g/cm^3"
cfg_rho_XY_n.cbar_fmt  = '%5.1f'
cfg_rho_XY_n.colorFunc = colorScalarLinear
cfg_rho_XY_n.postPlot  = plotSimParams
cfg_rho_XY_n.imgdir    = ssconf.dir + "/viz_rho_XY_n"

cfg_rho_XY_z = GIplotConfig()
cfg_rho_XY_z.ax        = axzoom
cfg_rho_XY_z.cmap      = "RdYlGn"
cfg_rho_XY_z.scal_min  = 0.
cfg_rho_XY_z.scal_max  = 6.
cfg_rho_XY_z.scal_name = "rho"
cfg_rho_XY_z.cbar_plot = True
cfg_rho_XY_z.cbar_sc   = 1.
cfg_rho_XY_z.cbar_ut   = "g/cm^3"
cfg_rho_XY_z.cbar_fmt  = '%5.1f'
cfg_rho_XY_z.colorFunc = colorScalarLinear
cfg_rho_XY_z.postPlot  = plotSimParams
cfg_rho_XY_z.imgdir    = ssconf.dir + "/viz_rho_XY_z"


cfg_p___XY_n = GIplotConfig()
cfg_p___XY_n.ax        = axnorm
cfg_p___XY_n.cmap      = "RdYlGn"
cfg_p___XY_n.scal_min  = 1.e9
cfg_p___XY_n.scal_max  = 1.e13
cfg_p___XY_n.scal_name = "p"
cfg_p___XY_n.cbar_plot = True
cfg_p___XY_n.cbar_sc   = 1.e-10
cfg_p___XY_n.cbar_ut   = "GPa"
cfg_p___XY_n.cbar_fmt  = '%4.1f'
cfg_p___XY_n.cbar_noticks = 5
cfg_p___XY_n.colorFunc = colorScalarLog10
cfg_p___XY_n.postPlot  = plotSimParams
cfg_p___XY_n.imgdir    = ssconf.dir + "/viz_P___XY_n"

cfg_p___XY_z = GIplotConfig()
cfg_p___XY_z.ax        = axzoom
cfg_p___XY_z.cmap      = "RdYlGn"
cfg_p___XY_z.scal_min  = 1.e9
cfg_p___XY_z.scal_max  = 1.e13
cfg_p___XY_z.scal_name = "p"
cfg_p___XY_z.cbar_plot = True
cfg_p___XY_z.cbar_sc   = 1.e-10
cfg_p___XY_z.cbar_ut   = "GPa"
cfg_p___XY_z.cbar_fmt  = '%4.1f'
cfg_p___XY_z.cbar_noticks = 5
cfg_p___XY_z.colorFunc = colorScalarLog10
cfg_p___XY_z.postPlot  = plotSimParams
cfg_p___XY_z.imgdir    = ssconf.dir + "/viz_P___XY_n"

vizcfgs = [ cfg_mat_XY_n, cfg_mat_YZ_n,
    cfg_cid_XY_n, cfg_cid_YZ_n,
    cfg_phs_XY_n, cfg_phs_XY_z,
    cfg_T___XY_n, cfg_T___XY_z,
    cfg_rho_XY_n, cfg_rho_XY_z,
    cfg_p___XY_n, cfg_p___XY_z ]

simadm = SimAdmin(ssconf)
sims = simadm._sims

simA = sims["mtar001.000_mimp001.000_impa85.0_vimp1.15"]
simB = sims["mtar000.100_mimp000.100_impa00.0_vimp1.25"]
simC = sims["mtar001.000_mimp000.100_impa45.0_vimp1.50"]

simsA = simadm.getSims( SimParam(1.0, 1.0, nan, 1.15) )

vizcfg = GIvizConfig()
vizcfg.subcmd = "qsub -M andreas.reufer@space.unibe.ch -b n -cwd -l h_cpu=4:00:00 -l h_vmem=512M -N $JOBNAME $JOBCMD"
vizcfg.scdir = ssconf.dir + "/vizscratch"

vizsim = GIviz(vizcfg)

def vizAll():
  for sim in sims.values():
    if sim.nodumps > 0:
      vizsim.gatherVizTasks(sim, vizcfgs)
      vizsim.runTasksSGE()

