#!/usr/bin/env ipython

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

from sim_admin import SimAdmin, SimParam, SimSetConfig
from gi_plot import *
from gi_viz  import GIviz, GIvizConfig
from clumps import ClumpsPlotConfig
import os
import stat
import numpy as np
import commands

nan = float('nan')
lem = 3.40000e41

sscfg = SimSetConfig()
sscfg.name = "m0"
sscfg.bodtol = 0.2
sscfg.T = 3.50e11
sscfg.tcolanim = 0.1

sscfg.attr.append( ("umin",1.e9) )
sscfg.auxfiles = "./moon_colls/aux/aneos.input ./moon_colls/aux/aneos_tables.hdf5"
sscfg.srcdir   = "/home0/areufer/repos/sphlatch/apps/simple_sph"

sscfg.maketarg = "simple_sph_GHUaE_"
sscfg.binary   = "simple_sph_GHUaE_"
sscfg.attr     = [ ("rmaxunbound",1.e10), ("rmaxbound",2.e10), ("umin", 1.e9), ("rhominclump", 1.9875) ]

#sscfg.subcmd  = "qsub -M andreas.reufer@space.unibe.ch -b y -cwd -l h_cpu=240:00:00 -l h_vmem=512M -pe smp $NOCPUS -R y -N $SIMNAME ./$BINARY initial.h5part $SAVETIME $STOPTIME $RUNARGS"
sscfg.subcmd  = "qsubPSMP $NOCPUS $SIMNAME ./$BINARY initial.h5part $SAVETIME $STOPTIME $RUNARGS"
sscfg.delcmd  = "qdel $JOBID"
sscfg.nocpus  = 8
sscfg.runargs = "dump 16"
sscfg.keepfiles = "clumps.h5part disk.hdf5 moons.pdf"

eVinK = 11604.505

#axnorm = [-8.0e9 , 8.0e9 , -4.5e9 , 4.5e9 ]
axwide = [-3.2e10, 1.6e10, -1.8e10, 0.9e10]
axmedi = [-1.6e10, 8.0e9 , -9.0e9 , 4.5e9 ]

simadm = SimAdmin(sscfg)
sims = simadm._sims


cfg_orb_XY_m.imgdir = sscfg.dir + "/viz_orb_XY_m"
cfg_ecc_m.imgdir = sscfg.dir + "/viz_ecc____m"

cfg_orb_XY_m.imgdir = sscfg.dir + "/viz_orb_XY_m"

cfg_orb_XY_n = GIplotConfig()
cfg_orb_XY_n.ax        = [-1.2e10, 0.4e10, -0.6e10, 0.3e10]
cfg_orb_XY_n.postPlot  = [plotSimParams, plotDiskMass]
cfg_orb_XY_n.colorFunc = colorBodyAndOrbAndMat
cfg_orb_XY_n.imgdir = sscfg.dir + "/viz2_orb_XY_n"
cfg_T___XY_n.imgdir = sscfg.dir + "/viz_T___XY_n"
cfg_T___XY_n.ax        = [-1.2e10, 0.4e10, -0.6e10, 0.3e10]


vizcfg = GIvizConfig()
vizcfg.subcmd = "qsubSALL $JOBNAME $JOBCMD"
vizcfg.scdir = sscfg.dir + "/vizscratch"
vizcfg.tasksperjob = 10


vizsim = GIviz(vizcfg)

def fillViz(cfgs):
  for sim in sims.values():
    if sim.nodumps > 0:
      vizsim.gatherVizTasks(sim, cfgs)
      vizsim.runTasksSGE()

def animViz(cfgs):
  for sim in sims.values():
    if sim.nodumps > 0:
      vizsim.animSim( sim, cfgs)

cfg_orb_XY_Z = GIplotConfig()
cfg_orb_XY_Z.ax = [-2.00e9 , 2.00e9 , -1.50e9 , 0.75e9 ]
cfg_orb_XY_Z.postPlot  = [plotSimParams]
cfg_orb_XY_Z.colorFunc = colorBodyAndMat
cfg_orb_XY_Z.imgdir = sscfg.dir + "/viz_orb_XY_Z"
cfg_orb_XY_Z.filt = "negzonly"


cfg_fat_XY_Z = GIplotConfig()
cfg_fat_XY_Z.ax = [-2.00e9 , 2.00e9 , -1.50e9 , 0.75e9 ]
cfg_fat_XY_Z.postPlot  = [plotSimParams]
cfg_fat_XY_Z.colorFunc = colorBodyAndMat
cfg_fat_XY_Z.imgdir = sscfg.dir + "/viz_fat_XY_Z"
cfg_fat_XY_Z.filt = "aux"
cfg_fat_XY_Z.filtFunc = filtLaterOrb
cfg_fat_XY_Z.orbtype = 2


def impaPicAll(trel):
  for sim in sims.values():
    impaPic(sim, trel)

def impaPic(sim, trel):
  tcol = sim.tcol
  for drec in sim.dumps:
    (dfile, dtime) = drec
    if abs( dtime / tcol  - trel ) < 1.e-3:
      vizsim.addTask(sim, drec, cfg_orb_XY_Z)
  vizsim.runTasksSGE()


def fatePic(sim, trel):
  (lfile, ltime) = sim.dumps[-1]
  labsfile = sim.dir + lfile

  tcol = sim.tcol
  for drec in sim.dumps:
    (dfile, dtime) = drec
    if abs( dtime / tcol  - trel ) < 1.e-3:
      (lfile, ltime) = sim.dumps[-1]
      labsfile = sim.dir + lfile
      
      cfg_fat_XY_Z.lfile = labsfile
      vizsim.addTask(sim, drec, cfg_fat_XY_Z)
  vizsim.runTasksSGE()



# clump_finder

