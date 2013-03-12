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
sscfg.name = "f2"
sscfg.dir = "./" + sscfg.name
sscfg.bodiesdb = "fluffy_colls/bodies"
sscfg.resultsdb = sscfg.dir + "/results"
sscfg.bodtol = 0.2
sscfg.T = 3.50e11
sscfg.tcolanim = 0.1
sscfg.relsep = 10.

sscfg.auxfiles = "./aux/aneos.input ./aux/aneos_tables.hdf5"
sscfg.srcdir   = "/home0/areufer/repos/sphlatch/apps/simple_sph"

sscfg.maketarg = "simple_sph_GHUaE_"
sscfg.binary   = "simple_sph_GHUaE_"
sscfg.attr     = [ ("rmaxunbound",5.e10), ("rmaxbound",5.e10), ("umin", 1.e9), ("rhominclump", 1.9875), ("diskrmax",5.e10) ]

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

cfg_orb_XY_m.imgdir = sscfg.dir + "/viz_orb_XY_m"
cfg_orb_XY_m.ax        = [-3.2e10, 1.6e10, -1.8e10, 0.9e10]

cfg_ecc_m.imgdir = sscfg.dir + "/viz_ecc____m"
cfg_eng_m.imgdir = sscfg.dir + "/viz_eng____m"

vizcfg = GIvizConfig()
vizcfg.subcmd = "qsubSALL $JOBNAME $JOBCMD"
vizcfg.scdir = sscfg.dir + "/vizscratch"
vizcfg.tasksperjob = 10


simadm = SimAdmin(sscfg)
sims = simadm._sims

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

cfg_grd_rho_XY_m.imgdir = sscfg.dir + "/viz_grd_rho_XY_m"
cfg_grd_p___XY_m.imgdir = sscfg.dir + "/viz_grd_p___XY_m"

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


def vizAll():
  for sim in sims.values():
    vizsim.doSims(sims.values(), [cfg_T___XY_n, cfg_grd_rho_XY_m, cfg_grd_p___XY_m])

def animAll():
  for sim in sims.values():
    vizsim.doSims(sims.values(), [cfg_T___XY_n, cfg_grd_rho_XY_m, cfg_grd_p___XY_m])


# clump_finder

