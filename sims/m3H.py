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
import numpy as np
import commands
import time

nan = float('nan')
lem = 3.40000e41

sscfg = SimSetConfig()
sscfg.name = "m3H"
#sscfg.attr     = [ ("umin", 1.e9), ("rhominclump", 2.0000), ("diskrmax",5.e10) ]
sscfg.attr     = [ ("umin", 1.e9), ("rhominclump", 2.0000), ("diskrmax",1.5e10), ("rmaxunbound", 1.5e10), ("rmaxbound", 1.5e10) ]
sscfg.nocpus  = 12
sscfg.runargs = "dump 32"
sscfg.desnop  = 1000000
sscfg.bodtol  = 0.22
sscfg.maketarg = "simple_sph_GMSHUmD_"
sscfg.binary   = "simple_sph_GMSHUmD_"

simadm = SimAdmin(sscfg)

sims = simadm._sims

eVinK = 11604.505

cfg_orb_XY_m.imgdir = sscfg.dir + "/viz_orb_XY_m"
#cfg_orb_XY_m.ax        = [-3.2e10, 1.6e10, -1.8e10, 0.9e10]
cfg_orb_XY_m.ax        = [-1.2e10, 0.4e10, -0.6e10, 0.3e10]

cfg_ecc_m.imgdir = sscfg.dir + "/viz_ecc____m"
cfg_eng_m.imgdir = sscfg.dir + "/viz_eng____m"
cfg_T___XY_n.imgdir = sscfg.dir + "/viz_T___XY_n"
cfg_mat_XY_n.imgdir = sscfg.dir + "/viz_mat_XY_n"

cfg_grd_rho_XY_m.imgdir = sscfg.dir + "/viz_rho_XY_m"
cfg_grd_p___XY_m.imgdir = sscfg.dir + "/viz_p___XY_m"
cfg_fat_XY_Z.imgdir = sscfg.dir + "/viz_fat_XY_Z"

cfg_T___XY_n.ax        = [-0.667e10, 0.4e10, -0.3e10, 0.3e10]
cfg_T___XY_n.axreluse  = False
cfg_mat_XY_n.ax        = [-0.667e10, 0.4e10, -0.3e10, 0.3e10]
cfg_mat_XY_n.axreluse  = False
cfg_mat_XY_n.postPlot  = [plotSimParams, plotDiskMass]

#cfg_mat_XY_n.ax        = [-16.0e9, 8.0e9 , -9.0e9 , 4.5e9 ]
#cfg_mat_XY_n.axreluse  = False
#cfg_mat_XY_n.postPlot  = [plotSimParams, plotDiskMass]


mimpa = [0.1, 0.2, 0.5, 1.0, 2.0]
mtara = [0.1, 0.2, 0.5, 1.0, 2.0]
impaa = [0.1, 15., 30., 45., 60., 75.]
vimpa = [1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4]
cond  = "mimp <= mtar"


vizcfg = GIvizConfig()
vizcfg.scdir = sscfg.dir + "/vizscratch"
vizcfg.tasksperjob = 1

vizsim = GIviz(vizcfg)

vizcfgs = [cfg_T___XY_n, cfg_orb_XY_m, cfg_mat_XY_n]

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
  vizsim.doSims(sims.values(), vizcfgs)

def animAll():
  vizsim.animSims(sims.values(), vizcfgs)

def housekeepingLoop():
  while True:
    print "simadm.update()"
    simadm.update()

    print "visualize ..."
    vizAll()
  
    print "sleep for half-an-hour ..."
    time.sleep(1800)
  
    print "animate ..."
    animAll()
    time.sleep(60)




