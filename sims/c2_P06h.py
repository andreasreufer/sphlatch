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
eVinK = 11604.505

sscfg = SimSetConfig()
sscfg.name = "c2_P06h"

sscfg.maketarg = "simple_sph_GCSHUmD_"
sscfg.binary   = "simple_sph_GCSHUmD_"
sscfg.attr     = [ ("umin", 1.e9), ("rhominclump", 1.5000), ("diskrmax",1.5e10), ("courant",0.2) ]

sscfg.nocpus  = 4
sscfg.runargs = "dump 6"

simadm = SimAdmin(sscfg)
sims = simadm._sims

cfg_orb_XY_m.imgdir = sscfg.dir + "/viz_orb_XY_m"

cfg_ecc_m.imgdir = sscfg.dir + "/viz_ecc____m"
cfg_eng_m.imgdir = sscfg.dir + "/viz_eng____m"

cfg_T___XY_n.imgdir = sscfg.dir + "/viz_T___XY_n"
cfg_T___XY_n.scal_min  =  150. / eVinK
cfg_T___XY_n.MminPlt = 1.e-5*ME
cfg_T___XY_n.MminPlt = 5.e-4*ME
cfg_T___XY_n.axreluse = True

cfg_mat_XY_n.imgdir = sscfg.dir + "/viz_mat_XY_n"
cfg_mat_XY_n.MminPlt = 1.e-5*ME
cfg_mat_XY_n.MminPlt = 5.e-4*ME
cfg_mat_XY_n.axreluse = True

cfg_grd_rho_XY_m.imgdir = sscfg.dir + "/viz_rho_XY_m"
cfg_grd_p___XY_m.imgdir = sscfg.dir + "/viz_p___XY_m"
cfg_fat_XY_Z.imgdir = sscfg.dir + "/viz_fat_XY_Z"

#mimpa = [0.1, 0.2, 0.5, 1.0, 2.0]
#mtara = [0.1, 0.2, 0.5, 1.0, 2.0]
mimpa = [0.007]
mtara = [0.01]

vimpimpaa = [ \
 #(15., 1.0), (15., 1.4), (15., 1.7), (15., 2.0), (15., 2.5), (15., 3.0), 
 #(15., 3.5), (15., 4.0), 
 #(22.5, 1.0), (22.5, 1.4), (22.5, 1.6), (22.5, 1.7), (22.5, 1.8), 
 #(22.5, 2.0), (22.5, 2.5), (22.5, 3.0), (22.5, 3.5), (22.5, 4.0), 
 #(30.0, 1.0), (30.0, 1.2), (30.0, 1.4), (30.0, 1.5), (30.0, 1.6), 
 (30.0, 2.0), (30.0, 2.5), (30.0, 3.0), (30.0, 3.5), (30.0, 4.0), 
 (37.5, 1.0), (37.5, 1.1), (37.5, 1.2), (37.5, 1.3), (37.5, 1.4), (37.5, 1.5),
 (37.5, 1.6), (37.5, 2.0), (37.5, 2.5), (37.5, 3.0), (37.5, 3.5), (37.5, 4.0), 
 (45.0, 1.0), (45.0, 1.05), (45.0, 1.1), (45.0, 1.15), (45.0, 1.2), (45.0, 1.3),
 (45.0, 1.4), (45.0, 1.7) , (45.0, 2.0), (45.0, 2.5) , (45.0, 3.0), (45.0, 3.5),
 (45.0, 4.0),
 (52.5, 1.0), (52.5, 1.05), (52.5, 1.1), (52.5, 1.15), (52.5, 1.2), (52.5, 1.3),
 (60.0, 1.0), (60.0, 1.05), (60.0, 1.1), (60.0, 1.15), (60.0, 1.2), (60.0, 1.3),
 (60.0, 1.4), (60.0, 1.7) , (60.0, 2.0), (60.0, 2.5) , (60.0, 3.0), (60.0, 3.5),
 (60.0, 4.0),
 (75.0, 1.0), (75.0, 1.05), (75.0, 1.1), (75.0, 1.15), (75.0, 1.2), (75.0, 1.3),
 (75.0, 1.4), (75.0, 1.7) , (75.0, 2.0), (75.0, 2.5) , (75.0, 3.0), (75.0, 3.5),
 (75.0, 4.0),
 (82.5, 1.0), (82.5, 1.05), (82.5, 1.1), (82.5, 1.15), (82.5, 1.2), (82.5, 1.3),
 (89.0, 1.0), (89.5, 1.05), (89.5, 1.1), (89.5, 1.15), (89.5, 1.2), (89.5, 1.3)
 ]

vimpimpaa2 = [ \
 (15., 1.0), (30., 1.0), (45., 1.0), (60., 1.0), (75., 1.0), (89.5, 1.0),
 (15., 1.05), (30., 1.05), (45., 1.05), (60., 1.05), (75., 1.05), (89.5, 1.05),
 (15., 1.1), (30., 1.1), (45., 1.1), (60., 1.1), (75., 1.1), (89.5, 1.1),
 (15., 1.2), (30., 1.2), (45., 1.2), (60., 1.2), (75., 1.2), (89.5, 1.2) ]
 
cond  = "mimp <= mtar"


vizcfg = GIvizConfig()
vizcfg.scdir = sscfg.dir + "/vizscratch"
vizcfg.tasksperjob = 100

vizsim = GIviz(vizcfg)


def vizAll():
  vizsim.clearScratch()
  vizsim.doSims(sims.values(), [cfg_T___XY_n, cfg_mat_XY_n])

def animAll():
  vizsim.animSims(sims.values(), [cfg_T___XY_n, cfg_mat_XY_n])

def findClumps():
  scrname = sscfg.name+"_do_clumps.sh"
  scrfile = open(scrname,"w")
  print >>scrfile, "#!/bin/bash"
  for sim in simadm.getSimsByState("finished"):
     print >>scrfile, sim._redoClumpsIfnodata()
  scrfile.close()
  os.chmod(scrname, stat.S_IRWXU)


def finalizeSims(sims):
  for sim in sims:
    vizsim.clearSim(sim, [cfg_T___XY_n, cfg_mat_XY_n])
  
  archname = sscfg.name + "_" + socket.gethostname() \
             + "_" + str(datetime.date.today())  + ".tar"
  simadm.archiveSims(archname, sims)

  for sim in sims:
    if not sim.state == "run":
      sim._delAll()

