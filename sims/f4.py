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
eVinK = 11604.505

sscfg = SimSetConfig()
sscfg.name = "f4"

sscfg.T = 3.50e11
sscfg.tcolanim = 0.1

sscfg.attr.append( ("umin",1.e9) )

sscfg.maketarg = "simple_sph_GCSHUmD_"
sscfg.binary   = "simple_sph_GCSHUmD_"
sscfg.attr     = [ ("umin", 1.e9), ("rhominclump", 0.7000), ("diskrmax",1.5e9), ("courant",0.2) ]

sscfg.nocpus  = 4
sscfg.runargs = "dump 6"
sscfg.keepfiles = "clumps.h5part disk.hdf5 moons.pdf"


simadm = SimAdmin(sscfg)
sims = simadm._sims

#mimpa = [0.2, 0.02, 0.002, 0.0002]
#mtara = [1.0, 0.10, 0.010, 0.0010]

mimpa = [00.2]
mtara = [01.0]

vimpimpaa = [ (15.0, 2.5), (15.0, 3.0), (15.0, 3.5), (15.0, 4.0), 
 (30.0, 2.0), (30.0, 2.5), (30.0, 3.0), (30.0, 3.5), (30.0, 4.0), 
 (45.0, 1.4), (45.0, 1.7) , (45.0, 2.0), (45.0, 2.5) , (45.0, 3.0), (45.0, 3.5),
 (45.0, 4.0),
 (60.0, 1.4), (60.0, 1.7) , (60.0, 2.0), (60.0, 2.5) , (60.0, 3.0), (60.0, 3.5),
 (60.0, 4.0) ]

#simadm.newSimMatrixVimpImpaaComb(mimpa, mtara, vimpimpaa, cond="True")

vizcfg = GIvizConfig()
vizcfg.scdir = sscfg.dir + "/vizscratch"
vizcfg.tasksperjob = 50

cfg_mat_XY_n.imgdir = sscfg.dir + "/viz_mat_XY_n"
cfg_mat_XY_n.ax = [-8.0e9, 8.0e9, -4.5e9, 4.5e9]


vizsim = GIviz(vizcfg)

def vizAll():
  vizsim.doSims(sims.values(), [cfg_mat_XY_n])

def animAll():
  vizsim.animSims(sims.values(), [cfg_mat_XY_n])


