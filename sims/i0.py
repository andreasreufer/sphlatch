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
sscfg.name = "i0"
sscfg.dir = "./" + sscfg.name
sscfg.bodiesdb = "./icy_colls/bodies"
sscfg.bodtol = 0.5
sscfg.T = 3.50e11
sscfg.tcolmax = 100.

sscfg.attr.append( ("umin",1.e9) )
sscfg.auxfiles = "./icy_colls/aux/aneos.input ./icy_colls/aux/aneos_tables.hdf5"
sscfg.srcdir   = "/home0/areufer/repos/sphlatch/apps/simple_sph"

sscfg.maketarg = "simple_sph_GHUaE_"
sscfg.binary   = "simple_sph_GHUaE_"
#sscfg.maketarg = "simple_sph_GHUAC_"
#sscfg.binary   = "simple_sph_GHUAC_"

#sscfg.subcmd  = "qsub -M andreas.reufer@space.unibe.ch -b y -cwd -l h_cpu=240:00:00 -l h_vmem=512M -pe smp $NOCPUS -R y -N $SIMNAME ./$BINARY initial.h5part $SAVETIME $STOPTIME $RUNARGS"
sscfg.subcmd  = "qsubPSMP $NOCPUS $SIMNAME ./$BINARY initial.h5part $SAVETIME $STOPTIME $RUNARGS"
sscfg.delcmd  = "qdel $JOBID"
sscfg.nocpus  = 8
sscfg.runargs = "dump 16"
sscfg.keepfiles = "clumps.h5part disk.hdf5 moons.pdf"

eVinK = 11604.505

axnorm = [-8.0e9 , 8.0e9 , -4.5e9 , 4.5e9 ]
axwide = [-3.2e10, 1.6e10, -1.8e10, 0.9e10]
axmedi = [-1.6e10, 8.0e9 , -9.0e9 , 4.5e9 ]

cfg_mat_XY_n = GIplotConfig()
cfg_mat_XY_n.ax = axnorm
cfg_mat_XY_n.postPlot  = plotSimParams
cfg_mat_XY_n.imgdir = sscfg.dir + "/viz_mat_XY_n"

cfg_orb_XY_n = GIplotConfig()
cfg_orb_XY_n.ax = [-1.2e10, 0.4e10, -0.6e10, 0.3e10]
cfg_orb_XY_n.postPlot  = [plotSimParams,  plotDiskMass]
cfg_orb_XY_n.colorFunc = colorBodyAndOrbAndMat
cfg_orb_XY_n.imgdir = sscfg.dir + "/viz_orb_XY_n"

cfg_orb_XY_w = GIplotConfig()
cfg_orb_XY_w.ax = axwide
cfg_orb_XY_w.postPlot  = [plotSimParams]
cfg_orb_XY_w.colorFunc = colorBodyAndOrbAndMat
cfg_orb_XY_w.imgdir = sscfg.dir + "/viz_orb_XY_w"

cfg_orb_XY_m.imgdir = sscfg.dir + "/viz2_orb_XY_m"

cfg_ecc_m.imgdir = sscfg.dir + "/viz_ecc____m"

cfg_orb_XY_Z = GIplotConfig()
cfg_orb_XY_Z.ax = [-2.00e9 , 2.00e9 , -1.50e9 , 0.75e9 ]
cfg_orb_XY_Z.postPlot  = [plotSimParams]
cfg_orb_XY_Z.colorFunc = colorBodyAndMat
cfg_orb_XY_Z.imgdir = sscfg.dir + "/viz_orb_XY_Z"
cfg_orb_XY_Z.filt = "negzonly"

vizcfg = GIvizConfig()
vizcfg.subcmd = "qsubSALL $JOBNAME $JOBCMD"
vizcfg.scdir = sscfg.dir + "/vizscratch"
vizcfg.tasksperjob = 10

vizsim = GIviz(vizcfg)

simadm = SimAdmin(sscfg)
sims = simadm._sims

def fillViz(cfgs):
  for sim in sims.values():
    if sim.nodumps > 0:
      vizsim.gatherVizTasks(sim, cfgs)
      vizsim.runTasksSGE()

def animViz(cfgs):
  for sim in sims.values():
    if sim.nodumps > 0:
      vizsim.animSim( sim, cfgs)

def impaPic(sim, trel):
  tcol = sim.tcol
  for drec in sim.dumps:
    (dfile, dtime) = drec
    if abs( dtime / tcol  - trel ) < 1.e-3:
      vizsim.addTask(sim, drec, cfg_orb_XY_Z)
  vizsim.runTasksSGE()

def impaPicAll(trel):
  for sim in sims.values():
    impaPic(sim, trel)

def redoClumps(sim):
  mclmp = 0.01*sim.mtot
  morbt = 0.10*sim.mtot

  cfile   = sim.dir + "clumps.h5part"
  dskfile = sim.dir + "disk.hdf5"

  print ("rm " + cfile)
  (exstat, out) = commands.getstatusoutput("rm " + cfile)
  jobname = sim.dir + "redoclumps.sh"
  scriptfile = open(jobname,"w")
  print >>scriptfile, "#!/bin/bash"
  
  for dump in sim.dumps:
    (dfilerel, dtime) = dump
    dfile = sim.dir + dfilerel
    cmd = "h5part_writeattr -i " + dfile + " -k rhominclump -v 1.9875 "
    print >>scriptfile, cmd
    
    cmd = "find_clumps_A_ " + dfile + " " + cfile + " " + str(mclmp) + " " + str(morbt) + " 1"
    print >>scriptfile, cmd
    
    cmd = "moon_pp " + dfile + " " + dskfile + " 1 0. 2.e10 100 0.1 1.5 1.e23"
    print >>scriptfile, cmd
  
  print >>scriptfile, "rm " + jobname
  scriptfile.close()

  os.chmod(jobname, stat.S_IRWXU)
  print jobname
  jobsubstr = vizcfg.subcmd
  jobsubstr = jobsubstr.replace("$JOBNAME", "redoclumps_"+sim.params.key)
  jobsubstr = jobsubstr.replace("$JOBCMD",  "redoclumps.sh")
  
  oldwd = os.getcwd()
  os.chdir(sim.dir)
  (exstat, out) = commands.getstatusoutput(jobsubstr)
  print jobsubstr
  os.chdir(oldwd)
  print out


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

