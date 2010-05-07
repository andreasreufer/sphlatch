import os
import os.path as path
import stat
import shelve

from gi_plot import GIplotConfig, GIplot

class GIvizTask(object):
  def __init__(self, pfile, cfile, ifile, plotcfg, ax, tscl):
    self.pfile = pfile
    self.cfile = cfile
    self.ifile = ifile

    self.plotcfg = plotcfg
    self.ax      = ax
    self.tscl    = tscl

  def execute(self):
    gipl = GIplot(self.pfile, self.cfile)
    gipl.plotParticles()
    gipl.plotClumps(self.tscl)
    gipl.storePlot(self.ifile, self.ax)


class GIviz(object):
  def __init__(self, cfgfile, plotcfg):
    self.cfg = {}
    if path.exists(cfgfile):
      execfile(cfgfile, self.cfg)
      print cfgfile + " loaded"
    
    self.dir = self.cfg["VIZDIR"] + "/"
    if not path.exists(self.dir):
      os.mkdir(self.dir)

    self.scdir = self.cfg["VIZSCRATCHDIR"] + "/"
    if not path.exists(self.scdir):
      os.mkdir(self.scdir)

    self.plotcfg = plotcfg

  def vizSim(self, sim, ax=[0., 0., 1., 1.]):
    scdir = self.scdir
    simkey = sim.params.key
    
    tasks = {}
    self.addPlotTasks(tasks, sim, ax)
    self.plotJobScript(tasks, simkey)


  def addPlotTasks(self, tasks, sim, ax=[0., 0., 1., 1.]):
    dumps = sim.dumps
    for dumprec in dumps:
      (dfile, dtime) = dumprec
      dprfx = dfile.replace('.h5part','')
      pfile = sim.dir + dfile
      cfile = sim.dir + "clumps/" + dfile

      key = sim.params.key + "_" + dprfx

      idir  = self.dir + sim.params.key + "/"
      ifile = idir + dprfx + self.plotcfg.imgext

      if not path.exists(idir):
        os.mkdir(idir)

      if path.exists(pfile) and path.exists(cfile) and not path.exists(ifile):
        open(ifile,'w').close()
        tasks[key] = GIvizTask(pfile, cfile, ifile, self.plotcfg, ax, sim.tscl)


def plotJobScript(self, tasks, jobname):
    ascdir = os.path.abspath( self.scdir ) + "/"
    drvname = ascdir + jobname + "_driver.sh"
    excname = ascdir + jobname + "_exec.py"
    sdbname = ascdir + jobname
    
    excstr = "#!/usr/bin/env python\n" +\
        "from gi_viz import GIvizTask\n" +\
        "from gi_plot import GIplot\n" +\
        "import shelve, sys\n" +\
        "taskkey = sys.argv[1]\n" +\
        "tasks = shelve.open(\"" + sdbname + "\")\n" +\
        "tasks[taskkey].execute()\n" +\
        "tasks.close()\n"

    excfile = open(excname,"w")
    print >>excfile, excstr
    excfile.close()
    os.chmod(excname, stat.S_IRWXU)

    taskdb = shelve.open(sdbname)
    taskdb.clear()
    drvstr = "#!/bin/bash \n"
    for tkey in tasks.keys():
      taskdb[tkey] = tasks[tkey]
      drvstr += "ipython " + excname + " " + tkey + "\n"
    taskdb.close()

    drvstr += "rm " + sdbname + "\n"
    drvstr += "rm " + excname + "\n"
    drvstr += "rm " + drvname + "\n"

    drvfile = open(drvname,"w")
    print >>drvfile, drvstr
    drvfile.close()
    os.chmod(drvname, stat.S_IRWXU)


