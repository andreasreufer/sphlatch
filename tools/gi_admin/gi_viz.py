import os
import os.path as path
import stat
import shelve
import commands
import sys
import shutil
import tempfile

from gi_plot import GIplotConfig, GIplot

class GIvizConfig(object):
  def __init__(self):
    self.dir = ""
    self.scdir = "./scdir"
    self.tasksperjob = 50
    self.subcmd = "qsub -cwd -l qname=all.q -N $JOBNAME -b n $JOBCMD"

class GIvizTask(object):
  def __init__(self, pfile, cfile, ifile, plotcfg):
    self.pfile = pfile
    self.cfile = cfile
    self.ifile = ifile

    self.plotcfg = plotcfg

  def execute(self):
    gipl = GIplot(plotcfg)
    gipl.plotParticles(pname, cname, ifile)


class GIviz(object):
  def __init__(self, cfg, plotcfg):
    self.cfg = cfg
    
    self.dir = os.path.abspath( self.cfg.dir ) + "/"
    if not path.exists(self.dir):
      os.mkdir(self.dir)

    self.scdir = os.path.abspath( self.cfg.scdir ) + "/"
    if not path.exists(self.scdir):
      os.mkdir(self.scdir)
    
  def vizSim(self, sim):
    scdir = self.scdir
    simkey = sim.params.key
    
    tasks = {}
    self.addPlotTasks(tasks, sim)
    self.plotJobScript(tasks, simkey)

  def animSim(self, sim):
    idir  = self.dir + sim.params.key + "/"
    iext  = self.plotcfg.imgext
    animext = ".mp4"
    adir = self.dir + "videos/"
    afile = adir + sim.params.key + animext
    ffmpegargs = '-b 20000k'

    def isimg(str):
      return os.path.splitext(str)[1] == iext

    imglist = filter( isimg, os.listdir( idir ) )
    imglist.sort()
      
    tmpdir = tempfile.mkdtemp( dir=idir )

    count = 0
    for img in imglist:
      os.symlink( idir + img, tmpdir + "/img" + '%06d' % count + iext )
      count += 1

    if not path.exists(adir):
      os.mkdir(adir)

    if path.exists(afile):
      os.remove(afile)

    cmd =  "ffmpeg -i " + tmpdir + "/img%06d" + iext + " " +\
        ffmpegargs + " " + afile
    print cmd
    (exstat, out) = commands.getstatusoutput(cmd)

    count = 0
    for img in imglist:
      os.remove(tmpdir + "/img" + '%06d' % count + iext)
      count += 1

    os.removedirs(tmpdir)


  def addPlotTasks(self, tasks, sim):
    dumps = sim.dumps
    for dumprec in dumps:
      (dfile, dtime) = dumprec
      dprfx = dfile.replace('.h5part','')
      pfile = sim.dir + dfile
      cfile = "clumps.h5part"

      key = sim.params.key + "_" + dprfx

      idir  = self.dir + sim.params.key + "/"
      ifile = idir + dprfx + self.plotcfg.imgext

      if not path.exists(idir):
        os.mkdir(idir)

      if path.exists(pfile) and path.exists(cfile) and not path.exists(ifile):
        open(ifile,'w').close()
        tasks[key] = GIvizTask(pfile, cfile, ifile, self.plotcfg)
        print "new plot task",pfile, cfile, ifile
  
  def plotJobScript(self, tasks, jobname):
    ascdir = self.scdir
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

    drvno   = 0
    taskno  = 0
    drvname = ascdir + jobname + "_driver" + ("%04d" % drvno ) + ".sh"
    drvhead = "#!/bin/bash \n"
    drvfoot = "rm " + drvname + "\n"
    drvstr  = drvhead

    taskdb = shelve.open(sdbname)
    taskdb.clear()
    for tkey in tasks.keys():
      taskdb[tkey] = tasks[tkey]

      drvstr += "python " + excname + " " + tkey + "\n"
      taskno += 1
      if taskno > self.cfg.tasksperjob:
        taskno = 0
        drvstr += drvfoot
    
        self.submitJob(drvstr, drvname, jobname)

        drvstr = drvhead
        drvno += 1
        drvname = ascdir + jobname + "_driver" + ("%04d" % drvno ) + ".sh"
    taskdb.close()
    drvstr += drvfoot
        
    self.submitJob(drvstr, drvname, jobname)
    
    drvstr += "rm " + drvname + "\n"

  def submitJob(self, scriptstr, scriptname, jobname):
    script = open(scriptname,"w")
    print >>script, scriptstr
    script.close()
    os.chmod(scriptname, stat.S_IRWXU)

    jobsubstr = self.cfg["JOBSUBMIT"]

    jobsubstr = jobsubstr.replace("$JOBNAME", "viz_" + jobname)
    jobsubstr = jobsubstr.replace("$JOBCMD",  scriptname)
    
    print jobsubstr
    oldwd = os.getcwd()
    #os.chdir(self.scdir)
    #(exstat, out) = commands.getstatusoutput(jobsubstr)
    #os.chdir(oldwd)


  def clearScratch(self):
    filelist = os.listdir( self.scdir )
    filelist.sort()

    for file in filelist:
      os.remove(self.scdir + file)



